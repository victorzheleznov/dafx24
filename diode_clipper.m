% -------------------------------------------------------------------------
% This script considers application of the AA-IIR method with linear or 
% cubic interpolation and a first-order Butterworth antialiasing filter
% for simulation of the diode clipper circuit.
%
% Input signal can be either a sine tone or a sine sweep. For a sine
% tone input signal, signal-to-noise ratio (SNR) and noise-to-mask ratio 
% (NMR) are computed for an output signal to measure aliasing.
%
% Circuit equation is discretized using the trapezoidal rule. Damped
% Newton-Raphson method is chosen to solve the nonlinear equation. For
% cubic AA-IIR, the trapezoidal quadrature is used for integral 
% calculation.
%
% Author: Victor Zheleznov
% -------------------------------------------------------------------------

clear all; close all;

% add libraries
addpath("interpolation");
addpath("metrics");
addpath("utils");

%% parameters
% global parameters
fs = 44.1e3;                    % base sampling rate [Hz]

% input signal parameters
input = "sine";                 % "sine", "sweep"
t_dur = 1.5;                    % duration [sec]
amp = 10;                       % amplitude [V]
if strcmp(input, "sine") == 1
    f0 = 2.^((83-69)/12) * 440; % input frequency [Hz]
elseif strcmp(input, "sweep") == 1
    f0 = 20;                    % lowest sweep frequency [Hz]
    f1 = 22e3;                  % highest sweep frequency [Hz]
end

% interpolation parameters
interp = "cubic";               % "linear", "cubic"

% antialiasing filter parameters (first-order Butterworth filter)
fc = 0.33;                      % normalised cutoff

% numerical integration (for cubic AA-IIR)
M = 8;                          % number of quadrature points

% Newton-Raphson parameters
rel_tol = 1e-12;                % relative tolerance
abs_tol = 1e-14;                % absolute tolerance
max_iter = 50;                  % maximum number of iterations
max_subiter = 5;                % maximum number of sub-iterations

% physical parameters
res = 1e3;                      % resistance [Ohms]
cap = 33e-9;                    % capacitance [F]
Is = 2.52e-9;                   % diode saturation current [A]
Vt = 25.83e-3;                  % diode thermal voltage [V]
Ni = 1.752;                     % diode ideality factor

%% pre-processing
% check stability
assert((strcmp(interp, "linear") && fc > 0) || (strcmp(interp, "cubic") && fc >= 0.327), "Compensation filter is unstable!");

% synthesise input
if strcmp(input, "sine") == 1
    t = (0:1/fs:t_dur).';
    u = amp.*sin(2*pi*f0*t);
elseif strcmp(input, "sweep") == 1
    [t,u] = gen_sine_sweep(t_dur, fs, f0, f1, amp);
end

% define differential equation (dv/dt = Av + Bu + Cf(v))
A = -1/(res*cap);
B = 1/(res*cap);
C = -1/cap;

% define coefficient and input gain for nonlinear function sinh()
coef = 2*Is;
gain = 1/(Ni*Vt);

% precompute discretization constants for trapezoidal rule
T = 1/fs;
Hm = (2/T) - A;
Hp = (2/T) + A;
K = C/Hm;

% calculate partial fraction decomposition of the antialiasing filter
wc = 2*pi*fc*fs;
[b,a] = butter(1, wc, 's');
[r,p] = residue(b, a);

% define compensation filter transfer function
E = exp(p*T);
if strcmp(interp, "linear") == 1
    % define linearization
    b_lin = (r/p^2)*[(1/T)*(E-1) - p, (p-1/T)*E + (1/T)];
    a_lin = [1, -E];
elseif strcmp(interp, "cubic") == 1
    % define trapezoidal quadrature for integral calculation
    tq = (T/M*(0:M)).';
    wq = (T/M)*[0.5; ones(length(tq)-2,1); 0.5];

    % calculate convolution integrals for each summand of polynomial
    S = r*[sum(wq.*exp(p*(T-tq)));                                      % bias
           sum(wq.*tq.*exp(p*(T-tq)));                                  % linear
           sum(wq.*tq.^2.*exp(p*(T-tq)));                               % quadratic
           sum(wq.*tq.^3.*exp(p*(T-tq)))];                              % cubic

    % define transfer function numerators for each polynomial coefficient
    gamma = sqrt(3) - 2;
    Hc = [[ 0,              1,                         -gamma  ];       % bias 
          [-sqrt(3)*gamma,  3-sqrt(3)+2*sqrt(3)*gamma,  3*gamma]*fs;    % linear
          [-3*gamma,        6*gamma,                   -3*gamma]*fs^2;  % quadratic
          [ gamma,         -2*gamma,                    gamma  ]*fs^3]; % cubic
    denom = [1, -gamma];

    % calculate linearization transfer function
    b_lin = sum(S.*Hc, 1);
    a_lin = conv(denom, [1 -E]);
end
b_comp = a_lin ./ b_lin(1);
a_comp = b_lin ./ b_lin(1);

% initialise system state
v1 = 0;
v = 0;

% initialise the AA-IIR method
y1 = 0;                                               % output of the regular AA-IIR method
f1 = 0;                                               % output of the compensation filter
zc1 = zeros(max(length(a_comp),length(b_comp))-1, 1); % delays of the compensation filter
if strcmp(interp, "cubic") == 1
    w1 = cspline2_init(v1);                           % output of the interpolation pre-filter 
end

% initialiase output
out = zeros(size(u));      % diode clipper output
num_iter = zeros(size(u)); % number of Newton-Raphson iterations
N = length(u);

%% processing
% time loop
for n = 2:N
    % compute known quantity
    q = (Hp*v1 + B*(u(n-1) + u(n)) + C*f1) / Hm;

    % initialise Newton-Raphson method
    iter = 0;     % iteration count
    stop = false; % criteria stopping flag
    v = v1;       % initial guess

    % apply Newton-Raphson method
    while (stop == false) && (iter < max_iter)
        % calculate step
        if strcmp(interp, "linear") == 1
            [f,df] = sinh_aaiir_linear(v, v1, y1, zc1, fs, p, r, coef, gain, b_comp, a_comp);
        elseif strcmp(interp, "cubic") == 1
            [f,df] = sinh_aaiir_cubic(v, v1, y1, zc1, w1, fs, p, r, coef, gain, b_comp, a_comp, tq, wq);
        end
        G = q + K*f - v;
        DG = K*df - 1;
        step = G / DG;
        
        % apply damping
        iter_damped = 0;
        Gnext = inf;
        while (abs(Gnext) > abs(G)) && (iter_damped < max_subiter)
            if iter_damped ~= 0
                step = 0.5*step;
            end
            v_next = v - step;
            if strcmp(interp, "linear") == 1
                [f_next,~,y,zc] = sinh_aaiir_linear(v_next, v1, y1, zc1, fs, p, r, coef, gain, b_comp, a_comp);
            elseif strcmp(interp, "cubic") == 1
                [f_next,~,y,zc,w] = sinh_aaiir_cubic(v_next, v1, y1, zc1, w1, fs, p, r, coef, gain, b_comp, a_comp, tq, wq);
            end
            Gnext = q + K*f_next - v_next;
            iter_damped = iter_damped + 1;
        end

        % move to the next iteration
        v = v_next;
        iter = iter + 1;
        num_iter(n) = num_iter(n) + iter_damped;

        % check stopping scriteria
        stop = (abs(step) <= (abs_tol + rel_tol*abs(v1)));
    end

    % write output
    out(n) = v;
    
    % shift system state
    v1 = v;

    % shift the AA-IIR method
    f1 = f;
    y1 = y;
    zc1 = zc;
    if strcmp(interp, "cubic") == 1
        w1 = w;
    end
end

%% output
% display number of Newton-Raphson iterations
disp("Number of Newton-Raphson iterations:");
disp("Mean = " + mean(num_iter));
disp("Max = " + max(num_iter));
fprintf("\n");

% plot output waveform
fig_t = figure; hold on;
plot(t, out, 'k');
xlim([0 4/f0]);
xlabel("Time [sec]", 'Interpreter', 'latex');
ylabel("Voltage [V]", 'Interpreter', 'latex');
title("Diode clipper waveform for AA-IIR with " + interp + " interpolation", 'Interpreter', 'latex');
ax_spec = fig_t.CurrentAxes;
set(ax_spec.XAxis, 'TickLabelInterpreter', 'latex');
set(ax_spec.YAxis, 'TickLabelInterpreter', 'latex');

if strcmp(input, "sine") == 1
    % calculate metrics
    [snr, nmr] = metrics_sin(out, fs, f0, true, true);
    cmd_str_metr = sprintf("SNR = %.2f dB\nNMR = %.2f dB\n", snr, nmr);
    title_str_metr = sprintf("SNR = $%.2f$ dB, NMR = $%.2f$ dB", snr, nmr);
    disp("Metrics:");
    disp(cmd_str_metr);

    % plot spectrum
    [spec, fig_spec] = myfft(out, fs, "blackman");
    add_harmonic_marks(spec, fig_spec, fs, f0, true);
    title("Diode clipper spectrum for AA-IIR with " + interp + " interpolation", 'Interpreter', 'latex');
    subtitle(title_str_metr, 'Interpreter', 'latex');
    ylim([-120 10]);
elseif strcmp(input, "sweep") == 1
    % plot spectrogram
    [spec, fig_spec] = myspec(out, fs, 1024, 0.9921875, "blackman");
    title("Diode clipper spectogram for AA-IIR with " + interp + " interpolation", 'Interpreter', 'latex');
end

%% FUNCTIONS
% calculate linear AA-IIR solution for sinh() nonlinearity
% input: 
%   x --- current input sample;
%   x1 --- previous input sample;
%   y1 --- previous output of the regular AA-IIR method;
%   zc1 --- previous delays of the compensation filter;
%   fs --- sample rate [Hz];
%   p --- single real pole value;
%   r --- residue value;
%   coef --- coefficient for sinh();
%   gain --- input gain for sinh();
%   b_comp --- numerator of the compensation filter;
%   a_comp --- denominator of the compensation filter.
% output:
%   y_comp --- output of the compensation filter;
%   dy_comp --- derivative of this output with respect to the current input sample;
%   y --- output of the regular AA-IIR method;
%   zc --- delays of the compensation filter.
function [y_comp,dy_comp,y,zc] = sinh_aaiir_linear(x, x1, y1, zc1, fs, p, r, coef, gain, b_comp, a_comp)
    % calculate parameters
    T = 1/fs;
    E = exp(p*T);

    % apply input signal gain
    x = x * gain;
    x1 = x1 * gain;
    
    % calculate AA-IIR step
    I = p*T^2/((p*T)^2 - (x-x1)^2)*(E*(sinh(x1) + (x-x1)/(p*T)*cosh(x1)) - (sinh(x) + (x-x1)/(p*T)*cosh(x)));
    I = coef * I;
    y = r*I + E*y1;

    % calculate derivative of the AA-IIR output with respect to the current input sample
    dI = 2*p*T^2 * (x-x1) * E * sinh(x1);
    dI = dI + T * E * ((p*T)^2 + (x-x1)^2) * cosh(x1);
    dI = dI + T * (x-x1) * ((x-x1)^2 - 2*p*T - (p*T)^2) * sinh(x);
    dI = dI + T * ((p*T-1)*(x-x1)^2 - (p*T+1)*(p*T)^2) * cosh(x);
    dI = dI / ((p*T)^2 - (x-x1)^2)^2;                              % derivative with respect to the amplified sample
    dI = dI * gain;                                                % adjust for input gain
    dI = dI * coef;                                                % adjust for function coefficient
    dy = r*dI;                                                     % multiply by residue
    
    % apply compensation filter
    [y_comp,zc] = filter(b_comp, a_comp, y, zc1);
    dy_comp = b_comp(1)*dy;
end

% calculate cubic AA-IIR solution for sinh() nonlinearity
% input: 
%   x --- current input sample;
%   x1 --- previous input sample;
%   y1 --- previous output of the regular AA-IIR method;
%   zc1 --- previous delays of the compensation filter;
%   w1 --- previous output of the interpolation pre-filter;
%   fs --- sample rate [Hz];
%   p --- single real pole value;
%   r --- residue value;
%   coef --- coefficient for sinh();
%   gain --- input gain for sinh();
%   b_comp --- numerator of the compensation filter;
%   a_comp --- denominator of the compensation filter;
%   tq --- quadrature nodes;
%   wq --- quadrature weights.
% output:
%   y_comp --- output of the compensation filter;
%   dy_comp --- derivative of this output with respect to the current input sample;
%   y --- output of the regular AA-IIR method;
%   zc --- delays of the compensation filter;
%   w --- output of the interpolation pre-filter
function [y_comp,dy_comp,y,zc,w] = sinh_aaiir_cubic(x, x1, y1, zc1, w1, fs, p, r, coef, gain, b_comp, a_comp, tq, wq)
    % calculate parameters
    T = 1/fs;
    E = exp(p*T);

    % calculate coefficients of interpolating polynomial
    [c,w] = cspline2_step(x, x1, fs, w1);

    % define the integrated expression
    f = @(t) coef * sinh(gain * horner(t, c));
    g = @(t) f(t).*exp(p*(T-t));

    % calculate AA-IIR step
    I = sum(wq.*g(tq));
    y = r*I + E*y1;

    % calculate derivative of the AA-IIR output with respect to the current input sample
    gamma = sqrt(3) - 2;
    dc = [gamma*fs^3, -3*gamma*fs^2, -sqrt(3)*gamma*fs, 0];             % partial derivative of polynomial coefficients to the current input sample
    df = @(t) coef * gain * cosh(gain * horner(t, c)) .* horner(t, dc); % partial derivative of the nonlinear function to the current input sample
    dg = @(t) df(t).*exp(p*(T-t));                                      % partial derivative of the integrated expression to the current input sample
    dI = sum(wq.*dg(tq));
    dy = r*dI;

    % apply compensation filter
    [y_comp,zc] = filter(b_comp, a_comp, y, zc1);
    dy_comp = b_comp(1)*dy;
end
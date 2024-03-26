% -------------------------------------------------------------------------
% This script considers application of the AA-IIR method with linear or 
% cubic interpolation and a Chebyshev type 1 antialiasing filter
% for the hard clipping nonlinearity.
%
% Input signal can be either a sine tone or a sine sweep. For a sine
% tone input signal, signal-to-noise ratio (SNR) and noise-to-mask ratio 
% (NMR) are computed for an output signal to measure aliasing.
%
% For cubic AA-IIR, the midpoint quadrature is used for integral 
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
input = "sweep";                % "sine", "sweep"
t_dur = 5;                      % signal duration [sec]
amp = 10;                       % amplitude
if strcmp(input, "sine") == 1
    f0 = 2.^((83-69)/12) * 440; % input sine frequency [Hz]
elseif strcmp(input, "sweep") == 1
    f0 = 20;                    % lowest sweep frequency [Hz]
    f1 = 22e3;                  % highest sweep frequency [Hz]
end

% interpolation parameters
interp = "cubic";               % "linear", "cubic"
L = 6;                          % look-ahead for cubic interpolation

% antialiasing filter parameters (Chebyshev type 1 filter)
K = 8;                          % order
fc = 0.4;                       % normalised cutoff
Rp = 0.05;                      % passband ripple [dB]

% numerical integration (for cubic AA-IIR)
M = 8;                          % number of quadrature points

% algorithm tolerance (for linear AA-IIR)
tol = 1e-12;

% output trimming to remove boundary effects (look-ahead, filtering etc.)
cut = 100;

%% pre-processing
% synthesise input
if strcmp(input, "sine") == 1
    t = (0:1/fs:t_dur).';
    u = amp.*sin(2*pi*f0*t);
elseif strcmp(input, "sweep") == 1
    [t,u] = gen_sine_sweep(t_dur, fs, f0, f1, amp);
end

% calculate partial fraction decomposition of the antialiasing filter
wc = 2*pi*fc*fs;
[b,a] = cheby1(K, Rp, wc, 's');
[r,p,k] = residue(b, a);
assert(length(unique(p)) == length(p), "Only distinct poles are supported!");

% initialise the AA-IIR method
y = zeros(size(u));              % output of the AA-IIR method
m = 1;                           % pole counter
if strcmp(interp, "cubic") == 1
    pp = cspline(t, u, fs, L-2); % construct interpolating piecewise cubic polynomial
end

%% processing
% iterate through poles of the antialiasing filter
while m <= length(p)
    if strcmp(interp, "linear") == 1
        % process linear AA-IIR
        y = y + hard_aaiir_linear(u, fs, p(m), r(m), tol);
    elseif strcmp(interp, "cubic") == 1
        % process cubic AA-IIR
        y = y + hard_aaiir_cubic(pp, fs, p(m), r(m), M);
    end

    % increment pole index
    m = m + 1 + (imag(p(m)) ~= 0);
end

% add constant term
if ~isempty(k)
    assert(length(k) == 1, 'Direct term of degree >= 1 is not supported!')
    y = y + k*hard(u);
end

%% output
% trim the output signal
t = t(cut:end-cut);
y = y(cut:end-cut);

% plot output waveform
fig_t = figure; hold on;
plot(t, y, 'k');
xlim([t(1), t(1)+4/f0]);
xlabel("Time [sec]", 'Interpreter', 'latex');
ylabel("Amplitude", 'Interpreter', 'latex');
title("Hard clipper waveform for AA-IIR with " + interp + " interpolation", 'Interpreter', 'latex');
ax_spec = fig_t.CurrentAxes;
set(ax_spec.XAxis, 'TickLabelInterpreter', 'latex');
set(ax_spec.YAxis, 'TickLabelInterpreter', 'latex');

if strcmp(input, "sine") == 1
    % calculate metrics
    [snr, nmr] = metrics_sin(y, fs, f0, true);
    cmd_str_metr = sprintf("SNR = %.2f dB\nNMR = %.2f dB\n", snr, nmr);
    title_str_metr = sprintf("SNR = $%.2f$ dB, NMR = $%.2f$ dB", snr, nmr);
    disp(cmd_str_metr);

    % plot spectrum
    [spec, fig_spec] = myfft(y, fs, "blackman");
    add_harmonic_marks(spec, fig_spec, fs, f0, true);
    title("Hard clipper spectrum for AA-IIR with " + interp + " interpolation");
    subtitle(title_str_metr, 'Interpreter', 'latex');
    ylim([-120 10]);
elseif strcmp(input, "sweep") == 1
    % plot spectrogram
    [spec, fig_spec] = myspec(y, fs, 1024, 0.9921875, "blackman");
    title("Hard clipper spectogram for AA-IIR with " + interp + " interpolation");
end

%% FUNCTIONS
% calculate hard clipper output
% input:
%   x - input.
% output:
%   y - output.
function y = hard(x)
    y = min(max(x, -1), 1);
end

% calculate linear AA-IIR solution for hard clipper
% input: 
%   x --- input vector;
%   fs --- sample rate [Hz];
%   p --- single pole value (real or complex);
%   r --- residue value;
%   tol --- tolerance value.
% output:
%   y --- output of the AA-IIR method.
function y = hard_aaiir_linear(x, fs, p, r, tol)
    % calculate parameters
    N = length(x);
    T = 1/fs;
    E = exp(p*T);
    F = (E - 1)/p;

    % check real or complex pole
    if imag(p) ~= 0
        r = 2*r;
    end

    % time loop
    y = zeros(N,1);
    x1 = x(1);
    for n = 2:N
        % calculate convolution integral (analytical solution)
        if (abs(x(n)-x1) < tol)
            I = hard(0.5*(x(n)+x1)) * F;
        else
            if x1 < x(n)
                if x(n) <= -1
                    I = (1 - E)/p;
                elseif x1 <= -1 && x(n) > -1 && x(n) <= 1
                    I = ((x(n)-x1)/T*(exp(p*T*(x(n)+1)/(x(n)-x1)) - 1) - p*(x(n) + E))/p^2;
                elseif x1 <= -1 && x(n) > 1
                    I = ((x(n)-x1)/T*(exp(p*T*(x(n)+1)/(x(n)-x1)) - exp(p*T*(x(n)-1)/(x(n)-x1))) - p*(1+E))/p^2;
                elseif x1 > -1 && x(n) <= 1
                    I = (E*((x(n)-x1)/T + p*x1) - (x(n)-x1)/T -p*x(n))/p^2;
                elseif x1 > -1 && x1 <= 1 && x(n) > 1
                    I = ((x1-x(n))/T*exp(p*T*(x(n)-1)/(x(n)-x1)) + E*((x(n)-x1)/T + p*x1) - p)/p^2;
                else
                    I = (E - 1)/p;
                end
            else
                if x1 <= -1
                    I = (1 - E)/p;
                elseif x1 > -1 && x1 <= 1 && x(n) <= -1
                    I = ((x1-x(n))/T*exp(p*T*(x(n)+1)/(x(n)-x1)) + E*((x(n)-x1)/T + p*x1) + p)/p^2;
                elseif x1 <= 1 && x(n) > -1
                    I = (E*((x(n)-x1)/T + p*x1) + (x1-x(n))/T - p*x(n))/p^2;
                elseif x1 > 1 && x(n) <= -1
                    I = ((x(n)-x1)/T*(exp(p*T*(x(n)-1)/(x(n)-x1)) - exp(p*T*(x(n)+1)/(x(n)-x1))) + p*(1+E))/p^2;
                elseif x1 > 1 && x(n) > -1 && x(n) <= 1
                    I = (exp(p*T*(x(n)-1)/(x(n)-x1))*(x(n)-x1)/T - p*x(n) - (x(n)-x1)/T + p*E)/p^2;
                else
                    I = (E - 1)/p;
                end
            end
        end
                    
        % calculate AA-IIR step
        y(n) = r*I + E*y(n-1);

        % shift input
        x1 = x(n);
    end

    % check imaginary part
    if imag(p) ~= 0
        y = real(y);
    end
end

% calculate cubic AA-IIR solution for hard clipper
% input: 
%   pp --- piecewise cubic polynomial interpolation of input signal;
%   fs --- sample rate [Hz];
%   p --- single pole value (real or complex);
%   r --- residue value;
%   M --- number of quadrature points for numerical integration.
% output:
%   y --- output of the AA-IIR method.
function y = hard_aaiir_cubic(pp, fs, p, r, M)
    % calculate parameters
    N = length(pp.breaks);
    T = 1/fs;
    E = exp(p*T);

    % define midpoint quadrature
    tq = (T/M)*((1:M)-0.5).';
    wq = (T/M)*ones(size(tq));

    % check real or complex pole
    if imag(p) ~= 0
        r = 2*r;
    end

    % define the integrated expression
    f = @(t,idx) hard(horner(t, pp.coefs(idx-1,:)));
    g = @(t,idx) f(t,idx).*exp(p*(T-t));

    % time loop
    y = zeros(N,1);
    for n = 2:N
        % calculate convolution integral (numerical integration)
        I = sum(wq.*g(tq,n));
        
        % calculate AA-IIR step
        y(n) = r*I + E*y(n-1);
    end

    % check imaginary part
    if imag(p) ~= 0
        y = real(y);
    end
end
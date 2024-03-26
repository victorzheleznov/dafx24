% calculate bandlimited version of a nonlinearly processed sine signal
% input:
%   sig --- test signal;
%   fs --- base sample rate [Hz];
%   f0 --- fundamental frequency of input sine [Hz];
%   IS_SYM --- symmetry flag for nonlinearity;
%   APPLY_LP --- low-pass filter flag (optional);
%   PLOT_SIG --- flag to plot input and bandlimited signals (optional).
% output:
%   sig --- test signal (truncated and DC compensated => aligned with 
%           bandlimited version);
%   sig_lim --- bandlimited version of a signal;
%   alias --- aliases from a signal.
% references:
%   [1] Carson, A. (2020). "Aliasing Reduction in Virtual Analogue Modelling". 
%       thesis: https://www.researchgate.net/publication/344362244_Aliasing_Reduction_in_Virtual_Analogue_Modelling
%       code: was provided via email request.
%   [2] Kahles, J., Esqueda, F., & Valimaki, V. (2019). "Oversampling for 
%       Nonlinear Waveshaping: Choosing the Right Filters". Journal of the 
%       Audio Engineering Society.
%       https://www.aes.org/tmpFiles/elib/20230722/20485.pdf
function [sig, sig_lim, alias] = bandlimit_signal(sig, fs, f0, IS_SYM, varargin)
    % parameters
    sidelobe_att = 120; % chebyshev window sidelobe attenuation [dB]
    dc_thres = 40;      % DC detection threshold [dB]

    % check filtering flag
    if length(varargin) >= 1
        APPLY_LP = varargin{1};
    else
        APPLY_LP = false;
    end

    % check plot flag for debugging
    if length(varargin) >= 2
        PLOT_SIG = varargin{2};
    else
        PLOT_SIG = false;
    end
    
    % filter input signal [1]
    if APPLY_LP == true
        % use the same filter as in decimate() for oversampled algorithms
        [b,a] = cheby1(8, 0.05, 0.8);
        sig = filter(b, a, sig);
    end

    % truncate signal to remove transient from initial condition [1]
    % and leave one-second fragment [2]
    if length(sig) > fs
        sig = sig(end-fs+1:end);
    else
        warning('Signal is too short for robust metrics calculation!');
    end
    N = length(sig);

    % check odd length
    assert(rem(N,2) == 0, 'Signal should have an odd length after truncation!');

    % calculate harmonics [1]
    num_harmonics = floor(0.5*fs/f0);
    if IS_SYM == true
        % even harmonics for symmetric nonlinearity
        harmonics = f0 * (1:2:num_harmonics).';
    else
        % even + odd harmonics for asymmetric nonlinearity
        harmonics = f0 * (1:num_harmonics).';
    end

    % Chebyshev window with 120 dB sidelobe attenuation [2]
    win = chebwin(N, sidelobe_att);

    % calculate single-sided spectrum
    S = fft(win.*sig, N);  % calculate FFT
    S = S(1:N/2+1);        % keep only the half for the real signal

    % adjust for scalloping loss (including window)
    bins_exact = (harmonics*N/fs) + 1; % exact bins for harmonics
    bins = round(bins_exact);          % discrete bins for harmonics
    d = bins_exact - bins;             % difference between exact and discrete bins
    idx = (0:N-1).';                   % indexes for summation  
    amps = zeros(size(bins));          % array for storing adjusted FFT data
    phases = zeros(size(bins));
    for i = 1:length(bins)
        % calculate adjusted FFT value at discrete bin
        tmp = S(bins(i))/sum(win.*exp(1j*2*pi*d(i)/N.*idx));
        
        % extract magnitude and phase
        amps(i) = 2*abs(tmp);   % multiplying by 2 for a real-valued signal
        phases(i) = angle(tmp);
    end

    % synthesise bandlimited signal [1]
    t = ((0:N-1)/fs).';
    sig_lim = zeros(size(t));
    for i = 1:length(harmonics)
        sig_lim = sig_lim + amps(i)*cos(2*pi*harmonics(i)*t + phases(i));
    end
    
    % remove DC if needed
    if 20*log10(abs(S(1))./max(abs(S))) > -dc_thres
        disp("bandlimit_signal.m: DC detected!");
        DC = abs(S(1)/sum(win));
        sig = sig - DC;
    end

    % calculate aliased components [1]
    alias = sig - sig_lim;

    % debug plots
    if PLOT_SIG == true
        fig_sig = figure; hold on;
        plot(sig, 'color', "#e6194b");
        plot(sig_lim, 'color', "#3cb44b");
        plot(alias, 'color', "#4363d8");
        xlim([-100, length(sig)+100]);
        xlabel('Samples [n]', 'Interpreter', 'latex');
        ylabel('Amplitude', 'Interpreter', 'latex');
        legend({'Input signal', 'Bandlimited signal', 'Alias signal'}, 'Interpreter', 'latex');
        title(sprintf("Input signal at $f_0 = %.2f$ Hz", f0), 'Interpreter', 'latex');
        ax_spec = fig_sig.CurrentAxes;
        set(ax_spec.XAxis, 'TickLabelInterpreter', 'latex');
        set(ax_spec.YAxis, 'TickLabelInterpreter', 'latex');
    end
end
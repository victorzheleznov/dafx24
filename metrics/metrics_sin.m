% calculate SNR and NMR for a nonlinearly processed sine signal
% input:
%   sig --- test signal;
%   fs --- sample rate [Hz];
%   f0 --- fundamental frequency of input sine [Hz];
%   IS_SYM --- symmetry flag for nonlinearity;
%   APPLY_LP --- low-pass filter flag (optional);
%   PLOT_SIG --- flag to plot input and bandlimited signals (optional).
% output:
%   snr --- signal-to-noise ratio;
%   nmr --- noise-to-mask ratio.
function [snr, nmr] = metrics_sin(sig, fs, f0, IS_SYM, varargin)
    % bit truncation for NMR calculation
    bit = 24;

    % obtain bandlimited version of a signal
    [sig, sig_lim, alias] = bandlimit_signal(sig, fs, f0, IS_SYM, varargin{:});

    % calclate SNR via Parseval's theorem
    snr = 10*log10(sum(sig_lim.^2)./sum(alias.^2));

    % calculate NMR
    nmr = calc_nmr(sig, sig_lim, fs, bit);
end
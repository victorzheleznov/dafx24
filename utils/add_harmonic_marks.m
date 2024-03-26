% add harmonic marks to a spectrum plot of a nonlinearly processed sine signal
% input:
%   Y --- FFT;
%   fig_fft --- FFT plot;
%   fs --- sampling rate [Hz];
%   f0 --- fundamental frequency [Hz];
%   IS_SYM --- symmetry flag for nonlinearity.
function [] = add_harmonic_marks(Y, fig_fft, fs, f0, IS_SYM)
    % calculate harmonics
    num_harmonics = floor(0.5*fs/f0);
    if IS_SYM == true
        % even harmonics for symmetric nonlinearity
        harmonics = f0 * (1:2:num_harmonics).';
    else
        % even + odd harmonics for asymmetric nonlinearity
        harmonics = f0 * (1:num_harmonics).';
    end
    
    % extract magnitude and phase
    NFFT_2 = length(Y);
    NFFT = 2*(NFFT_2-1);
    bins = round(harmonics*NFFT/fs) + 1;

    % plot
    figure(fig_fft); hold on;
    fvec = (0:fs/NFFT:fs/2).';
    plot(fvec(bins), 20*log10(abs(Y(bins))./max(abs(Y))), 'kx', 'Linewidth', 0.5);
end
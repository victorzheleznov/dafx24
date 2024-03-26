% calculate and plot FFT
% input:
%   y --- input signal;
%   Fs --- sampling rate [Hz];
%   window --- window;
%   NORM_MAG --- normalise magnitude (optional);
%   fig_fft --- figure handle (optional);
%   clr --- line color (optional).
% output:
%   Y --- FFT;
%   fig_fft --- FFT plot.
function [Y, fig_fft] = myfft(y, Fs, window, varargin)
    % check for 0 dB normalisation
    if length(varargin) >= 1
        NORM_MAG = varargin{1};
    else
        NORM_MAG = true;
    end

    % check prespecified figure
    if length(varargin) >= 2
        fig_fft = varargin{2};
        USE_FIG = true;
    else
        fig_fft = figure;
        USE_FIG = false;
    end

    % check prespecified color
    if length(varargin) >= 3
        clr = varargin{3};
    end
    USE_COLOR = exist('clr', 'var');

    % calculate FFT parameters
    N = length(y);
	NFFT = 2^(ceil(log(N)/log(2)));
	NFFT_2 = NFFT / 2 + 1;

    % generate window
    if strcmp(window, "hann")
        win = hann(N, 'periodic');
    elseif strcmp(window, "hamming")
        win = hamming(N, 'periodic');
    elseif strcmp(window, "blackman")
        win = blackman(N, 'periodic');
    end
        
    % calculate FFT
    Y = fft(y.*win, NFFT);
	Y = Y(1:NFFT_2);

    % plot FFT
    figure(fig_fft);
    fvec = (0:Fs/NFFT:Fs/2).';
    if NORM_MAG == true
        YdB = 20*log10(abs(Y)./max(abs(Y)));
        max_dB = 0;
    else
        YdB = 20*log10(abs(Y));
        max_dB = max(YdB);
    end
    plt = plot(fvec, YdB, 'k', 'Linewidth', 0.5);
    if USE_COLOR == true
        set(plt, 'Color', clr);
    end
    if USE_FIG == false
        xlim([0 Fs/2]);
        ylim([-60+max_dB, 5+max_dB]);
        xlabel('Frequency [Hz]', 'Interpreter', 'latex');
        ylabel('Magnitude [dB]', 'Interpreter', 'latex');
        title('Spectrum', 'Interpreter', 'latex');
        ax_spec = fig_fft.CurrentAxes;
        set(ax_spec.XAxis, 'TickLabelInterpreter', 'latex');
        set(ax_spec.YAxis, 'TickLabelInterpreter', 'latex');
    end
end
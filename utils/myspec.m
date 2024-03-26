% create a spectogram plot of an input signal
% input:
%   x --- mono input signal;
%   Fs --- sampling rate [Hz];
%   N --- frame length;
%   O --- overlap factor (between 0 and 1);
%   window --- window;
%   PLOT_SPEC --- plot flag (optional).
% output:
%   STFT --- short time Fourier transform;
%   fig_spec --- figure handle (optional).
function [STFT, varargout] = myspec(x, Fs, N, O, window, varargin)
    % check plot flag
    if length(varargin) >= 1
        PLOT_SPEC = varargin{1};
    else
        PLOT_SPEC = true;
    end

    % find hop size
    HA = round(N - O*N);

    % generate window
    if strcmp(window, "hann")
        win = hann(N, 'periodic');
    elseif strcmp(window, "hamming")
        win = hamming(N, 'periodic');
    elseif strcmp(window, "blackman")
        win = blackman(N, 'periodic');
    end

    % calculate number of frames
    L = length(x);
    NF = ceil(L/HA);
    x = [x; zeros((NF-1)*HA+N-L,1)];
    
    % STFT size
    NFFT = 2^(ceil(log(N)/log(2))); % next power of 2
    NFFT_2 = NFFT / 2 + 1;

    % calculate STFT
    STFT = zeros(NFFT_2, NF);
    for m = 0:NF-1
        x_frame = win.*x((1:N).'+m*HA);
        X = fft(x_frame, NFFT);
        STFT(:,m+1) = X(1:NFFT_2);
    end
    
    % plot spectogram
    if PLOT_SPEC == true
        fig_spec = figure;
        t = ((0:NF-1).*HA/Fs).';
        freq = (0:Fs/NFFT:Fs/2).';
        STFT_dB = 20*log10(abs(STFT)./max(abs(STFT)));
        imagesc(t, freq, STFT_dB, 'CDataMapping', 'scaled');
        c = colorbar;
        c.Label.String = 'dB';
        colormap hot
        clim([-60 0]);
        xlim([0 t(end)]);
        ylim([0 freq(end)]);
        ax_spec = fig_spec.CurrentAxes;
        set(ax_spec, 'YDir', 'normal');
        set(ax_spec, 'YTick', 0:5000:Fs/2);
        set(ax_spec, 'YTickLabel', 0:5:Fs/2/1e3);
        xlabel('Time [s]', 'interpreter', 'latex');
        ylabel('Frequency [kHz]', 'interpreter', 'latex');
        title_str = sprintf("Spectogram with frame length = $%d$ ms and overlap factor = $%d$\\%%", floor((N/Fs)*1e3), O*1e2);
        title(title_str, 'interpreter', 'latex');
        set(ax_spec.XAxis, 'TickLabelInterpreter', 'latex');
        set(ax_spec.YAxis, 'TickLabelInterpreter', 'latex');
        set(c, 'TickLabelInterpreter', 'latex');
        set(c.Label, 'Interpreter', 'latex');
        varargout{1} = fig_spec;
    end
end
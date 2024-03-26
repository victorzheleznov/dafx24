% calculate noise-to-mask ratio (NMR) using basic version of PEAQ [1]
% input: 
%   sig --- signal at test;
%   sig_ref --- reference signal (without noise);
%   fs --- sample rate [Hz].
% output:
%   nmr_total --- total NMR.
% 
% Implementation is largely based on the code provided in [2] but is given 
% in a less computationally optimised form to facilitate readability. In
% addition, the function is meant to be used at different sample rates 
% compared to [2].
%
% Data boundaries aren't checked since clean data is assumed for reference
% and aliased signals (plus signals are truncated so zero padding doesn't 
% affect the ratio). Time spreading isn't implemented since test signals
% don't have transients.
%
% references:
% [1] BS.1387-2 (05/2023)
%     Method for objective measurements of perceived audio quality
%     https://www.itu.int/rec/R-REC-BS.1387
% [2] Peter Kabal, "An examination and interpretation of ITU-R BS.1387: 
%     Perceptual evaluation of audio quality," Tech. Rep., Department of 
%     Electrical & Computer Engineering, McGill University, 2003
%     report: https://mmsp.ece.mcgill.ca/Documents/Reports/2002/KabalR2002v2.pdf
%     code: https://www-mmsp.ece.mcgill.ca/Documents/Downloads/PQevalAudio/PQevalAudio-v1r0.tar.gz
function [nmr_total, nmr] = calc_nmr(sig, sig_ref, fs, bit)
    % parameters
    N = 2048;                   % frame size
    O = 0.5;                    % overlap factor
    Emin = 1e-12;               % energy threshold

    % calculate STFT parameters
    win = hann(N, 'symmetric'); % window ([1,2] use symmetric window)
    HA = round(N - O*N);        % hop size
    L = length(sig);            % input vector length
    NF = floor((L-N)/HA);       % number of frames (truncating signal)
    
    % outter and middle ear response
    w_sq = ears_response(fs, N);

    % critical bands parameters
    [Nc, fc, fl, fu, dz] = critical_band_param();
    U = group_matrix(fs, N, fl, fu);

    % internal noise
    Ein = internal_noise(fc);

    % normalazing factor for frequency spread
    Bs = freq_spread_fast(ones(Nc,1), ones(Nc,1), fc, dz);

    % quantize audio to signed integer
    Amax = 2^(bit-1);
    sig = round(sig*Amax);
    sig_ref = round(sig_ref*Amax);
    
    % loudness calibration and window scaling compensation
    ftest = 1019.5;
    Lp = 92;
    GL = scaling_factor(N, Amax, ftest/fs, Lp);
    sig = sig*GL;
    sig_ref = sig_ref*GL;

    % frames loop
    Eb_err = zeros(Nc,1);
    Es_ref = zeros(Nc,1);
    for i = 1:NF
        % get current frame
        idx = (1:N).'+(i-1)*HA;
        X = get_stft_frame(sig(idx), win);
        X_ref = get_stft_frame(sig_ref(idx), win);

        % square magnitude
        X_sq = abs(X).^2;
        X_sq_ref = abs(X_ref).^2;

        % outer and middle ear filtering
        Xw_sq = w_sq.*X_sq;
        Xw_sq_ref = w_sq.*X_sq_ref;

        % difference magnitude signal
        Xw_sq_err = Xw_sq - 2*sqrt(Xw_sq.*Xw_sq_ref) + Xw_sq_ref;

        % group into critical bands
        Eb_err(:,i) = max(U*Xw_sq_err, Emin);
        Eb_ref = max(U*Xw_sq_ref, Emin); 

        % add internal noise
        E_ref = Eb_ref + Ein;

        % frequency spreading
        Es_ref(:,i) = freq_spread_fast(E_ref, Bs, fc, dz);
    end

    % calculate noise-to-mask
    gm = mask_offset(Nc, dz);
    M = gm.*Es_ref;           % masking threshold
    nmr = Eb_err./M;
    nmr_total = 10*log10(mean(nmr,'all'));
end

% get STFT frame
% input:
%   x --- current input;
%   win --- window;
% output:
%   X --- FFT of current frame.
function X = get_stft_frame(x, win)
    N = length(win);
    x_frame = win.*x;
    X = fft(x_frame, N);
    X = X(1:N/2+1);
end

% calculate mask offset (`PQ_MaskOffset` from [2])
% input:
%   Nc --- number of bands;
%   dz --- Bark scale band step.
% output:
%   gm --- mask offset.
function gm = mask_offset(Nc, dz)
    % the amount in dB by which the masking threshold lies below the time-frequency spread Bark energy
    mdb = @(k) (k<=12/dz).*3 + (k>12/dz).*(0.25*k*dz);
    
    % weighting vector for masking threshold in energy units
    idx = (0:Nc-1).';
    gm = 10.^(-mdb(idx)/10);
end

% calculate parameters for time domain spreading (`PQtConst` from [2])
% input:
%   fs --- sample rate [Hz];
%   N --- FFT size;
%   fc --- center band frequencies (column).
% output:
%   alpha --- parameters for time filtering in each band.
function alpha = time_param(fs, N, fc)
    % frame rate (frames are updated every N/2 samples for assumed 50% overlap)
    fss = fs / (N/2);

    % time parameters for smallest case and 100 Hz
    t100 = 3e-2;
    tmin = 8e-3;

    % time parameters for every band
    t = tmin + (100./fc).*(t100-tmin);
    alpha = exp(-1./(fss*t));
end

% calculate frequency spread effect across critical bands using regular
% spreading function definition
% input:
%   E --- energy in critical bands ("pitch patterns" [2]);
%   Bs --- normalising factor (calculated for E = 1 across all bands);
%   fc --- center band frequencies;
%   dz --- Bark scale band step.
% output:
%   Es --- energy in critical bands after frequency spreading effect
%          ("unsmeared (in time) excitation patterns" [2])
function Es = freq_spread(E, Bs, fc, dz)
    % number of bands
    Nc = length(E);

    % power law for addition of spreading
    e = 0.4;

    % spreading function (indexes are assumed in range from 0 to Nc-1)
    Sdb = @(i,l,E) ((i<=l).*27 + (i>l).*(-24-230./fc(l+1)+2*log10(E))).*(i-l)*dz; % spreading function in dB scale
    S = @(i,l,E) 10.^(Sdb(i,l,E)/10);                                             % regular spreading function
    A = @(l,E) sum(S(0:Nc-1,l,E));                                                % normalisation factor
    Sn = @(i,l,E) S(i,l,E)./A(l,E);                                               % normalised spreading function

    % bands loop
    Es = zeros(Nc,1);
    for i = 0:Nc-1
        % spreading loop
        for l = 0:Nc-1
            Es(i+1) = Es(i+1) + (E(l+1)*Sn(i,l,E(l+1)))^e;
        end
        Es(i+1) = Es(i+1)^(1/e);
        Es(i+1) = Es(i+1)/Bs(i+1);
    end
end

% calculate frequency spread effect across critical bands using recursive
% implementation (`PQ_SpreadCB` from [2])
% input:
%   E --- energy in critical bands ("pitch patterns" [2]);
%   Bs --- normalising factor (calculated for E = 1 across all bands);
%   fc --- center band frequencies;
%   dz --- Bark scale band step.
% output:
%   Es --- energy in critical bands after frequency spreading effect
%          ("unsmeared (in time) excitation patterns" [2])
function Es = freq_spread_fast(E, Bs, fc, dz)
    % number of bands
    Nc = length(E);

    % power law for addition of spreading
    e = 0.4;

    % allocate storage
    aUCEe = zeros(Nc,1);
    Ene = zeros(Nc,1);
    Es = zeros(Nc,1);

    % calculate energy dependent terms
    aL = 10^(-2.7 * dz);
    for m = 0:Nc-1
        aUC = 10^((-2.4 - 23 / fc(m+1)) * dz);
        aUCE = aUC * E(m+1)^(0.2 * dz);
        gIL = (1 - aL^(m+1)) / (1 - aL);
        gIU = (1 - aUCE^(Nc-m)) / (1 - aUCE);
        En = E(m+1) / (gIL + gIU - 1);
        aUCEe(m+1) = aUCE^e;
        Ene(m+1) = En^e;
    end

    % lower spreading
    Es(Nc-1+1) = Ene(Nc-1+1);
    aLe = aL^e;
    for m = Nc-2:-1:0
        Es(m+1) = aLe * Es(m+1+1) + Ene(m+1);
    end

    % upper spreading i > m
    for m = 0:Nc-2
        r = Ene(m+1);
        a = aUCEe(m+1);
        for i = m+1:Nc-1
            r = r * a;
            Es(i+1) = Es(i+1) + r;
        end
    end

    for i = 0:Nc-1
        Es(i+1) = (Es(i+1))^(1/e) / Bs(i+1);
    end
end

% calculate intertal noise in the ear (`PQIntNoise` from [2])
% input:
%   fc --- center band frequencies (column).
% output:
%   Ein --- internal noise value for each band (column).
function Ein = internal_noise(fc)
    f = fc / 1000; % [kHz]
    Edb = 1.456*f.^(-0.8);
    Ein = 10.^(Edb/10);
end

% calculate grouping matrix from FFT bins to critical bands
% input:
%   fs --- sample rate [Hz];
%   N --- FFT length;
%   fl ---  lower frequency edges (column);
%   fu ---  upper frequency edges (column);
% output:
%   U --- grouping matrix (number of bands x number of real FFT bins)
function U = group_matrix(fs, N, fl, fu)
    df = fs/N;
    idx = (0:N/2);
    ku = (2*idx+1)/2;
    kl = (2*idx-1)/2;
    U = max(0, min(fu, ku*df) - max(fl, kl*df)) ./ df; % eq. (13) in [2]
end

% calculate squared outer and middle ears frequency response (`PQWOME` from [2])
% input:
%   fs --- sample rate [Hz];
%   N --- FFT length.
% output:
%   w_sq --- squared ears response (column).
function w_sq = ears_response(fs, N)
    f = (linspace(0, fs/2, N/2+1) / 1000).'; % [kHz]
    Adb = -2.184*f.^(-0.8) + 6.5*exp(-0.6*(f - 3.3).^2) - 0.001*f.^(3.6);
    w_sq = 10.^(Adb/10);
end

% calculate parameters for critical bands (`PQCB` from [2])
% output:
%   Nc --- number of bands;
%   fc --- center frequencies (column);
%   fl --- lower frequency edges (column);
%   fu --- upper frequency edges (column);
%   dz --- band step in Bark scale.
function [Nc, fc, fl, fu, dz] = critical_band_param()
    dz = 1/4;
    fL = 80;    % frequency bnd start [Hz]
    fU = 18000; % frequency band stop [Hz]

    % Bark scale conversion
    B = @(f) 7*asinh(f/650);
    Binv = @(z) 650*sinh(z/7);

    % compute number of bands
    zL = B(fL);
    zU = B(fU);
    Nc = ceil((zU - zL) / dz);
    
    % obtain frequency bands using Bark scale conversion
    %zl = zL + (0:Nc-1) * dz;
    %zu = min(zL + (1:Nc)*dz, zU);
    %zc = 0.5*(zl + zu);
    %fl = Binv(zl);
    %fc = Binv(zc);
    %fu = Binv(zu);

    % frequency bands from BS.1387-2 standard [1,2]
    fl = [  80.000,   103.445,   127.023,   150.762,   174.694, ...
           198.849,   223.257,   247.950,   272.959,   298.317, ...
           324.055,   350.207,   376.805,   403.884,   431.478, ...
           459.622,   488.353,   517.707,   547.721,   578.434, ...
           609.885,   642.114,   675.161,   709.071,   743.884, ...
           779.647,   816.404,   854.203,   893.091,   933.119, ...
           974.336,  1016.797,  1060.555,  1105.666,  1152.187, ...
          1200.178,  1249.700,  1300.816,  1353.592,  1408.094, ...
          1464.392,  1522.559,  1582.668,  1644.795,  1709.021, ...
          1775.427,  1844.098,  1915.121,  1988.587,  2064.590, ...
          2143.227,  2224.597,  2308.806,  2395.959,  2486.169, ...
          2579.551,  2676.223,  2776.309,  2879.937,  2987.238, ...
          3098.350,  3213.415,  3332.579,  3455.993,  3583.817, ...
          3716.212,  3853.817,  3995.399,  4142.547,  4294.979, ...
          4452.890,  4616.482,  4785.962,  4961.548,  5143.463, ...
          5331.939,  5527.217,  5729.545,  5939.183,  6156.396, ...
          6381.463,  6614.671,  6856.316,  7106.708,  7366.166, ...
          7635.020,  7913.614,  8202.302,  8501.454,  8811.450, ...
          9132.688,  9465.574,  9810.536, 10168.013, 10538.460, ...
         10922.351, 11320.175, 11732.438, 12159.670, 12602.412, ...
         13061.229, 13536.710, 14029.458, 14540.103, 15069.295, ...
         15617.710, 16186.049, 16775.035, 17385.420 ];
    fc = [  91.708,   115.216,   138.870,   162.702,   186.742, ...
           211.019,   235.566,   260.413,   285.593,   311.136, ...
           337.077,   363.448,   390.282,   417.614,   445.479, ...
           473.912,   502.950,   532.629,   562.988,   594.065, ...
           625.899,   658.533,   692.006,   726.362,   761.644, ...
           797.898,   835.170,   873.508,   912.959,   953.576, ...
           995.408,  1038.511,  1082.938,  1128.746,  1175.995, ...
          1224.744,  1275.055,  1326.992,  1380.623,  1436.014, ...
          1493.237,  1552.366,  1613.474,  1676.641,  1741.946, ...
          1809.474,  1879.310,  1951.543,  2026.266,  2103.573, ...
          2183.564,  2266.340,  2352.008,  2440.675,  2532.456, ...
          2627.468,  2725.832,  2827.672,  2933.120,  3042.309, ...
          3155.379,  3272.475,  3393.745,  3519.344,  3649.432, ...
          3784.176,  3923.748,  4068.324,  4218.090,  4373.237, ...
          4533.963,  4700.473,  4872.978,  5051.700,  5236.866, ...
          5428.712,  5627.484,  5833.434,  6046.825,  6267.931, ...
          6497.031,  6734.420,  6980.399,  7235.284,  7499.397, ...
          7773.077,  8056.673,  8350.547,  8655.072,  8970.639, ...
          9297.648,  9636.520,  9987.683, 10351.586, 10728.695, ...
         11119.490, 11524.470, 11944.149, 12379.066, 12829.775, ...
         13294.850, 13780.887, 14282.503, 14802.338, 15341.057, ...
         15899.345, 16477.914, 17077.504, 17690.045 ];
    fu = [ 103.445,   127.023,   150.762,   174.694,   198.849, ...
           223.257,   247.950,   272.959,   298.317,   324.055, ...
           350.207,   376.805,   403.884,   431.478,   459.622, ...
           488.353,   517.707,   547.721,   578.434,   609.885, ...
           642.114,   675.161,   709.071,   743.884,   779.647, ...
           816.404,   854.203,   893.091,   933.113,   974.336, ...
          1016.797,  1060.555,  1105.666,  1152.187,  1200.178, ...
          1249.700,  1300.816,  1353.592,  1408.094,  1464.392, ...
          1522.559,  1582.668,  1644.795,  1709.021,  1775.427, ...
          1844.098,  1915.121,  1988.587,  2064.590,  2143.227, ...
          2224.597,  2308.806,  2395.959,  2486.169,  2579.551, ...
          2676.223,  2776.309,  2879.937,  2987.238,  3098.350, ...
          3213.415,  3332.579,  3455.993,  3583.817,  3716.212, ...
          3853.348,  3995.399,  4142.547,  4294.979,  4452.890, ...
          4643.482,  4785.962,  4961.548,  5143.463,  5331.939, ...
          5527.217,  5729.545,  5939.183,  6156.396,  6381.463, ...
          6614.671,  6856.316,  7106.708,  7366.166,  7635.020, ...
          7913.614,  8202.302,  8501.454,  8811.450,  9132.688, ...
          9465.574,  9810.536, 10168.013, 10538.460, 10922.351, ...
         11320.175, 11732.438, 12159.670, 12602.412, 13061.229, ...
         13536.710, 14029.458, 14540.103, 15069.295, 15617.710, ...
         16186.049, 16775.035, 17385.420, 18000.000 ];

    % transpose
    fl = fl.';
    fc = fc.';
    fu = fu.';
end

% calculate sclaing for loudness and Hann window (`PQ_GL` from [2])
% input:
%   N --- frame length;
%   Amax --- maximum amplitude;
%   fn --- normalised frequency;
%   Lp --- sound pressure level.
% output:
%   GL --- scaling constant.
function GL = scaling_factor(N, Amax, fn, Lp)
    W = N - 1;
    gp = peak_factor(fn, N, W);
    GL = 10^(Lp/20) / (gp * Amax/4 * W);
end

% calculate peak factor (`PQ_gp` from [2])
% (since sinusoid can fall between FFt bins the largest bin value will be 
% the peak factor times the peak of the continuous response)
% input:
%   fn --- normalised input frequency;
%   N --- frame length.
% output:
%   gp --- peak factor.
function gp = peak_factor(fn, N, W)
    % distance to the nearest DFT bin
    df = 1 / N;
    k = floor(fn / df);
    dfN = min((k+1)*df - fn, fn - k*df);

    dfW = dfN * W;
    gp = sin(pi * dfW) / (pi * dfW * (1 - dfW^2));
end
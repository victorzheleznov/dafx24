% generate linear sine sweep
% input:
%   t_dur --- sweep duration [sec];
%   fs --- sample rate [Hz];
%   f0 --- lowest sweep frequency [Hz];
%   f1 --- highest sweep frequency [Hz];
%   amp --- maximum amplitude of the sweep.
% output:
%   t --- time vector;
%   x --- sweep signal.
function [t,x] = gen_sine_sweep(t_dur, fs, f0, f1, amp)
    % calculate sweep parameters
    t = (0:1/fs:t_dur).';

    % generate sweep
    sweep_arg = 2*pi*f0*t + 2*pi*t.*(t/t_dur)*(f1-f0)/2;
    x = amp*sin(sweep_arg);
    len_x = length(x);

    % add fade-in/out
    t_fade = 0.1;
    len_fade = t_fade*fs;
    fade = 0.5*(1 - cos(2*pi.*(0:len_fade-1)./(2*len_fade))).';
    fade = [fade; ones(len_x-2*len_fade, 1); flipud(fade)];
    x = x.*fade;
end
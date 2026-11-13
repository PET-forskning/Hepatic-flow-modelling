function c_out = convexvarc(t_out, par, t, f)
%CONVEXVARC  Convolve a sampled curve with a mono-exponential and frame-average.
%
%   c_out = CONVEXVARC(t_out, par, t, f) computes
%
%       c(t) = f(t) (*) par(1)*exp(-par(2)*t)
%
%   where (*) denotes convolution, and returns c(t) averaged over the PET
%   frames whose mid-times are given in t_out.
%
%   Inputs
%     t_out  Frame mid-times of the output curve (min). The frame boundaries
%            are reconstructed from the mid-times assuming contiguous frames
%            with the first frame starting at t = 0.
%     par    [a b ...]: amplitude a and rate constant b (1/min) of the
%            exponential. Further elements are ignored.
%     t      Sample times of the input curve f (min). A time shift (e.g. a
%            delay) can be applied by passing t + delay.
%     f      Input curve (e.g. arterial activity concentration).
%
%   Output
%     c_out  Column vector with the frame-averaged convolution, one value
%            per element of t_out.
%
%   The convolution is evaluated numerically on a fine, regular time grid
%   (0.01 min, or 0.6 when the time axis appears to be in seconds) after
%   linear interpolation of the input curve, and then averaged within each
%   frame.

a = par(1);
b = par(2);

t = t(:);
f = f(:);

% Make sure the input curve starts at t = 0 with zero activity
if min(t) > 0
    t = [0; t];
    f = [0; f];
end

% Fine time grid for the numerical convolution
if max(t) < 100
    step_size = 0.01;   % time in minutes
else
    step_size = 0.6;    % time in seconds
end
t0 = 0:step_size:max(t);

% Reconstruct frame boundaries from the mid-times (first frame starts at 0)
t_frame = zeros(length(t_out)+1, 1);
for k = 1:length(t_out)
    t_frame(k+1) = 2*t_out(k) - t_frame(k);
end

% Interpolate input onto the fine grid, convolve, and frame-average
f0    = interp1(t, f, t0, 'linear', 'extrap');
c_fine = conv(f0, a*exp(-b*t0));
c_fine = c_fine(1:length(t0));
c_out  = frame_average(t_frame, t0, c_fine)' * step_size;

end


function c_frame = frame_average(t_frame, t_in, c_in)
% Mean of c_in over the samples t_in that fall within each frame
% [t_frame(i), t_frame(i+1)). Frames without samples are set to 0.
c_frame = zeros(1, length(t_frame)-1);
for i = 1:length(c_frame)
    % For increasing frame times this selects t_frame(i) <= t_in < t_frame(i+1)
    index = (t_in >= t_frame(i)) == (t_in < t_frame(i+1));
    if any(index)
        c_frame(i) = sum(c_in(index)) / sum(index);
    end
end
end

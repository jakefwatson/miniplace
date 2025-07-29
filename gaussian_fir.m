%% gaussian_fir.m Generate FIR coefficients for a Gaussian low-pass filter.
%
% Inspired by: William R. Rose's "1D Gaussian lowpass filter" (2006),
% MATLAB Central File Exchange (File ID: 12606)
% https://www.mathworks.com/matlabcentral/fileexchange/12606-1d-gaussian-lowpass-filter

function b = gaussian_fir(sr, fco)
%   b = gaussian_fir(sr, fco) returns FIR filter coefficients for
%   a Gaussian low-pass filter with sample rate sr (Hz) and −3 dB
%   cut-off frequency fco (Hz). The resulting FIR filter is symmetric
%   and normalised to have unity gain at DC.
%
%   Example:
%       b = gaussian_fir(1000, 50);
%       y = filter(b, 1, x);  % Apply to signal x
%

    % Compute standard deviation of Gaussian
    sigma = 1 / (2 * pi * fco);     % seconds

    % Number of samples to cover ±3 sigma
    N = ceil(3 * sigma * sr);
    n = -N:N;                       % Sample indices
    t = n / sr;                     % Time vector (in seconds)

    % Compute Gaussian impulse response
    g = exp(-(t.^2) / (2 * sigma^2));

    % Normalize for unit gain (DC response = 1)
    b = g / sum(g);
end

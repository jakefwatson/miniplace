%% noisegener.m Function for generation of normally distributed noise for simulated mini analysis
% This script can be used to generate realistic simulated recording noise for 
% miniature PSC (mini/mPSC) data simulation. 
% It was written for Greger & Watson 2024 (doi: https://doi.org/10.1101/2024.10.26.620084 )
%
%     noisegener.m Copyright (C) 2025 Jake Watson
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%     You should have received a copy of the GNU General Public LicenseAdd commentMore actions
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% This function generates simulated recording noise of a defined length,
% with some control over noise properties.
%
% Function INPUTS
%
% dur - Duration of required simulated recording (seconds)
% sf - Sampling frequency of simulated data (Hz)
% type - Type of noise required from recording. Options are 'white'
%       'filtered' and 'mixed'.
%       'white' - uses simple white noise, a random signal with equal
%           intensity across all frequencies.
%       'filtered' - uses white noise filtered with a low pass filter set
%       to a defined frequency (filt).
% filt - Sets the frequency of low pass filtering for recording noise (Hz; default =
%       1000 Hz).  If no filtering is required, enter 0
% 
% Function OUTPUTS
%
% noise - is a vector of simulated noise
% time - is a vector of time corresponding to 'noise'
%
% This code uses '1D Gaussian lowpass filter, gaussfiltcoef' by William Rose (2006)
% William Rose (2006). 1D Gaussian lowpass filter (https://www.mathworks.com/matlabcentral/fileexchange/12606-1d-gaussian-lowpass-filter), MATLAB Central File Exchange.
%
% Pink noise was introduced following Spectral Audio Signal Processing by Julius O. Smith
% Smith, J.O. Spectral Audio Signal Processing, http://ccrma.stanford.edu/~jos/sasp/,
% Online book, 2011 edition.
%
% Written by Jake Watson for Greger & Watson 2024 (doi: doi.org/10.1101/2024.10.26.620084)

function [noise,time] = noisegener(dur,sf,type,filt)
samples = dur .* sf;
time = linspace(0,dur,samples)';
if strcmp(type, 'white') == 1
    % Change the standard deviation of white noise here. 1.4 is default.
    sd = 1.4; 
    noise = 0 + sd.*randn(samples,1);

elseif strcmp(type, 'filtered') == 1
    % Note - SD values reported here are pre-filtering and do not reflect
    % the final SD of simulated data. Actual simulated recording SD should
    % be measured from the output of noisegener.
    if sf == 100000
        sd = 12;
    elseif sf == 10000 
        sd = 3; % 3 gives standard SD of resulting simulation (final SD = approx. 1.4 pA)
    else
        % Noise has not been optimised for sf values other than 10 kHz 
        % It is recommended that you adjust the sd of noise for your own
        % simulation purposes.
        sd = 12;
    end
    white = 0 + sd.*randn(samples,1);
    filtset = gaussian_fir(sf,filt);
    noise = filter(filtset,1, white);  
    
elseif strcmp(type,'mixed') == 1
    % Recommended values for simulated noise levels are:
    %   Normal = 1* Medium-high = 1.5* High = 2*
    white = 1.5*randn(samples,1); % Change scaling factor to alter noise level.
    
    % The following follows Smith 2011 (see above)
    B = [0.049922035 -0.095993537 0.050612699 -0.004408786];
    A =  [1 -2.494956002   2.017265875  -0.522189400];
    nT60 = round(log(1000)/(1-max(abs(roots(A)))));
    v = randn(1,samples+nT60);
    x = filter(B,A,v);
    x = x(nT60+1:end);
    mixer = (x-mean(x)).*3;
    
    combo = mixer' + white;
    
    filtset = gaussian_fir(sf,filt);
    combo = filter(filtset,1, combo);  
    noise = combo .* (10/3.5);
    
end

end

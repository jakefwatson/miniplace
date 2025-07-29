%% miniranscaler.m Function for generation of simulated mini data.
% This script can be used to generate realistic simulated miniature PSC
% (mini/mPSC) data.  It was written for Greger & Watson 2024
% (doi: https://doi.org/10.1101/2024.10.26.620084 )
%
%     miniranscaler.m Copyright (C) 2025 Jake Watson
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
% This script can be used for simulation of mini/mPSC/sPSC data for
% exploration of the effect of event parameters on detection and resulting
% data distributions. 
%
% Without further  modification, input events are simulated using a
% biexponential function, formed of a rising exponential (τrise) and a
% decaying exponential (τdecay), filling a 70 ms window. Rise and decay
% taus are randomly selected from a lognormal distribution of realistic
% values (τrise μ: 0.2, σ: 0.25; τdecay μ: 1.7, σ: 0.4) with with τrise
% constrained between 0.3 < τ < 2.5 ms, and τdecay between 1 < τ < 25 ms.
% Resulting curves are scaled to a peak amplitude randomly sampled from a
% lognormal distribution of realistic peak amplitudes (peak μ: 1.6, σ: 0.8).
% For scaled amplitude datasets, scaling factors (amplitude addition or
% multiplication) were applied to target peak amplitudes prior to curve scaling.
%
% Function INPUTS
%
% dur - Duration of required simulated recording (seconds)
% freq - Required frequency of events to be input (Hz)
% sf - Sampling frequency of simulated data (Hz)
% peaklimit - The minimum limit of peak amplitude for event simulation (pA).
%       Note: The event distribution will be cut at this value, therefore using
%       this function will alter the distribution of input events.
% type - Type of noise required from recording. Options are 'white'
%       'filtered' and 'mixed'. See noisegener.m for further information
% filt - Sets the frequency of low pass filtering for recording noise (Hz; default =
%       1000 Hz).
% modifier - Determines the input event modifier required.
% scal - Defines the scaling magnitudes required. Input values as a vector
% in one row (e.g. [ 1 2 3 4 5]).
% The modifier options and appropriate scal values are:
%   'none' - No modifier applied. Individual simulated recordings produced
%       For 'none', scal should be set to [0];
%   'add' - Parallel simulated recordings are generated, where events in 
%       each simulation have been increased by addition of a defined magnitude
%       For 'add', scal values denote the amplitude values to be added to
%       all event amplitudes. 0 is required for basal amplitude runs.
%       e.g. [0 2 4] produces three parallel simulated recordings with embedded
%       events of randomly selected amplitudes (0), all increased by 2 pA
%       and 4 pA in separate simulations.       
%   'scale' - Parallel simulated recordings are generated, where the
%       amplitude of events in each simulation have been scaled by a defined
%       multiplication factor.
%       For 'scale', scal denotes the multiplication factors applied to
%       embedded event amplitudes. 1 is required for basal amplitudes.
%       e.g. [1 1.2 1.4] produces three parallel simulated recordings with
%       embedded events of randomly selected amplitudes (0), all scaled by
%       1.2-fold and 1.4-fold in separate simulations.
%   'biol' - Parallel simulated recordings are generated, where events have
%       been scaled with a biologically-inspired scaling ratio. Smaller events
%       are scaled to a larger degree than larger events, according to the
%       equation: Scaling factor = A(s.0.9^A + 1); where A is the initial
%       amplitude, and s is the magnitude of scaling.
%       For 'biol', scal sets the values of s. s = 0 denotes no
%       potentiation (baseline distribution). Eecommended values of s
%       are 1-3. e.g. scal = [0 1 2 3];
%   'kin' - Input event kinetics are fixed to set values. Using other
%       modifier values, kinetics are randomly selected from a distribution for
%       each event. Using 'kin', all events have fixed kinetics of:
%           τrise = x ms
%           τdecay = 5 + 10x ms
%       where x is a user input value/s from 'scal'.
%
% Function OUTPUTS
%
% recording - contains simulated recording data (rows are individual data points,
%   columns contain a recording for each 'scal'.
% noisefree - contains simulated recording (as per 'rec'), but without noise
%   embedding, allowing visulisation of encoded events.
% time - contains a vector of time values (units: seconds) for 'rec' and 'noisefree'
% events - contains vectors of input events without noise
% minilocs - contains a list of sample numbers at which events have been
%   included in the recording
% noise - contains the recording noise without events embedded.
%
% Written by Jake Watson for Greger & Watson 2024 (doi: doi.org/10.1101/2024.10.26.620084)
%

function [noisefree,recording,time,events,minilocs,noise] = miniranscaler(dur,freq,sf,peaklimit,type,filt,modifier,scal)

outs = size(scal,2);
n = dur .* freq;

%-------- Generate noise and simulated minis----------
[noise] = noisegener(dur,sf,type,filt);
[events] = minigener(n,peaklimit,sf,modifier,scal);
minidur = size(events,1);

%-------- Plot positions of minis----------

minitimes = (dur - 0.14) .* rand(n,1); % Events cannot be simulated in the last 140 ms to prevent false detection due to cutoff
minilocs = floor(minitimes.*sf);

%------ Noisefree recording generation ----------

sampno = size(noise,1);
noisefree = zeros(sampno,outs);

for g = 1:outs
    for j = 1:n
        noisefree(minilocs(j):(minilocs(j) + minidur-1),g) = noisefree(minilocs(j):(minilocs(j) + minidur -1),g) + events(:,j,g);
    end
end
        
%-------- Superimpose noise ----------

recording = noisefree;
for i = 1:outs
        recording(:,i) = recording(:,i) + noise(:,1);
end
time = (1/sf):(1/sf):dur;
time = time';
end

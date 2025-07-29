%% minisimulation.m - Script for generation and analysis of simulated mini data
% This script is can be used for simulation of mini data with control over
% input event parameters and simulated recording noise. It was written for
% Greger & Watson 2024 (doi: https://doi.org/10.1101/2024.10.26.620084 )
%
%     minisimulation.m Copyright (C) 2025 Jake Watson
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
% data distributions. It offers multiple means to modify input event
% amplitudes or kinetics in resulting data. Simulated recordings are
% subject to event detection using minidet.m, developed by Alois Schlögl
% and available as part of the Biosig project (https://biosig.sourceforge.net/)
%
% Without further script modification, input events are simulated using a
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
% This script requires user input in the initial section to provide:
% cells - The number of repetitions to be performed. i.e. simulated 'cells'
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
% dur - Duration of required simulated recording (seconds)
% freq - Required frequency of events to be input (Hz; default = 1 Hz)
% sf - Sampling frequency of simulated data (Hz; default = 10000 Hz).
%       Values other than 10 kHz have not been tested or optimised.
% peaklimit - The minimum limit of peak amplitude for event simulation.
%       Note: The event distribution will be cut at this value, therefore using
%       this function will alter the distribution of input events (pA).
% filt - Sets the frequency of filtering for recording noise (Hz; default =
%       1000 Hz).
%
% Further parameters such as the properties of recording noise, the
% distributions of input event amplitudes and kinetics etc. must be changed
% at specific locations in the respective scripts.
%
% The outputs of simulation can be found in the following variables:
% MINIRESULTS - A structure containing parameters of input and detected
% events from simulation and minidet.m detection.
% MINIRESULTS.Freq - Lists the frequency of detected events in each
% simulated recording. Each row contains data from one 'scal' value, and
% each column contains the individual 'cells' (repeats) for that condition.
% MINIRESULTS.Amplitude - Lists the measured amplitude of detected events.
% Each column contains a list of amplitudes for each 'scal' value.
% Data for different 'cells' are included in the third dimension.
% MINIRESULTS.MeanAmp - Calculates the mean peak amplitude from MINIRESULTS.Amplitude
% MINIRESULTS.InputAmps - Lists the peak amplitudes of input events,
% independent of whether these events were detected. This variable is
% structured equivalently to MINIRESULTS.Amplitude.
%
% rec - contains simulated recording data (rows are individual data points,
% columns contain a recording for each 'scal'. 'cells' are repeated in the
% third dimension).
% 
% noisefree - contains simulated recording (as per 'rec'), but without noise
% embedding, allowing visulisation of encoded events.
%
% time - contains a vector of time values (units: seconds) for 'rec' and
% 'noisefree'
%
% minilocs - contains a list of encoded event times. These are the same for
% all 'scal' values, therefore each row contains event times, and columns
% show events for each repetition ('cells').
%
% Written by Jake Watson for Greger & Watson 2024 (doi: doi.org/10.1101/2024.10.26.620084)
%

%% USER INPUT SECTION
cells = 1;          % number of repetitions. Default = 1
modifier = 'add';   % modifier (see header). Default = 'none'
scal = [ 0 2 ];         % scaling values for 'modifier'.
dur = 100;            % Duration of simulated recording (s)
freq = 1;               % Frequency of encoded events (Hz)
sf = 10000;              % Sampling frequency (Hz)
peaklimit = 0;            % Minimum peak amplitude (pA)
filt = 1000;               % Frequency of noise filtering (Hz)

%% Generation of simulated recordings (miniranscaler.m)
% Recording noise is automatically set as 'mixed noise' (pink noise). This
% can be altered in the input to miniranscaler, below. Alternative options
% are 'white' and 'filtered'. See 'noisegener.m' for more information.

scales = size(scal,2);
for j = 1:cells
     [noisefree(:,:,j),rec2(:,:,j),time,fakes(:,:,:,j),minilocs(:,j)] = miniranscaler(dur,freq,sf,peaklimit,'mixed',filt,modifier,scal);
end
clear j

%% Detection of simulated events (minidet.m)
% Simulated recordings are next subject to event detection using minidet.m
% For more information on minidet, please consult the Biosig project (see header).

EventLocs = nan(1,scales,cells);
EventPlaces = nan(round(dur*freq),scales,cells);
MINIRESULTS.Freq = nan(scales,cells);

for i = 1:cells
    for j = 1:scales
        lpdata = rec(:,j,i); % Data input
        refract = 0.004*sf;  % Refractory period (Default set at 0.004 s)
        rCrit   = 0.5;       % Template search parameters.
        aCrit   = 10;
    % It is recommended that these parameters be adjusted for optimal detection of recorded events.
    % Default values are 'rCrit = 0.5; aCrit = 10'
        res     = minidet(lpdata,[],rCrit,aCrit,refract);
        pos1    = res.tEventList;
    % 'res' is the output of minidet for simulated data with noise
        nfdata  = noisefree(:,j,i);
        resnf   = minidet(nfdata,[],rCrit,aCrit,refract);
        pos2    = resnf.tEventList;
    % minidet is also run on noisefree recording, with output 'resnf'
        if size(EventLocs,1) < size(pos1,1)
            diff = size(pos1,1) - size(EventLocs,1);
            EventLocs = [EventLocs ; nan(diff,scales,cells)];
        end
    % EventLocs records the list of all detected events, embedded in a NaN array.
        EventLocs(1:size(pos1,1),j,i) = pos1;
        EventPlaces(1:size(pos2,1),j,i) = pos2;
        MINIRESULTS.Freq(j,i) = size(pos1,1)./dur;
    % EventPlaces records the equivalent of EventLocs for noisefree detection.
    end
end

%% Measuring event amplitudes (eventfitter.m)
% Next, events are extracted and fit with a biexponential curve to estimate
% peak amplitude.

rangepull = [-0.005 0.03]; % Range of events to be fitted with biexponential.
% Event two times, start and end of extracted window, in seconds, relative to
% the mid point of detected event rise. Recommended values are [-0.005 0.03]

EventExtract = nan(round((rangepull(2)-rangepull(1)).*sf)+1,size(EventLocs,1),scales,cells);
MINIRESULTS.Amplitude = nan(size(EventLocs,1),scales,cells);

failcount = 0; % failcount will report the number of events for which curve
% fitting was not possible. Should curve fitting fail, fit parameters need
% to be adjusted, or an alternative peak amplitude measurement method should
% be applied.
for i = 1:cells
for j = 1:scales
for k = 1:size(EventLocs,1)
    if isnan(EventLocs(k,j,i)) % No event may be present if conditions have different numbers of events.
        continue
    end
    
    samplepull = [EventLocs(k,j,i)+(rangepull(1).*sf) EventLocs(k,j,i)+(rangepull(2).*sf)];
    EventExtract(:,k,j,i) = rec(round(samplepull(1)):round(samplepull(2)),j,i);
    
    try
        MINIRESULTS.Amplitude(k,j,i) = eventfitter(time(1:size(EventExtract,1)),EventExtract(:,k,j,i));
    catch
        fprintf('Curve fitting failed for event %d, scale %d, cell %d.\n', k, j, i);
        failcount = failcount +1;
        % Alternative peak amplitude measurement for failed curve fitting
        % can be implemented within catch (e.g. below), but must be considered with
        % caution not to affect results. The example below finds the mean
        % of the five datapoints surrounding the minimum value within 5 ms
        % of the event location (midpoint of rise phase). This is
        % reliant on the peak of a detected event occurring within 5 ms of
        % the rising phase, and to contain the minimum value, which is not
        % always the case for small events with surrounding noise.
        %
        % v = EventExtract(round(-rangepull(1)*sf+1):round((-rangepull(1)+0.005)*sf), k, j, i);
        % MINIRESULTS.Amplitude(k,j,i) = mean(v(max(1,find(v==min(v),1)-2):min(end,find(v==min(v),1)+2)));
    end
end
end
fprintf('CELL %d DONE \n', i)
end

MINIRESULTS.MeanAmp = squeeze(nanmean(MINIRESULTS.Amplitude,1));
MINIRESULTS.InputAmps = nan(size(fakes,2),size(fakes,3),size(fakes,4));
MINIRESULTS.InputAmps = -squeeze(min(fakes,[],1))';

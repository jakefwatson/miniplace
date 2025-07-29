%% minigener.m Generates simulated mPSCs with random parameters
% This script generates simulated mPSC events with amplitude, rise, and
% decay kinetics randomly selected from defined distributions. It is used
% in conjunction with miniranscaler.m and minisimulation.m for data
% simulation.
%
%     minigener.m Copyright (C) 2025 Jake Watson
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
% Function INPUTS
%   n - number of minis desired
%   peaklimit - peak of smallest possible mini
%   sf - sampling frequency of "recording"
%   modifier - Determines the input event modifier required.
%   scales - Defines the scaling magnitudes required. Input values as a vector
%       in one row (e.g. [ 1 2 3 4 5]).
%   The modifier options and appropriate scal values are:
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
%   events - An array of n columns of minis, 70 ms long
%   time - A corresponding time vector   
%          
% Written by Jake Watson for Greger & Watson 2024 (doi: doi.org/10.1101/2024.10.26.620084)

function [events,time] = minigener(n,peaklimit,sf,modifier,scales)
%% --- Initialising modifiers

if strcmp(modifier , 'none') == 1
    outs = 1;
    modcode = 1;
elseif strcmp(modifier , 'add') == 1
    outs = size(scales , 2);
    modcode = 2;
elseif strcmp(modifier , 'scale') == 1
    outs = size(scales , 2);
    modcode = 3;
elseif strcmp(modifier , 'biol') ==1
    outs = size(scales , 2);
    modcode = 4;
elseif strcmp(modifier , 'kin') == 1
    outs = size(scales , 2);
    modcode = 5;
end

%% ----- Generating event risetimes ---------
fakerise = lognrnd(0.2,0.25,[100*n 1]); % Distribution of rise time constants can be altered as desired.
fakerise(fakerise>2.5) = []; % Maximal rise time constant cutoff
fakerise(fakerise<0.3) = []; % Minimum rise time constant cutoff
fakeriseset = fakerise(1:n,1);

%% ----- Generating event decay times -------
fakedecay = lognrnd(1.7,0.4,[100*n 1]); % Distribution of decay time constants can be altered as desired.
fakedecay(fakedecay<1) = []; % Minimum decay time constant cutoff
fakedecay(fakedecay>25) = []; % Maximum decay time constant cutoff
fakedecayset = fakedecay(1:n,1);
if modcode == 5 % If modifier is 'kim' kinetics - values are set rather than chosen from a distribution.
    fakeriseset(:) = 0;
    fakedecayset(:) = 5;
end

%% ------ Generating event peak amplitudes ------ 
fakepeak = lognrnd(1.6,0.8,[100*n 1]);
fakepeak(fakepeak < peaklimit) = []; % Accounting for peaklimit
%fakepeak(fakepeak > 30) = []; % A maximum amplitude can also be set.
fakepeakset = fakepeak(1:n,1);
fakepeakscaled = zeros(size(fakepeakset,1),outs);

if modcode == 1
    fakepeakscaled(:,1) = fakepeakset(:,1);

elseif modcode == 2 % Implementing additive scaling of peak amplitudes.
    for h = 1:outs
        fakepeakscaled(:,h) = fakepeakset(:,1) + scales(1,h);
    end
    
elseif modcode == 3 % Implementing multiplicative scaling of peak amplitudes.
    for g = 1:outs
        fakepeakscaled(:,g) = fakepeakset(:,1) .* scales(1,g);
    end
elseif modcode == 4 % Implementing biologically inspired scaling of peak amplitudes.
    f1 = @(a,b,c,x) (c.*a.^x+b);
    a = 0.9;
    b = 1;
     for g = 1:outs
         c = scales(1,g);
         fakepeakscaled(:,g) = fakepeakset(:,1).*f1(a,b,c,fakepeakset(:,1));
     end

elseif modcode == 5
     fakepeakscaled(:,1) = fakepeakset(:,1);
end

%% ------ Generating the simulated events --------
samples = floor(sf.* 0.07); % Should events be simulated over a longer time period, this can be altered here.
events = zeros(samples,n);

if modcode ~= 5
T = (1000./sf):(1000./sf):70;
for j = 1:outs
    for i = 1:n
        tau_rise = fakeriseset(i);
        tau_decay = fakedecayset(i);
        
        % The following modification may be required depending on the
        % distribution of rise and decay time constants used:

        %if tau_rise >= tau_decay
            %tau_decay = tau_decay + tau_rise;
        %end
        
        Template = ((1-exp(-T/tau_rise)).*exp(-T/tau_decay));
        TemplatePeak = max(Template);
        Template = Template/TemplatePeak;
        Template = (-fakepeakscaled(i,j).*Template)';

        events(:,i,j) = Template;
    end
end
else
    T = (1000./sf):(1000./sf):140;
    events = zeros(samples*2,n);
    for j = 1:outs
        for i = 1:n
            tau_rise = fakeriseset(i) + scales(1,j);
            tau_decay = fakedecayset(i) + scales(1,j)*10;
            Template = ((1-exp(-T/tau_rise)).*exp(-T/tau_decay));
            TemplatePeak = max(Template);
            Template = Template/TemplatePeak;
            Template = (-fakepeakscaled(i).*Template)';

            events(:,i,j) = Template;
        end
    end
end


%% ------ Generating time --------
time = linspace(0,0.07,(size(events,1)));

end
%% Limiter - calculates the detection limit from two scaled mEPSC datasets
% This script can be used to estimate the mPSC detection limit from two datasets 
% of detected events recorded at different holding potentials.
% It was written for Greger & Watson 2024 (doi: https://doi.org/10.1101/2024.10.26.620084)
%
%     limiter.m Copyright (C) 2025 Jake Watson
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
% By analysing two datasets (data1 and data2) of mPSC peak amplitudes
% recorded at two holding potentials (hold1 and hold2 - in mV) the detection
% limit cutoff amplitude (limit) is returned. This function uses the different
% driving force for ion flow at altered holding potentials to 'visualise'
% mEPSC events below the detection threshold.
%
% INPUTS    - data1 and data2 are two single column vectors containing a list of
%               mEPSC peak amplitudes (units: pA).
%           - hold1 and hold2 are the holding potentials corresponding to data1 and data2
%
% OUTPUT    - limit, the estimated detection limit from broken stick fitting (in pA).
%
% This analysis uses a 'sliding bin' to scan datasets 1 and 2 at a
% resolution of 'fineness' (default 0.1 pA). This style of binning should
% not be used for representation of data, as individual events are
% included in multiple bins.
%
% Limiter estimates the detection limit by comparing the number of events
% detected per bin at hold1 (in data1), with the equivalent mEPSC bins post
% after scaling (data2). For example, mEPSCs between 5 and 6 pA at -60 mV
% holding potential would occur between 7.5 and 9 pA at -90 mV holding potential.
% Given the scaling factor (1.5 in above example) allows detection of
% events at hold2 that are below the detection limit at hold1, the number
% of 'recorded events' (hold1) is compared to 'expected events' (hold2).
%
% 'Recorded' and 'Expected' events follow y = x away from the detection
% limit (high mEPSC amplitudes), however linearity is lost when moving bin
% includes values below the detection limit. Limiter calulates this value
% as the break point of a broken stick regression fit of 'expected - recorded'
% mEPSCs.
%
% The detection limit is not a strict cutoff in biological data analysis;
% smaller events are detected with decreasing liklihood. We define the
% detection limit as the point at which events are evidently lost from
% recorded distributions (i.e. where false negatives are introduced).
% In addition, imperfect voltage clamp in complex neuronal morphologies
% will preclude perfect manipulation of event amplitudes using holding
% potential scaling. Therefore, this method provides an estimate of the
% detection limit, necessary for rigorous data analysis and interpretation,
% but not a definitive measurement. 
%
% Written by Jake Watson for Greger & Watson 2024 (doi: doi.org/10.1101/2024.10.26.620084)

function [limit] = limiter(data1,data2,hold1,hold2)

data1 = abs(data1);
data2 = abs(data2);

scaler = hold2 ./ hold1 ;

low = 0; % set lowest bin (default value: 0 pA)
high = 50; % set highest bin (default value: 50 pA)
binwidth = 5; %change the size of bins (default value: 5 pA)
fineness = 0.1; % change the fineness of scanning (default value: 0.1 pA)
fitwidth = 25; % Fit range for broken stick (default: 25 pA)

scans = binwidth/fineness;
binslow = zeros(scans , (high ./ binwidth) + 1);

for i = 1:scans
    g = (i-1) .* fineness;
    
    binslow(i,:) = linspace((low + g) , (high + g) , (high./binwidth)+1);
      
end
binshigh = binslow .* scaler;

histlow = zeros(scans , high ./ binwidth );
histhigh = histlow;

% ----- Generation of binned data ---------
for i = 1:scans
    histlow(i,:) = histcounts(data1 , binslow(i,:));
    histhigh(i,:) = histcounts(data2 , binshigh(i,:));
end

bins = reshape(binslow(:,1:end-1),[],1);

hold1vals = reshape(histlow,[],1);
hold2vals = reshape(histhigh,[],1);

% ----- Plotting scanned histogram ---------
figure(1)
plot(bins,hold1vals)
hold on
plot(bins,hold2vals)
xlabel('Peak amplitude (pA)')
ylabel('Number of values')
title('Moving bin scan of scaled data')
hold off

differ = (hold2vals - hold1vals);

%% Broken Stick plotting

[~,peakIdx] = max(hold2vals); % Dataset for broken stick fitting removing downslope of curves at 0 - peak pA
maxamp = bins(peakIdx);
maxidx = peakIdx;
endidx = min(length(bins), maxidx + round(fitwidth / fineness));
cutbins = bins(maxidx:endidx);
xdata = cutbins;
ydata = differ(maxidx:endidx);

coefs = @(xbreak) (min(xdata, xbreak) - xbreak) \ ydata;
errfun = @(xbreak) norm((min(xdata, xbreak) - xbreak) * coefs(xbreak) - ydata);

xbreak = fminbnd(errfun, maxamp, high);
C = coefs(xbreak);

yfun = @(x) (min(x, xbreak) - xbreak) * C;

figure(2)
plot(xdata,ydata,'o')
hold on
fplot(yfun, [cutbins(1), cutbins(end)])
xlabel('Amplitude (pA)')
ylabel('Frequency difference')
title('Broken stick estimation of the detection limit')
hold off

limit = xbreak;

%% Plot histograms of data above detection limit

finalbins = xbreak:1:floor(high);

figure(3)
histogram(data1,finalbins)
hold on
histogram(data2,finalbins)
xlabel('Amplitude (pA)')
ylabel('Frequency')
title('Histogram truncated at the detection limit')
hold off

end

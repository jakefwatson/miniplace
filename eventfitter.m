%% eventfitter.m - Script for measurement of mini amplitudes by curve fitting
% This function can be used to fit detected minis (mPSC/sPSC events) with a
% biexponential function for measurement of peak amplitudes. The advantage
% of this approach over simple 'min/max value' measurement is that for
% small events with a relatively large noise component, peak amplitudes are
% less likely to be overestimated due to noise. It was written for
% Greger & Watson 2024 (doi: https://doi.org/10.1101/2024.10.26.620084 )
%
%     eventfitter.m Copyright (C) 2025 Jake Watson
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
% Function INPUTS
% A single recorded event should be input with x (time) and y (current)
% data as ex and why respectively.
%
% Function OUTPUTS
% This function outputs the measured peak amplitude ('peak'), and the
% parameters of fit (a vector 'output').
%
% Written by Jake Watson for Greger & Watson 2024 (doi: doi.org/10.1101/2024.10.26.620084)

function [peak,output] = eventfitter(ex,why) 

x  = ex;
y = why;

f = @(a,b,c,d,x) (a.*(-exp(-(x-d)./b)+exp(-(x-d)./c)).*(x>=d)); 
offered = [-8; 0.001; 0.015; 0.005];
lowlim  = [-1000 ; 0 ; 0.002 ; 0];
uplim   = [0 ; 0.015 ; 0.1 ; x(end)];

f1 = fit(x,y,f,'Start', offered,'lower',lowlim,'upper',uplim); 
output = coeffvalues (f1);

peak = min(f(output(1), output(2), output(3), output(4), x));

end


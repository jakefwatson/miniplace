# miniplace
These scripts can be used for simulation of mini (mPSC) data. It was written for Greger & Watson 2024 [doi:10.1101/2024.10.26.620084](https://doi.org/10.1101/2024.10.26.620084)

Included files:

minisimulation.m - Main script for event simulation

miniranscaler.m - Function for recording simulation, reliant on minigener and noisegener

minigener.m - Generates simulated events

noisegener.m - Generates simulated recording noise

eventfitter.m - Fits detected events with a biexponential function for peak amplitude measurement

gaussian_fir.m - Returns FIR filter coefficients for a Gaussian low-pass filter

limiter.m - function for broken stick estimation of the detection limit from mini datasets recorded at two holding potentials

 
This project uses material from:

Julius O. Smith III, *Spectral Audio Signal Processing*.  
https://ccrma.stanford.edu/~jos/sasp/  

and based some analysis on:
William Rose, *1D Gaussian lowpass filter*.
https://de.mathworks.com/matlabcentral/fileexchange/12606-1d-gaussian-lowpass-filter
MATLAB File Exchange 12606

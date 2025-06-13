% initializeGUIdefaults: Sets the default values for the processing options
% and flag values that influence what the GUI initially shows when
% initialized.
%
%   INPUTS:     NONE - Edit this file to adjust the default values!
%
%   OUTPUTS:
%       procflgs    -   Struct containing values of logical flags for
%                       data processing options
%       params      -   Struct containing numerical values/strings for data
%                       processing specifications
%
function [procflgs,params]=initializeGUIdefaults
procflgs.override = false; %if true, user is able to override acq/proc  
    %parameters pulled in from the dataset
params.pw90 = Inf; %true 90 pulse width to use to calculate saturation amplitudes
procflgs.topproc = false; %if true, data are loaded with processing previously  
    %done in Topspin program
params.zf = 16; %factor by which data will be zerofilled
params.filter = 'exponential'; %type of FID weighting (current options: 
    %['exponential','gaussian']
params.ap = 100; %spectral apodization (Hz)
params.edge = 40; %Gaussian filter attenuation at ends of FID (dB)
% params.zeropts = 200; %number of points to set to zero at beginning + end of filter
params.ppmwdw=8; %adjusts ppm window of z-spectra: spectral values kept  
    %within +/-ppmwdw
procflgs.norm = true; %if true, performs spectral normalization - this 
    %helps reduce signal fluctuations between B1 values, due to stability 
    %issues e.g.
params.normtype='ppmval'; %type of normalization to use: if 'nosatfit', a fitting
    %routine is performed on the non-saturated spectra to determine scaling 
    %parameters, which are then performed on spectra segregated into 2 
    %groups using k-means. If 'ppmval', raw spectral data are normalized 
    %to a specified ppm signal value (params.ppmnorm)    
params.ppmnorm=12; %ppm value to use for normalization (if params.normtype 
    %set to 'ppmval')
procflgs.peakfit=true; %if true, perform z-peak fitting and use fit 
    %amplitudes instead of MTRasym(ppm value selected) for QUESP
params.pools={'water','OH','amine','amide','Trp'}; %names of pools for fitting  
procflgs.water1st=true; %if true, water (and NOE, if specified) fitted first 
    %using negative ppm values, then all other peaks
params.peaktype='Pseudo-Voigt'; %type of peaks to fit z-spectral lines to
procflgs.PVcharconstr=true; %if true, all fitted Pseudo-Voigt peaks will 
    %have the same Lorentzian-Gaussian character (except water and MT)
procflgs.ppmwt=false; %if true, focus peak fitting on a certain ppm value
params.ppmwt='5.35'; %ppm value(s) to focus peak fitting on, as string
procflgs.fix=true; %if true, peak offsets at specified power index used 
    %for fitting for all other powers
params.fixind=2; %index in power list of spectrum to use for fixing 
    %fitted peak offsets 
procflgs.MTsuperLorentz=true; %if true, MT pool will be fitted using a 
    %super-Lorentzian lineshape
end
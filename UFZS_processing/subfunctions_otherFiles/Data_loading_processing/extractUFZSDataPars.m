% extractUFZSDataPars: Takes in acquisition/processing parameters extracted
% from parameter files and organizes them into output variables to be used
% later on in UFZS processing
%
%   INPUTS:
%       results     -   Struct containing processed results
%       datadirs    -   Cell array of strings specifying the 1D data
%                       directories in the main study directory to be
%                       processed (used only if 1D processing)
%       cfg         -   Struct containing subfields describing user
%                       configuration settings: paths to load/save 
%                       folders, file/parameter identifiers, etc.
%       procflgs    -   Struct containing values of logical flags for
%                       data processing options
%       params      -   Struct containing numerical/string processing 
%                       parameters 
%       pars        -   Struct containing acquisition/processing
%                       parameter values extracted from associated
%                       parameter files
%       np          -   Value of number of points that each spectrum 
%                       contains along the frequency axis
%
%   OUTPUTS:
%       results     -   Struct containing processed results, here updated
%                       to include plotting/legend labels and axes values
%                       as well as the saturation amplitudes in units of uT
%                       and Hz
%       timing      -   Struct containing pulse duration and recovery delay
%                       for UFZS acquisition (both in units of s)
%       nosatidx    -   Vector containing the indices of the spectra
%                       corresponding with zero saturation amplitude. These
%                       correspond with the first dimension of
%                       results.spec.
%
function [results,timing,nosatidx]=extractUFZSDataPars(results,datadirs,cfg,...
    procflgs,params,pars,np)
if procflgs.jeol
    sw_Hz=pars.x_X_SWEEP;
    omega_0_MHz=pars.x_X_FREQ/1e6;
    sw_ppm=sw_Hz/omega_0_MHz;
else
    sw_ppm=pars.sw_p; %spectral width, in ppm
    omega_0_MHz=pars.sfo1; %1H Larmor frequency, in MHz
end
results.specppm=linspace(sw_ppm/2,-sw_ppm/2,np);

%Get saturation dB values for data labels
if procflgs.proc1dflg %1D: pull in values from each dataset
    satdB=zeros(length(datadirs),1);
    for i=1:length(datadirs)
        dBtemp=pars(i).plw;
        satdB(i)=-log10(dBtemp(cfg.plval+1))*10;
        if isinf(satdB(i)) %catch very low powers
            satdB(i)=1000;
        end
        excdB=-log10(dBtemp(2))*10; %PLDB1, to use for nutation freq calc
    end
    nosatidx=find(satdB==1000); %get indices of spectra w/o saturation
elseif procflgs.jeol %2D, JEOL: pull in Hz values from parameters
    satHz=pars.x_SAT_B1;
    nosatidx=find(satHz==0); %get indices of spectra w/o saturation
else %2D, Bruker: pull in dB values from valist
    satdB=pars.valist(~isnan(pars.valist)); %remove NaN in 1st index (dB or W)
    dBtemp=pars.plw;
    excdB=-log10(dBtemp(2))*10; %PLDB1, to use for nutation freq calc 
    nosatidx=find(satdB==1000); %get indices of spectra w/o saturation
end

%Calculate nutation frequencies based upon 90deg pulse calibration (p1,
%pl1), or inputted pw90 value
disp('Calculating saturation amplitudes (in Hz) using P1 and PLDB1 parameters...')
if procflgs.override && ~isinf(params.pw90) && ~isnan(params.pw90)
    pw90 = params.pw90 * 1e-6; %90deg pulse width (s)
elseif procflgs.jeol
    pw90 = pars(1).x_H1_90REF;
else
    pw90 = pars(1).p(2) * 1e-6; %90deg pulse width (s)
end

if procflgs.jeol
    pwExp = pars(1).x_H1_90REF;
    satHz = satHz * pwExp/pw90; 
else
    B1_90 = 1 / 4 / pw90; %B1 nutation frequency corresponding with PLW1 (Hz)
    satHz = B1_90 * 10 .^ ((excdB - satdB) / 10 / 2); %B1 saturation frequencies (Hz)
end
satT = satHz' ./ gamma_;% * 1e6; %B1 saturation amplitudes (uT)
results.speclabels = cell(length(satT),1);
for i = 1:length(satHz)
    results.speclabels{i} = [num2str(satT(i),'%2.1f') ' \muT'];
end

results.satHz=satHz(satHz>1e-9); %remove 0 Hz entries 
results.satT=satT(satT>1e-9); %remove 0 uT entries

if procflgs.jeol
    timing.tp=pars(1).x_SAT_DELAY; %pulse duration for saturation
    timing.rd=pars(1).x_CEST_RELAXATION_DELAY; %recovery delay following pulse sequence    
else
    timing.tp=pars(1).d(cfg.dsatval+1); %pulse duration for saturation
    timing.rd=pars(1).d(cfg.drecval+1); %recovery delay following pulse sequence
    if strcmp(cfg.pptarget1d,'zgpr_llc_iyz') %old pulse seq used d1*2 for recovery
        timing.trec=timing.trec*2;
    end
end 

results.omega_0_MHz=omega_0_MHz;
end
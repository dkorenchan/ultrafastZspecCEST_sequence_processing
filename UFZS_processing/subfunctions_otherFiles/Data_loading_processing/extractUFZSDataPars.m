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

%% PARAMETER READING
% Pull in general NMR acquisition parameters
if procflgs.jeol
    sw_Hz=pars.x_X_SWEEP;
    omega_0_MHz=pars.x_X_FREQ/1e6;
    sw_ppm=sw_Hz/omega_0_MHz;
else
    sw_ppm=pars.sw_p; %spectral width, in ppm
    omega_0_MHz=pars.sfo1; %1H Larmor frequency, in MHz
end
results.specppm=linspace(sw_ppm/2,-sw_ppm/2,np);

if procflgs.procConvflg
    results.satppm=pars.x_SAT_OFF_RWAT;
end

% Pull in timing parameters: saturation pulse duration, recovery delay
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

% Pull in 90deg pulse calibration (for JEOL: also pull in actual "90deg" 
% pulse used)
if procflgs.override && ~isinf(params.pw90) && ~isnan(params.pw90)
    pw90 = params.pw90 * 1e-6; %90deg pulse width (s)
    pwExp = pw90;
elseif procflgs.jeol
    pw90 = pars(1).x_H1_90REF;
    pwExp = pars(1).x_H1_90REF;
else
    pw90 = pars(1).p(2) * 1e-6; %90deg pulse width (s)
    pwExp = pw90;
end

% Pull in saturation amplitude values. If in dB, convert to Hz based upon 
% 90deg pulse calibration (p1, pl1), or inputted pw90 value
if procflgs.proc1dflg %1D: pull in values from each dataset
    satdB=zeros(length(datadirs),1);
    for ii=1:length(datadirs)
        dBtemp=pars(ii).plw;
        satdB(ii)=-log10(dBtemp(cfg.plval+1))*10;
        if isinf(satdB(ii)) %catch very low powers
            satdB(ii)=1000;
        end
        excdB=-log10(dBtemp(2))*10; %PLDB1, to use for nutation freq calc
        disp('Calculating saturation amplitudes (in Hz) using P1 and PLDB1 parameters...')
        B1_90 = 1 / 4 / pw90; %B1 nutation frequency corresponding with PLW1 (Hz)
        satHz = B1_90 * 10 .^ ((excdB - satdB) / 10 / 2); %B1 saturation frequencies (Hz)
    end
%     nosatidx=find(satdB==1000); %get indices of spectra w/o saturation
elseif procflgs.jeol %2D, JEOL: pull in Hz values from parameters
    satHz=pars.x_SAT_B1;
%     nosatidx=find(satHz==0); %get indices of spectra w/o saturation
else %2D, Bruker: pull in dB values from valist
    satdB=pars.valist(~isnan(pars.valist)); %remove NaN in 1st index (dB or W)
    dBtemp=pars.plw;
    excdB=-log10(dBtemp(2))*10; %PLDB1, to use for nutation freq calc 
    disp('Calculating saturation amplitudes (in Hz) using P1 and PLDB1 parameters...')
    B1_90 = 1 / 4 / pw90; %B1 nutation frequency corresponding with PLW1 (Hz)
    satHz = B1_90 * 10 .^ ((excdB - satdB) / 10 / 2); %B1 saturation frequencies (Hz)    
%     nosatidx=find(satdB==1000); %get indices of spectra w/o saturation
end

satHz = satHz * pwExp/pw90; 


%% PARAMETER PARSING
% Find non-saturated indices and pout together labels for spectra
% depending on what was arrayed: B1, or Tsat
if numel(satHz)>1 && numel(timing.tp)==1 %B1 varies
    disp('Saturation pulse amplitudes are arrayed in the experiment(s).')
    nosatidx=find(satHz<1e-3); %get indices of spectra w/ low B1 saturation
    % Use amplitudes in uT and rad/s
    satT = satHz' ./ gamma_;% * 1e6; %B1 saturation amplitudes (uT)
    satRadS = satHz' * 2 * pi;
    results.speclabels = cell(length(satT),1);
    for ii = 1:length(satHz)
        results.speclabels{ii} = [num2str(satT(ii),'%2.1f') ' \muT (' ...
            num2str(satRadS(ii),'%4.0f') ' rad/s)'];
    end
elseif numel(timing.tp)>1 && numel(satHz)==1 %Tsat varies
    disp('Saturation pulse durations are arrayed in the experiment(s).')
    nosatidx=find(timing.tp<1e-3); %get indices of spectra w/ super short saturation
    results.speclabels = cell(length(timing.tp),1);
    for ii = 1:length(timing.tp)
        results.speclabels{ii} = [num2str(timing.tp(ii),'%2.3f') ' s'];
    end
    satT = satHz ./ gamma_;
    satRadS = satHz * 2 * pi;
elseif numel(satHz)==1 && numel(timing.tp)==1 %only 1 B1 and Tsat combination
    disp('Only one saturation pulse amplitude and duration used.')
    nosatidx=[];
    results.speclabels = cell(0);
else
    error(['Cannot process data! Either the pulse amplitude or duration ' ...
        'must be arrayed, not both!'])
end

results.satHz=satHz(satHz>1e-9); %remove 0 Hz entries 
results.satT=satT(satHz>1e-9); %remove 0 Hz entries 
results.satRadS=satRadS(satHz>1e-9); %remove 0 Hz entries 

results.omega_0_MHz=omega_0_MHz;
end
% PROC_UFZS_DATA.m
%
% Takes specified Topspin datasets with ultrafast Z-spectra, plots the 
% spectra overlay, and then calculates MTR asymmetry
%
% INPUTS (within first section of code)
%   prevresults (optional)  -   Structure obtained from previous run; if 
%                               used as input, will be plotted    
%
% OUTPUTS
%   results     -   Structure containing spectra, z-spectra, and MTR 
%                   asymmetry, along with corresponding ppm values and data 
%                   labels   
% 
% REQUIRED TOOLBOXES:       
% Curve Fitting                     Optimization      Image Processing    
% Bioinformatics (for JEOL data)
%
% USAGE (processing mode):
%
%                

function results = PROC_UFZS_DATA(prevresults)
%% FUNCTION INITIALIZATION

[cfg,opts]=initUserSettingsUFZS();

% initializing graphic parameters
set(0, 'DefaultAxesLineWidth',1,'DefaultAxesFontSize',22,...
          'DefaultAxesFontWeight','normal','DefaultAxesFontName','Arial',...
          'DefaultLineLineWidth',2,'DefaultLineMarkerSize',14,...
          'DefaultAxesXColor','black','DefaultAxesYColor','black',...
          'DefaultAxesZColor','black','DefaultAxesGridColor','black',...
          'DefaultFigureColor',[1 1 1],'DefaultAxesTitleFontSizeMultiplier',1.4,...
          'DefaultAxesLabelFontSizeMultiplier',1.2);   

if nargin < 1
disp('No previous data specified. Starting processing pipeline...')


%% INTERACTIVE PARAMETER SETTINGS, DISPLAY NOTIFICATIONS

% Show interactive window for selecting processing parameters
[procflgs,params]=procOptionsGUI;


%% DATA PROCESSING: FOLDER SPECIFICATION
%
% Prompt user to specify folder containing data directories, then
% automatically detect which ones pertain to ultrafast z-spectroscopy
% datasets. Have the user specify which one(s) to use for processing.
[pathname,datadirs,procflgs]=findUFZSdatasets(cfg,procflgs);


%% DATA PROCESSING: LOADING
%
% Load in raw FIDs or data processed in TopSpin and located in pdata
% folder, depending on procflgs.topproc value and whether 1D or 2D processing is
% being performed (proc1dflg)
[results.spec,pars]=Load_Preprocess_Spectra(pathname,datadirs,procflgs,params);

% Parse out extracted data parameters into condensed variables to use later
if procflgs.procConvflg
    np=size(results.spec,3);
else
    np=size(results.spec,2);
end
[results,timing,nosatidx]=extractUFZSDataPars(results,datadirs,cfg,procflgs,params,...
    pars,np);


%% DATA PROCESSING: Z-SPECTRAL NORMALIZATION, CALCULATION
% If procflgs.norm == true, normalize spectra by one of the two procedures
% described to the user in the GUI
normpars=struct;
if procflgs.norm && ~procflgs.procConvflg
    [results,normpars]=normalizeAllSpectra(results,nosatidx,params);
end

% Calculate all z-spectra + MTR asymmetry
results=calcZspecMTRasym(results,procflgs,params,nosatidx,normpars);

% Plot data thus far
zspecPlot(results,procflgs.procConvflg,params.ppmwdw);


%% Z-SPECTRAL FITTING, QUESP ANALYSIS
% Perform z-spectral peak fitting, if indicated by procflgs.peakfit
if numel(timing.tp)>1
    disp(['Cannot currently perform spectral fitting when saturation ' ...
        'pulse duration is arrayed! Skipping...'])
    procflgs.peakfit=false;
end
if procflgs.peakfit
    results=zPeaks_fit_display(results,params,procflgs);
else
    params.pools={'MTRasym'}; %this is used later on
end

% Perform QUESP fitting and plotting
results=QUESP_fit_display(results,params,timing,opts);


else
disp('Previous processed data specified. Skipping to plotting...')
results = prevresults;
zspecPlot(results);
end

disp('...Function finished!')

end


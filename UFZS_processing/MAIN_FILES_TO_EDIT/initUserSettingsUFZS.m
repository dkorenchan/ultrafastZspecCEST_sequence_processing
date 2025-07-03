% initUserSettings: Sets fields of the struct settings to include
% information pertaining to the user's computer as well as processing
% preferences
%
%   INPUTS:     None
%
%   OUTPUTS:    
%       configs     -   Struct containing subfields describing user
%                       configuration settings: paths to load/save 
%                       folders, file/parameter identifiers, etc.
%       opts        -   Struct containing user start, upper, and lower
%                       bounds for QUESP fitting
%
function [configs,opts]=initUserSettingsUFZS()
disp('Loading user configuration settings from file initUserSettingsUFZS.m...')
%% USER CONFIGURATION SETTINGS
% Paths and strings to identify datasets to load
configs.load_path = '/Users/dk384/Documents/Laboratory/MGH_Farrar-Rosen/Data/14T_Bruker/NMR'; 
    %path to data
configs.pptarget1d = 'ufzs1d'; %string to use for finding 1D Z-spectral datasets (DK 
    %pulse sequence file - set plval = 9, dsatval = 9)
configs.pptarget2d = 'ufz'; %string to use for finding 2D Z-spectral datasets (DK 
    %pulse sequence file - set plval = 9, dsatval = 9, drecval = 1)

configs.pptarget2djeol = 'ufz'; %string to use for finding 2D UFZ JEOL datasets
configs.pptarget3djeol = 'cest'; %string to use for finding JEOL CEST datasets
    %(NOTE: 3D datasets are currently saved as separate 2D datasets;
    % if a conventional z-spectroscopy dataset (not ultrafast), then 
    % the file name will match .pptarget3djeol but NOT .pptarget2djeol!)

% Values for identifying parameters to load in, based upon pulse sequence
% file
configs.plval=9; %index of PLDB parameter array used for saturation power
    %(ex. if PL9 used in pulse program, set plval=9)
configs.dsatval=9; %index of D parameter array used for saturation duration
    %(ex. if D9 used in pulse program, set dsatval=9)
configs.drecval=1; %index of D parameter array used for relaxation delay
    %(ex. if D1 used in pulse program, set drecval=1)
% sgnlcutoff=20; %signal ratio (relative to max), used to cut out 
%     %low-signal regions of z-spectra (lower value, stricter cutoff)

% % Values for old sequence
% pptarget1d = 'zgpr_llc_iyz'; %string to use for finding 1D Z-spectral datasets
%     %(previous pulse sequence file - set plval = 9, dsatval = 1, drecval = 1)
% pptarget1d='zgpr_llc_iyz_fast_mse';
% plval=9;
% dsatval=1;
% drecval=1
% ppmwdw=6; 

% for QUESP fitting
% opts - Fit start value and boundaries
% fs (M)| ksw (Hz)
estConc = 11 * 1.62/(110e3);
% estConc = 3 * 50/(110e3);
opts.Lower      =   [estConc*.0001  50];
opts.StartPoint =   [estConc*.1     1000];
opts.Upper      =   [estConc*100    10000];
% % Previously used for Trp peptide fitting:
% estConc = 11 * 1.62/(110e3);
% opts.Lower      =   [estConc*.0001  100];
% opts.StartPoint =   [estConc*.01    1000];
% opts.Upper      =   [estConc*25     30000];
end
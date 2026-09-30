% procOptionsGUI: Creates a graphical user interface where the user can
% specify processing options for z-spectral data, fitting, etc.
%
%   INPUTS:     NONE
%
%   OUTPUTS:    
%       procflgs    -   Struct containing values of logical flags for
%                       data processing options
%       params      -   Struct containing numerical/string processing 
%                       parameters
%
function [procflgs,params]=procOptionsGUI
% Initialize default settings for GUI to display when loaded
[procflgs,params]=initializeGUIdefaults;
flgnames = fieldnames(procflgs);

paramfig = uifigure('Name','Ultrafast z-spectroscopy: processing parameters');
paramfig.Position(3) = 700;
paramfig.Position(4) = 700;

uicheckbox(paramfig,'Position',[20 605 800 20],'Tag','override',...
    'Value',procflgs.override,'ValueChangedFcn',@toggleFlag,...
    'Text','Override a parameter read from the data files (set to Inf to use default)');
% Items in this buttongroup disappear/reappear based upon procflgs.override
ui.override.bgrp = uibuttongroup(paramfig,'Position',[40 575 620 30],...
    'Visible',procflgs.override,'SelectionChangedFcn',@changeButton);
ui.override.pw90set = uieditfield(ui.override.bgrp,'numeric','Position',[10 5 50 20],...
    'Tag','pw90','Value',params.pw90,'ValueChangedFcn',@setParam);
ui.override.pw90lab = uilabel(ui.override.bgrp,'Position',[65 5 500 20],...
    'Text','90 degree pulse width (us)');

uicheckbox(paramfig,'Position',[20 550 800 20],'Tag','topproc',...
    'Value',procflgs.topproc,'ValueChangedFcn',@toggleFlag,...
    'Text','Load in z-spectra already processsed in TopSpin');
% Items in this buttongroup disappear/reappear based upon procflgs.topproc
ui.topproc.bgrp = uibuttongroup(paramfig,'Position',[40 470 620 80],...
    'Visible',~procflgs.topproc,'SelectionChangedFcn',@changeButton);
ui.topproc.zfset = uieditfield(ui.topproc.bgrp,'numeric','Position',[10 50 50 20],'Tag','zf',...
    'Value',params.zf,'ValueChangedFcn',@setParam);
ui.topproc.zflab = uilabel(ui.topproc.bgrp,'Position',[65 50 500 20],...
    'Text','Zerofilling factor (integer value >0)');
ui.topproc.filterlab = uilabel(ui.topproc.bgrp,'Position',[10 30 100 20],...
    'Text','FID filter type:');
ui.topproc.expset = uiradiobutton(ui.topproc.bgrp,'Position',[110 30 100 20],...
    'Value',strcmp(params.filter,'exponential'),'Text','Exponential',...
    'Parent',ui.topproc.bgrp);
ui.topproc.exp.apset = uieditfield(ui.topproc.bgrp,'numeric','Tag','ap',...
    'Position',[110 10 50 20],'Visible',strcmp(params.filter,'exponential'),...
    'Value',params.ap,'ValueChangedFcn',@setParam);
ui.topproc.exp.aplab = uilabel(ui.topproc.bgrp,'Position',[165 10 600 20],...
    'Visible',strcmp(params.filter,'gaussian'),...
    'Text','Exponential decay constant (Hz)','Visible',strcmp(params.filter,'exponential'));
ui.topproc.gaussset = uiradiobutton(ui.topproc.bgrp,'Position',[220 30 100 20],...
    'Value',strcmp(params.filter,'gaussian'),'Text','Gaussian',...
    'Parent',ui.topproc.bgrp);
ui.topproc.gauss.edgeset = uieditfield(ui.topproc.bgrp,'numeric','Tag','edge',...
    'Position',[110 10 50 20],'Visible',strcmp(params.filter,'gaussian'),...
    'Value',params.edge,'ValueChangedFcn',@setParam);
ui.topproc.gauss.edgelab = uilabel(ui.topproc.bgrp,'Position',[160 10 600 20],...
    'Visible',strcmp(params.filter,'gaussian'),...
    'Text','Gaussian filter attenuation at FID endpoints (dB)');
ui.topproc.rectchk = uicheckbox(ui.topproc.bgrp,'Position',[330 50 130 20],...
    'Tag','rect','Value',procflgs.rect,'ValueChangedFcn',@toggleFlag,...
    'Text','Rect window:');
ui.topproc.rectpctset = uieditfield(ui.topproc.bgrp,'numeric','Tag','rectpct',...
    'Position',[460 50 40 20],'Visible',procflgs.rect,...
    'Limits',[0 50],'Value',params.rectpct,'ValueChangedFcn',@setParam);
ui.topproc.rectpctlab = uilabel(ui.topproc.bgrp,'Position',[505 50 110 20],...
    'Visible',procflgs.rect,'Text','% zeroed each end');

uilabel(paramfig,'Position',[20 442 800 20],...
    'Text','Window size for displaying z-spectra (+/- ppm):');
uieditfield(paramfig,'numeric','Position',[280 442 40 20],...
    'Tag','ppmwdw','Value',params.ppmwdw,'ValueChangedFcn',@setParam);

uicheckbox(paramfig,'Position',[20 420 800 20],'Tag','norm',...
    'Value',procflgs.norm,'ValueChangedFcn',@toggleFlag,...
    'Text','Normalize z-spectra across all saturation B1 values');
% Items in this buttongroup disappear/reappear based upon procflgs.norm
ui.norm.bgrp = uibuttongroup(paramfig,'Position',[40 360 620 60],...
    'Visible',procflgs.norm,'SelectionChangedFcn',@changeButton);
ui.norm.nosatset = uiradiobutton(ui.norm.bgrp,'Position',[10 30 600 20],...
    'Value',strcmp(params.normtype,'nosatfit'),'Parent',ui.norm.bgrp,...
    'Text',['Use non-saturated spectrum scaling + clustering to scale '...
    'raw spectra + calculate z-spectra']);
ui.norm.ppmvalset = uiradiobutton(ui.norm.bgrp,'Position',[10 10 400 20],...
    'Value',strcmp(params.normtype,'ppmval'),'Parent',ui.norm.bgrp,...
    'Text','Scale all raw spectra to be equal at this ppm value:');
ui.norm.pnset = uieditfield(ui.norm.bgrp,'numeric','Position',[310 10 40 20],...
    'Tag','ppmnorm','Value',params.ppmnorm,'ValueChangedFcn',@setParam,...
    'Enable',strcmp(params.normtype,'ppmval'));

uicheckbox(paramfig,'Position',[20 330 800 20],'Tag','peakfit',...
    'Value',procflgs.peakfit,'ValueChangedFcn',@toggleFlag,...
    'Text','Use multi-pool lineshape fitting, rather than MTR asymmetry, for QUESP analysis');
% Items in this buttongroup disappear/reappear based upon procflgs.peakfit
ui.peakfit.bgrp = uibuttongroup(paramfig,'Position',[40 50 620 280],...
    'Visible',procflgs.peakfit,'SelectionChangedFcn',@changeButton);
ui.peakfit.poollab = uilabel(ui.peakfit.bgrp,'Position',[10 250 400 20],...
    'Text','Pools to fit (water always included):');
ui.peakfit.poolamide = uicheckbox(ui.peakfit.bgrp,'Position',[220 250 100 20],...
    'Tag','amide','Value',sum(strcmp(params.pools,'amide')),...
    'ValueChangedFcn',@setPools,'Text','amide');
ui.peakfit.poolamine = uicheckbox(ui.peakfit.bgrp,'Position',[320 250 100 20],...
    'Tag','amine','Value',sum(strcmp(params.pools,'amine')),...
    'ValueChangedFcn',@setPools,'Text','amine');
ui.peakfit.poolOH = uicheckbox(ui.peakfit.bgrp,'Position',[420 250 100 20],...
    'Tag','OH','Value',sum(strcmp(params.pools,'OH')),...
    'ValueChangedFcn',@setPools,'Text','OH');
ui.peakfit.poolTrp = uicheckbox(ui.peakfit.bgrp,'Position',[220 230 100 20],...
    'Tag','Trp','Value',sum(strcmp(params.pools,'Trp')),...
    'ValueChangedFcn',@setPools,'Text','Trp');
ui.peakfit.poolppm4pt4 = uicheckbox(ui.peakfit.bgrp,'Position',[320 230 100 20],...
    'Tag','ppm4pt4','Value',sum(strcmp(params.pools,'ppm4pt4')),...
    'ValueChangedFcn',@setPools,'Text','4.4 ppm');
ui.peakfit.poolppm7pt3 = uicheckbox(ui.peakfit.bgrp,'Position',[420 230 100 20],...
    'Tag','ppm7pt3','Value',sum(strcmp(params.pools,'ppm7pt3')),...
    'ValueChangedFcn',@setPools,'Text','7.3 ppm');
ui.peakfit.poolppm9pt8 = uicheckbox(ui.peakfit.bgrp,'Position',[220 210 100 20],...
    'Tag','ppm9pt8','Value',sum(strcmp(params.pools,'ppm9pt8')),...
    'ValueChangedFcn',@setPools,'Text','9.8 ppm');
ui.peakfit.poolNOE = uicheckbox(ui.peakfit.bgrp,'Position',[320 210 100 20],...
    'Tag','NOE','Value',sum(strcmp(params.pools,'NOE')),...
    'ValueChangedFcn',@setPools,'Text','rNOE');
ui.peakfit.poolMT = uicheckbox(ui.peakfit.bgrp,'Position',[420 210 100 20],...
    'Tag','MT','Value',sum(strcmp(params.pools,'MT')),...
    'ValueChangedFcn',@setPools,'Text','MT');
ui.peakfit.peaktypelab = uilabel(ui.peakfit.bgrp,'Position',[10 190 80 20],...
    'Text','Peak fit type:');
ui.peakfit.Lorentzset = uiradiobutton(ui.peakfit.bgrp,'Position',[210 190 90 20],...
    'Value',strcmp(params.peaktype,'Lorentzian'),'Text','Lorentzian',...
    'Parent',ui.peakfit.bgrp);
ui.peakfit.PseudoVoigtset = uiradiobutton(ui.peakfit.bgrp,'Position',[100 190 90 20],...
    'Value',strcmp(params.peaktype,'Pseudo-Voigt'),'Text','Pseudo-Voigt',...
    'Parent',ui.peakfit.bgrp);
ui.peakfit.Rexset = uiradiobutton(ui.peakfit.bgrp,'Position',[305 190 55 20],...
    'Value',strcmp(params.peaktype,'Rex'),'Text','Rex',...
    'Parent',ui.peakfit.bgrp);
ui.peakfit.PVcharconstr = uicheckbox(ui.peakfit.bgrp,'Position',[100 170 450 20],...
    'Visible',strcmp(params.peaktype,'Pseudo-Voigt'),'Tag','PVcharconstr',...
    'Value',procflgs.PVcharconstr,'ValueChangedFcn',@toggleFlag,...
    'Text','Keep Gaussian-Lorentzian character same across all peaks (except water + MT)');
ui.invZfit.bgrp = uibuttongroup(ui.peakfit.bgrp,'Position',[30 3 590 22],...
    'Visible',~procflgs.invZfit,'BorderType','none');
ui.peakfit.fixset = uicheckbox(ui.peakfit.bgrp,'Position',[10 150 600 20],'Tag','fix',...
    'Value',procflgs.fix,'ValueChangedFcn',@toggleFlag,...
    'Text','Fit lineshape offsets for selected z-spectrum only, then use for all other spectra');
ui.peakfit.fixind = uieditfield(ui.peakfit.bgrp,'numeric',...
    'Position',[30 130 50 20],'Tag','fixind','Visible',procflgs.fix,...
    'Value',params.fixind,'ValueChangedFcn',@setParam);
ui.peakfit.fixlab = uilabel(ui.peakfit.bgrp,'Position',[90 130 400 20],...
    'Visible',procflgs.fix,'Text','Spectrum index to use to fix peak offsets');
ui.peakfit.ppmwt = uicheckbox(ui.peakfit.bgrp,'Position',[10 110 600 20],'Tag','ppmwt',...
    'Value',procflgs.ppmwt,'ValueChangedFcn',@toggleFlag,...
    'Text','Focus fitting on a specific ppm value (e.g. for a small peak)');
ui.peakfit.ppmwtset = uieditfield(ui.peakfit.bgrp,'text',...
    'Position',[30 80 90 20],'Tag','ppmwt','Visible',procflgs.ppmwt,...
    'Value',params.ppmwt,'ValueChangedFcn',@setParam);
ui.peakfit.ppmwtlab = uilabel(ui.peakfit.bgrp,'Position',[130 70 350 40],...
    'Visible',procflgs.ppmwt,'WordWrap','on', ...
    'Text',['Spectral ppm value(s) to focus on while fitting, separated by '...
    'commas or semicolons']);
ui.peakfit.MTlineshapelab = uilabel(ui.peakfit.bgrp,'Position',[10 53 145 16],...
    'Text','MT pool lineshape:','Enable',logical(sum(strcmp(params.pools,'MT'))));
ui.peakfit.MTlineshapebgrp = uibuttongroup(ui.peakfit.bgrp,'Position',[160 48 405 22],...
    'Enable',logical(sum(strcmp(params.pools,'MT'))),...
    'SelectionChangedFcn',@changeButton,'BorderType','none');
ui.peakfit.MTlorentz = uiradiobutton(ui.peakfit.MTlineshapebgrp,'Position',[5 1 105 20],...
    'Value',strcmp(params.MTlineshape,'lorentzian'),'Text','Lorentzian');
ui.peakfit.MTgauss = uiradiobutton(ui.peakfit.MTlineshapebgrp,'Position',[115 1 105 20],...
    'Value',strcmp(params.MTlineshape,'gaussian'),'Text','Gaussian');
ui.peakfit.MTsuperLorentz = uiradiobutton(ui.peakfit.MTlineshapebgrp,'Position',[225 1 175 20],...
    'Value',strcmp(params.MTlineshape,'superlorentzian'),'Text','Super-Lorentzian');
ui.peakfit.invZfit = uicheckbox(ui.peakfit.bgrp,'Position',[10 26 540 20],...
    'Tag','invZfit','Value',~procflgs.invZfit,'ValueChangedFcn',@toggleFlag,...
    'Text','Use conventional Z fitting instead of 1/Z (NOT RECOMMENDED)');
ui.peakfit.water1st = uicheckbox(ui.invZfit.bgrp,'Position',[2 2 580 20],'Tag','water1st',...
    'Value',procflgs.water1st,'ValueChangedFcn',@toggleFlag,...
    'Text','Fit water (and rNOE/MT, if specified) first, using negative ppm values');

uibutton(paramfig,'Position',[20 10 140 20],'Text','Continue',...
    'ButtonPushedFcn',@endParamSet);
uibutton(paramfig,'Position',[170 10 220 20],'Text','Restore default settings',...
    'ButtonPushedFcn',@restoreDefaults);
uicheckbox(paramfig,'Position',[400 10 390 20],'Tag','saveDefaults',...
    'Value',procflgs.saveDefaults,'ValueChangedFcn',@toggleFlag,...
    'Text','Save current selections as defaults for next run');
waitfor(paramfig);

% Save as defaults if requested
if procflgs.saveDefaults
    defaultsFile = fullfile(fileparts(mfilename('fullpath')), 'UFZS_GUI_defaults.mat');
    procflgs_saved = rmfield(procflgs, 'saveDefaults');
    params_saved   = params;
    save(defaultsFile, 'procflgs_saved', 'params_saved');
    disp('GUI selections saved as defaults for next run.')
end

% Display parameter settings based on GUI selections
displayGUIselections(procflgs,params);

% Adjust relevant parameters based upon settings
if ~procflgs.ppmwt %no ppm value to focus fitting on
    params.ppmwt = NaN;
else
    params.ppmwt = str2num(params.ppmwt); %convert to numeric array
end


%% INTERNAL CALLBACK FUNCTIONS FOR UI ELEMENTS
%
% toggleFlag: changes value of processing flag + interactive ability of 
% other parameter settings based upon checkbox value
function toggleFlag(src,~)
% Change flag value
procflgs.(src.Tag) = ~procflgs.(src.Tag);
% Change visibility of button groups
for iii = 1:numel(flgnames)
    fnam = flgnames{iii};
    if isfield(ui,fnam)
        if strcmp(fnam,'topproc') || strcmp(fnam,'invZfit')
            set(ui.(fnam).bgrp,'Visible',~procflgs.(fnam));
        else
            set(ui.(fnam).bgrp,'Visible',procflgs.(fnam));
        end
    end
end
% Change visibility of other flags
set(ui.topproc.rectpctset,'Visible',procflgs.rect);
set(ui.topproc.rectpctlab,'Visible',procflgs.rect);
set(ui.peakfit.fixind,'Visible',procflgs.fix);
set(ui.peakfit.fixlab,'Visible',procflgs.fix);
set(ui.peakfit.ppmwtset,'Visible',procflgs.ppmwt);
set(ui.peakfit.ppmwtlab,'Visible',procflgs.ppmwt);
end


% setParam: updates parameter based upon value input in edit field
function setParam(src,~)
params.(src.Tag) = src.Value;    
end


% changeButton: updates button group if new radio button selected
function changeButton(~,~)
% Change value of params.filter based on radio button selection
if ui.topproc.expset.Value
    params.filter='exponential';
elseif ui.topproc.gaussset.Value
    params.filter='gaussian';
end
% Change value of params.normtype based on radio button selection
if ui.norm.nosatset.Value
    params.normtype='nosatfit';
elseif ui.norm.ppmvalset.Value
    params.normtype='ppmval';
end
% Change value of params.peaktype based on radio button selection
if ui.peakfit.Lorentzset.Value
    params.peaktype='Lorentzian';
elseif ui.peakfit.PseudoVoigtset.Value
    params.peaktype='Pseudo-Voigt';
elseif ui.peakfit.Rexset.Value
    params.peaktype='Rex';
end
% Update MT lineshape selection and keep MTsuperLorentz flag in sync
if ui.peakfit.MTlorentz.Value
    params.MTlineshape = 'lorentzian';       procflgs.MTsuperLorentz = false;
elseif ui.peakfit.MTgauss.Value
    params.MTlineshape = 'gaussian';         procflgs.MTsuperLorentz = false;
elseif ui.peakfit.MTsuperLorentz.Value
    params.MTlineshape = 'superlorentzian';  procflgs.MTsuperLorentz = true;
end
% Update visibility of auxiliary parameters
tpfields = {'exp','gauss'};
for iii = 1:numel(tpfields)
    tpf = tpfields{iii};
    tpsubfields = fieldnames(ui.topproc.(tpf));
    for jjj = 1:numel(tpsubfields)
        tpsf = tpsubfields{jjj};
        if strcmp(tpf,'exp')
            set(ui.topproc.(tpf).(tpsf),'Visible',strcmp(params.filter,'exponential'))
        elseif strcmp(tpf,'gauss')
            set(ui.topproc.(tpf).(tpsf),'Visible',strcmp(params.filter,'gaussian'))
        end
    end
end
set(ui.norm.pnset,'Enable',strcmp(params.normtype,'ppmval'));
set(ui.peakfit.PVcharconstr,'Visible',strcmp(params.peaktype,'Pseudo-Voigt'));
end


% setPools: adjusts params.pools to include pools specified in GUI
function setPools(src,~)
% See if pool name already included in params.pools - if there, remove; if
% not there, add
if sum(strcmp(params.pools,src.Tag)) > 0
    params.pools = params.pools(~strcmp(params.pools,src.Tag));
else
    params.pools = [params.pools {src.Tag}];
end
% Also, verify that water is included!
if sum(strcmp(params.pools,'water')) < 1
    params.pools = [{'water'} params.pools];
end
% Reorder based upon basePools, which contains the correct order. Go
% through basePools, find each entry in params.pools, then move to end of 
% params.pools. The result will be params.pools ordered the same as
% basePools! 
basePools={'water','NOE','MT','OH','amine','amide','ppm4pt4','Trp','ppm7pt3',...
    'ppm9pt8'};
for iii=1:numel(basePools)
    params.pools=[params.pools(~strcmp(params.pools,basePools{iii})), ...
        params.pools(strcmp(params.pools,basePools{iii}))];
end

% Change visibility of other flags
mtEnabled = sum(strcmp(params.pools,'MT')) > 0;
set(ui.peakfit.MTlineshapebgrp,'Enable',mtEnabled);
set(ui.peakfit.MTlineshapelab,'Enable',mtEnabled);
end


% endParamSet: closes parameter setting panel, continues with function
function endParamSet(~,~)
close(paramfig)
end


% restoreDefaults: resets all GUI elements to the hardcoded defaults
function restoreDefaults(~,~)
    [procflgs, params] = initializeGUIdefaults(true);
    set(findobj(paramfig,'Tag','override'), 'Value', procflgs.override);
    set(ui.override.bgrp,  'Visible', procflgs.override);
    set(ui.override.pw90set,'Value', params.pw90);
    set(findobj(paramfig,'Tag','topproc'), 'Value', procflgs.topproc);
    set(ui.topproc.bgrp, 'Visible', ~procflgs.topproc);
    set(ui.topproc.zfset,'Value', params.zf);
    set(ui.topproc.expset,   'Value', strcmp(params.filter,'exponential'));
    set(ui.topproc.gaussset, 'Value', strcmp(params.filter,'gaussian'));
    set(ui.topproc.exp.apset,    'Value', params.ap,   'Visible', strcmp(params.filter,'exponential'));
    set(ui.topproc.exp.aplab,    'Visible', strcmp(params.filter,'exponential'));
    set(ui.topproc.gauss.edgeset,'Value', params.edge, 'Visible', strcmp(params.filter,'gaussian'));
    set(ui.topproc.gauss.edgelab,'Visible', strcmp(params.filter,'gaussian'));
    set(ui.topproc.rectchk,      'Value', procflgs.rect);
    set(ui.topproc.rectpctset,   'Value', params.rectpct, 'Visible', procflgs.rect);
    set(ui.topproc.rectpctlab,   'Visible', procflgs.rect);
    set(findobj(paramfig,'Tag','ppmwdw'), 'Value', params.ppmwdw);
    set(findobj(paramfig,'Tag','norm'), 'Value', procflgs.norm);
    set(ui.norm.bgrp, 'Visible', procflgs.norm);
    set(ui.norm.nosatset,  'Value', strcmp(params.normtype,'nosatfit'));
    set(ui.norm.ppmvalset, 'Value', strcmp(params.normtype,'ppmval'));
    set(ui.norm.pnset, 'Value', params.ppmnorm, 'Enable', strcmp(params.normtype,'ppmval'));
    set(findobj(paramfig,'Tag','peakfit'), 'Value', procflgs.peakfit);
    set(ui.peakfit.bgrp, 'Visible', procflgs.peakfit);
    set(ui.peakfit.poolamide,  'Value', sum(strcmp(params.pools,'amide')));
    set(ui.peakfit.poolamine,  'Value', sum(strcmp(params.pools,'amine')));
    set(ui.peakfit.poolOH,     'Value', sum(strcmp(params.pools,'OH')));
    set(ui.peakfit.poolTrp,    'Value', sum(strcmp(params.pools,'Trp')));
    set(ui.peakfit.poolppm4pt4,'Value', sum(strcmp(params.pools,'ppm4pt4')));
    set(ui.peakfit.poolppm7pt3,'Value', sum(strcmp(params.pools,'ppm7pt3')));
    set(ui.peakfit.poolppm9pt8,'Value', sum(strcmp(params.pools,'ppm9pt8')));
    set(ui.peakfit.poolNOE,    'Value', sum(strcmp(params.pools,'NOE')));
    set(ui.peakfit.poolMT,     'Value', sum(strcmp(params.pools,'MT')));
    set(ui.peakfit.Lorentzset,    'Value', strcmp(params.peaktype,'Lorentzian'));
    set(ui.peakfit.PseudoVoigtset,'Value', strcmp(params.peaktype,'Pseudo-Voigt'));
    set(ui.peakfit.Rexset,        'Value', strcmp(params.peaktype,'Rex'));
    set(ui.peakfit.PVcharconstr,  'Value', procflgs.PVcharconstr, ...
        'Visible', strcmp(params.peaktype,'Pseudo-Voigt'));
    set(ui.peakfit.fixset,'Value', procflgs.fix);
    set(ui.peakfit.fixind,'Value', params.fixind, 'Visible', procflgs.fix);
    set(ui.peakfit.fixlab,'Visible', procflgs.fix);
    set(ui.peakfit.ppmwt,   'Value', procflgs.ppmwt);
    set(ui.peakfit.ppmwtset,'Value', params.ppmwt, 'Visible', procflgs.ppmwt);
    set(ui.peakfit.ppmwtlab,'Visible', procflgs.ppmwt);
    mtEnabled = logical(sum(strcmp(params.pools,'MT')));
    set(ui.peakfit.MTlineshapebgrp,'Enable', mtEnabled);
    set(ui.peakfit.MTlineshapelab, 'Enable', mtEnabled);
    set(ui.peakfit.MTlorentz,     'Value', strcmp(params.MTlineshape,'lorentzian'));
    set(ui.peakfit.MTgauss,       'Value', strcmp(params.MTlineshape,'gaussian'));
    set(ui.peakfit.MTsuperLorentz,'Value', strcmp(params.MTlineshape,'superlorentzian'));
    set(ui.peakfit.invZfit, 'Value', ~procflgs.invZfit);
    set(ui.invZfit.bgrp,    'Visible', ~procflgs.invZfit);
    set(ui.peakfit.water1st,'Value', procflgs.water1st);
    set(findobj(paramfig,'Tag','saveDefaults'),'Value', procflgs.saveDefaults);
end


end
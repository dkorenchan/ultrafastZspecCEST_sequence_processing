% findUFZSdatasets: Prompts the user to specify the main directory
% containing multiple NMR datasets, then searches for ultrafast
% z-spectroscopy datasets and prompts user to select the dataset(s) to
% process.
%
%   INPUTS:
%       cfg         -   Struct containing subfields describing user
%                       configuration settings: paths to load/save folders, 
%                       file/parameter identifiers, etc.
%       procflgs    -   Struct containing values of logical flags for
%                       data processing options
%
%   OUTPUTS:
%       pathname    -   String containing path to study directory
%       dirs        -   Cell array of strings specifying the data
%                       directory/ies in the main study directory to be
%                       processed
%       procflgs    -   Struct containing values of logical flags for
%                       data processing options, updated by this function
%
function [pathname,dirs,pflgs]=findUFZSdatasets(cfg,pflgs)
% First, have the user specify the study directory to use
[pathname,dirfind]=loadDirectories(cfg.load_path);

dirs1d=cell(length(dirfind),1);
dirs2d=dirs1d;
dirsConv=dirs1d;
ctr1d=1;
ctr2d=1;

% Identify the data directories containing ultrafast Z-spectra based upon 
% the PULPROG value in acqus, segregate into 1D or 2D data accordingly
for i = 1:length(dirfind)
    if ~isnan(str2double(dirfind{i})) %only take directories that are a number
        ppnam=readParsTopspin(fullfile(pathname,dirfind{i}),'acqus',...
            {'##$PULPROG'});
        if contains(ppnam,cfg.pptarget1d,'IgnoreCase',true)
            dirs1d(ctr1d)=dirfind(i);
            ctr1d=ctr1d+1;
        elseif contains(ppnam,cfg.pptarget2d,'IgnoreCase',true)
            dirs2d(ctr2d)=dirfind(i);
            ctr2d=ctr2d+1;        
        end
    end
end
dirs1d=dirs1d(1:ctr1d-1); %remove extra entries
dirs2d=dirs2d(1:ctr2d-1); %remove extra entries

% See if JEOL datasets exist (as .jdx files) within the study directory, 
% and if so add them either to dirs2d (if UFZS datasets) or to dirsConv (if
% conventional ZS datasets)
jeoldirs = dir(fullfile(pathname,'**',strcat('*',cfg.pptarget3djeol,'*.jdx'))); 
    %this searches all subdirectories
if ~isempty(jeoldirs)
    disp('JEOL datasets discovered! Adding to list...')
    for i=1:numel(jeoldirs)
        jdirnam=extractAfter(jeoldirs(i).folder,pathname);
        jdirfile=jeoldirs(i).name;
        % Divide into 2D UFZS datasets vs 3D conventional ZS datasets
        if contains(jdirfile,cfg.pptarget2djeol,'IgnoreCase',true) %UFZS dataset
            dirs2d(numel(dirs2d)+1)={fullfile(jdirnam,jdirfile)};
        else %conventional ZS dataset
            dirsConv(numel(dirsConv)+1)={fullfile(jdirnam,jdirfile)};
        end
    end
end

% Check to see if only 1D data, or only 2D data, or both, were found. If
% both, have user specify which to process
if ~isempty(dirs1d) && ~isempty(dirs2d) %both 1D and 2D datasets found: 
    %have user choose which to process
    prompt='Choose whether to process 1D or 2D data:';
    choices={'1D','2D'};
    answer = listdlg('ListString',choices,'SelectionMode','single',...
        'ListSize',[200 50],'PromptString',prompt);
    if answer == 1
        disp('1D processing selected. User will select all to process....')
        pflgs.proc1dflg=true;
        dirs=dirs1d;
    else
        disp('2D processing selected. User will process only one....')
        pflgs.proc1dflg=false; 
        dirs=dirs2d;        
    end
    pflgs.procConvflg=false;
elseif ~isempty(dirs1d) && isempty(dirs2d) %only 1D datasets found
    disp('Only 1D datasets found. User will select all to process....')
    pflgs.proc1dflg=true;
    pflgs.procConvflg=false;
    dirs=dirs1d;    
elseif isempty(dirs1d) && ~isempty(dirs2d) %only 2D datasets found
    disp('Only 2D datasets found. User will process only one....')
    pflgs.proc1dflg=false;
    pflgs.procConvflg=false;
    dirs=dirs2d;    
elseif ~isempty(dirsConv)
    disp(['No UFZS datasets found - only conventional z-spectroscopy data. ' ...
        'All 2D slices will be processed....'])
    pflgs.proc1dflg=false;
    pflgs.procConvflg=true;
    dirs=dirsConv;
else %no datasets found
    error(['No datasets found! Check that pptarget fields in struct config '...
        'in InitUserSettingsUFZS.m matches (part of) the pulse sequence '...
        'file name'])
end
dirs=num2cell(sort(str2double(dirs))); %sort in ascending order
dirs=cellfun(@num2str,dirs,'UniformOutput',false); %convert back to cell
    %array of strings

% This adds back in any JEOL directores, which would have been converted to 'NaN'
if pflgs.procConvflg
    dirs(strcmp(dirs,'NaN'))=dirsConv;     
else
    dirs(strcmp(dirs,'NaN'))=dirs2d(contains(dirs2d,'.jdx')); 
end

% Have user choose which datasets to process (if not just one 2D one)
%
if (pflgs.proc1dflg || numel(dirs) ~= 1) && ~pflgs.procConvflg
    choices=dirs;
    if pflgs.proc1dflg %set selection option to multiple for 1D processing
        prompt='Choose 1D experiment numbers to load:';
        selmod='multiple';
    else %set selection option to single for 2D processing
        prompt='Choose one 2D experiment number to load:';    
        selmod='single';
    end
    answer = listdlg('ListString' , choices , ...
        'SelectionMode' , selmod , ...
        'ListSize' , [200 500] , ...
        'PromptString' , prompt);
    dirs=dirs(answer);
elseif ~pflgs.procConvflg
    disp(['Only one 2D dataset detected: ' dirs{1} ' - processing...'])
end

% Detect whether the selection is a JEOL dataset
if contains(dirs,'.jdx')
    pflgs.jeol=true;
else
    pflgs.jeol=false;
end
end
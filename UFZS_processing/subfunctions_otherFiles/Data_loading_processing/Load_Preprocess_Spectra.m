% Load_Preprocess_Spectra: Loads in both the FID/spectral data and
% associated data parameters. If FID data is read in, the data are 
% pre-processed and Fourier transformed according to user specifications.
%
%   INPUTS:
%       pathn       -   String containing path to study directory
%       ddirs       -   Cell array of strings specifying the data
%                       directory/ies in the main study directory to be
%                       processed
%       pflgs       -   Struct containing values of logical flags for
%                       data processing options
%       ppars       -   Struct containing numerical/string processing 
%                       parameters 
%
%   OUTPUTS:
%       spec        -   2D array containing all processed spectra, with the
%                       frequency axis along the 2nd dimension
%       pars        -   Struct containing acquisition/processing
%                       parameter values extracted from associated
%                       parameter files
%
function [spec,pars]=Load_Preprocess_Spectra(pathn,ddirs,pflgs,ppars)
disp('Loading in data...')

if pflgs.topproc %pull in processed data from pdata folders
%     np = str2double(readParsTopspin(fullfile(pathname,datadirs{end},'pdata','1'),...
%         'procs',{'##$SI'})); %number of spectral points (should be equal for 
%         %all 1D datasets)
    if pflgs.proc1dflg %load in 1D datasets with loop
        % Read in 1st file, then fill in rest of 1st dimension with
        % zeros, then read in the rest
        [spec,pars] = ReadTopspinDataPars(fullfile(pathn,ddirs{1}),'1r');
        spec(2:numel(ddirs),:)=zeros(numel(ddirs)-1,size(spec,2));
        for i=2:numel(ddirs)
            [spec(i,:),pars(i)]=ReadTopspinDataPars(fullfile(pathn,ddirs{i}),'1r');
        end
    else %load in 2D dataset
        [spec,pars] = ReadTopspinDataPars(fullfile(pathn,ddirs{1}),'2rr');
    end
    spec = flip(spec,2); %flip so as to get ppm scale correct
else %pull in data from fid (1D) or ser (2D) file(s)
    if pflgs.jeol
        [fid,pars]=ReadJEOLjdx(fullfile(pathn,ddirs{end}));
        % Center FIDs via circular shifting (using 1st FID)
        maxval=max(abs(fid(1,:)));
        maxind=find(abs(fid(1,:))==maxval);
        np=size(fid,2);
        shift=ceil(np/2)-maxind;
        maxleft=abs(fid(1,maxind-1));
        maxright=abs(fid(1,maxind+1));
        if maxleft>maxright
            fid=circshift(fid,[0 shift+1]);
        else
            fid=circshift(fid,[0 shift]);
        end
        fid=flip(fid,2);
    else
%         np = str2double(readParsTopspin(fullfile(pathname,datadirs{end}),...
%             'acqus',{'##$TD'}))/2; %number of complex FID points (should be 
%             %equal for all 1D datasets)
        if pflgs.proc1dflg %load in 1D datasets with loop
            % Read in 1st file, then fill in rest of 1st dimension with
            % zeros, then read in the rest
            [fid,pars] = ReadTopspinDataPars(fullfile(pathn,ddirs{1}),'fid');
            fid(2:numel(ddirs),:)=zeros(numel(ddirs)-1,size(fid,2));
            for i=2:numel(ddirs)
                %NOTE: THIS MAY NOT WORK WITH TOPSPIN 4 YET
                [fid(i,:),pars(i)] = ReadTopspinDataPars(fullfile(pathn,...
                    ddirs{i}),'fid');
            end    
        else %load in 2D dataset
            [fid,pars] = ReadTopspinDataPars(fullfile(pathn,ddirs{1}),...
                'ser');
        end
    end
    
    % Apodize, zerofill, Fourier transform in spectral dimension only
    spec=ProcessFIDdata(fid,ppars,pflgs,pars);
end
end
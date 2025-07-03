% ProcessFIDdata: Takes 2D NMR data and performs apodization/filtering,
% zero-filling, and Fourier transform based upon specified parameters.
%
%   INPUTS:
%       fid         -   2D array containing FIDs to be processed, with the 
%                       time axis along the 2nd dimension
%       ppars       -   Struct containing numerical/string processing 
%                       parameters 
%       pflgs       -   Struct containing values of logical flags for data
%                       processing options and specifications
%       datapars    -   Struct containing values extracted from data
%                       parameter files
%
%   OUTPUTS:
%       spec        -   2D array containing processed spectra, with the
%                       frequency axis along the 2nd dimension
%
function spec=ProcessFIDdata(fid,ppars,pflgs,datapars)
% ID FID data sizes, depending on whether ultrafast or conventional
% Z-spectroscopy. At the same time, truncate the points prior to the FID
% max for all points
if pflgs.procConvflg %3D dataset!
    % First, lop off the first FID points
    truncPts=35; %# of points to truncate
    fid(:,:,1:truncPts)=[];

    fSize(1)=size(fid,2);
    fSize(2)=size(fid,3);
    fSize(3)=size(fid,1);
else
    fSize(1)=size(fid,1);
    fSize(2)=size(fid,2);
end

% Calculate the filter to apply to the FID
if strcmp(ppars.filter,'exponential')
    % First, we'll need the dwell time (in s) to generate a time array
    if pflgs.jeol
        np=datapars(1).x_X_CURR_POINTS;
        aqStart=datapars(1).x_X_START;
        aqEnd=datapars(1).x_X_STOP;
        dw=(aqEnd-aqStart)./(np-1);
    else
        dw = 1/datapars(1).sw/2; %dwell time 
    end
    
    % If ultrafast z-spec: max is centered at middle point
    % If conventional z-spec: max is at very first point
    if pflgs.procConvflg
        t=(0:(fSize(2)-1))*dw;
        apdvec = exp(-t*ppars.ap*pi);
    else
        t=(0:ceil(fSize(2)/2-1))*dw;
        apdvec = exp((t-max(t))*ppars.ap*pi);
        apdvec(floor(fSize(2)/2+1):fSize(2)) = exp(-t*ppars.ap*pi);
    end
elseif strcmp(ppars.filter,'gaussian')
    %% NEED TO FIX FOR CONVENTIONAL Z-SPEC! (I THINK CHANGE gmean=1 IS ENOUGH?)
    % Calculate parameters for Gaussian filter, generate
    gmean = fSize(2)/2 + 0.5;
    sigma = sqrt(5/log(10)*(gmean - 1)^2 / ppars.edge); %calculate sigma values for Gaussian
    apdvec = exp(-1 * ((1:fSize(2)) - gmean).^2 ./ 2 ./ (sigma^2));
end

if pflgs.procConvflg %need to work with a 3D dataset
    apdvec2(1,1,:)=apdvec;
    apdmat=repmat(apdvec2,[fSize(3) fSize(1) 1]);
    padSizes=[0,0,floor(((fSize(2)+truncPts) * (ppars.zf - 1)))+truncPts]; 
        %add zeros for the truncated points, too!
    padDir='post';
    specDim=3;
else
    apdmat=repmat(apdvec,[fSize(1) 1]);
    padSizes=[0,floor((fSize(2) * (ppars.zf - 1))/2)];
    padDir='both';
    specDim=2;
end

% Apodize, zerofill, Fourier transform in spectral dimension only
fidap=fid.*apdmat; %apodize
%     fidap(:,[1:ppars.zeropts (size(fid,2)-ppars.zeropts+1):size(fid,2)]) = 0; 
%         %zero out selected points
fidzf=padarray(fidap,padSizes,padDir); 
    %zerofill in spectral dimension
spec=fftshift(fft(fftshift(fidzf,specDim),[],specDim),specDim); 
    %cut out 1st points for group delay compensation
end
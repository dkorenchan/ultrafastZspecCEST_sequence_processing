% calcZspecMTRasym: Calculated MTR, as well as MTR asymmetry. 
%
%   INPUTS:
%       results     -   Struct containing raw spectral data and other
%                       parameters describing the data
%       pflgs       -   Struct containing logicals indicating processing 
%                       information. Subfields relevant to this function 
%                       are:
%                           .procConvflg    -   whether data are from
%                                               conventional Z-spectroscopy
%                                               instead of ultrafast
%                           .norm           -   whether spectral normalization  
%                                               was performed prior to this 
%                                               function
%       params      -   Struct containing numerical/string processing 
%                       parameters 
%       nosatidx    -   Vector containing the indices of the spectra
%                       corresponding with zero saturation amplitude. These
%                       correspond with the first dimension of
%                       results.spec.
%       normpars    -   Struct containing parameters relevant to the
%                       'nosatfit' normalization, if performed
%
%   OUTPUTS:
%       results     -   Struct containing processed results, now updated  
%                       with calculated z-spectra and MTR asymmetry
%
function results=calcZspecMTRasym(results,pflgs,params,nosatidx,normpars)
% Determine whether data are from conventional or ultrafast Z-spectroscopy
if pflgs.procConvflg
    % Reorder FIDs and results.specppm to be monotonically decreasing
    [results.satppm,sortIdx]=sort(results.satppm,'descend');
    results.spec=results.spec(:,sortIdx,:);

    % Integrate absolute value spectra
    integRegionPpm=[-1 1]; %values to integrate between, in ppm
    integIdx=(results.specppm>=integRegionPpm(1)) & ...
        (results.specppm<=integRegionPpm(2));
    results.specintegrals=squeeze(sum(abs(results.spec(:,:,integIdx)),3));

    % Calculate Z using largest abs value ppm value to normalize
    refIdx=find(results.satppm==max(abs(results.satppm)));
    results.zspec=results.specintegrals./results.specintegrals(:,refIdx);

    % Remove largest abs val ppm value from data
    results.zspec(:,refIdx)=[]; 
    results.zspecppm=results.satppm;
    results.zspecppm(refIdx)=[]; 

    results.zspeclabels = results.speclabels; 

    npz=size(results.zspec,2);

    results.zasymppm = flip(results.zspecppm(1:ceil(end/2)));
else
    % ID which non-saturation spectrum to use for
    % calculating the z-spectra
    if length(nosatidx) > 1 %use the spectrum with the higher amplitude
        nosatmax = max(abs(results.spec(nosatidx,:)')); %ID max amplitude for each
        refind = nosatidx(nosatmax == max(nosatmax));
        refind = refind(1); %in case somehow they are equal!
    %     refind = nosatidx(end); %use final acquisition for MTR calc
    else
        refind = nosatidx;
    end
    
    % Constrain spectral window to values determined by params.ppmwdw
    wdw=find(abs(results.specppm) - params.ppmwdw < 0); %indices of z-spectrum to keep
    
    % If using 'nosatfit' normalization: calculate z-spectra by clustering 
    % subspectra into groups, each with its
    % own associated non-saturated spectrum. Then use the appropriate
    % non-saturated reference spectrum with each to calculate the z-spectrum.
    % Otherwise, use the largest-amplitude non-saturated spectrum.
    if pflgs.norm && strcmp(params.normtype,'nosatfit') && length(nosatidx)>=2
    %     segVec=kmeans(abs(results.spec(:,results.specppm<-params.ppmwdw|...
    %         results.specppm>params.ppmwdw)),length(nosatidx)); %k-means clustering 
    %         %of spectra into groups equal to the # of non-saturated spectra
        results.zspec=zeros(size(results.spec,1),length(wdw));
        if length(unique(normpars.segVec(nosatidx)))<2 %if non-saturated spectra clustered 
            %together, take the mean and apply to all z-spectra
            results.zspec = abs(results.spec(:,wdw)) ./ repmat(...
                abs(mean(results.spec(nosatidx,wdw),1)),size(results.spec,1),1);    
        else
            for i=1:length(normpars.segVec)
                results.zspec(i,:) = abs(results.spec(i,wdw)) ./ ...
                    abs(results.spec(nosatidx(...
                    normpars.segVec(nosatidx)==normpars.segVec(i)),wdw));
                    %divide by the appropriate non-saturated profile, using segVec
            end
        end
    else
        results.zspec = abs(results.spec(:,wdw)) ./ repmat(...
            abs(squeeze(results.spec(refind,wdw))),size(results.spec,1),1);
    end
    
    results.zspec(nosatidx,:) = []; %remove z-spectra for datasets 
        %without saturation
    npz = length(wdw); %update np again 

    results.zspecppm = results.specppm(wdw);

    results.zspeclabels = results.speclabels;
    results.zspeclabels(nosatidx) = []; %eliminate legend entries for 
        %no-saturation datasets

    awdw = wdw(1:ceil(npz/2));
    results.zasymppm = flip(results.specppm(awdw));    
end

% Calculate MTR asymmetry
zright = results.zspec(:,ceil(npz/2):-1:1); %remember, data are flipped when plotted!
zleft = results.zspec(:,floor(npz/2+1):end);
results.zasym = zleft - zright;
end
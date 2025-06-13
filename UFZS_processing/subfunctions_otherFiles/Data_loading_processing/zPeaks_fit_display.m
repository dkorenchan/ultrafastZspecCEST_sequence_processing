% zPeaks_fit_display: Coordinating function to fit peaks in the MTR 
% z-spectral plots to the desired peak-fitting function, then plot the 
% results
%
%   INPUTS:
%       results     -   Struct containing processed z-spectral data as well
%                       as other information about the dataset
%       pars        -   Struct containing other processing parameters,
%                       including pool names for z-peak fitting
%       pflgs       -   Struct containing values of logical flags for
%                       data processing options
%
%   OUTPUTS:
%       results     -   Struct comprised of the input results but now also
%                       containing results from z-peak fitting
%
function results=zPeaks_fit_display(results,ppars,pflgs)
set(0,'DefaultFigureWindowStyle','docked') %this will keep all fitting 
    %figs together
pvals=1:size(results.zspec,1); %list of (other) power indices to fit

% Set up fitting parameters, based upon which type of peaks
if strcmp(ppars.peaktype,'Pseudo-Voigt')
    disp('Performing Pseudo-Voigt peak fitting of z-spectra...')
    npar=6; % # of parameters defining each Pseudo-Voigt peak
    pfitvals=setPVPeakBounds;
elseif strcmp(ppars.peaktype,'Lorentzian')
    disp('Performing Lorentzian peak fitting of z-spectra...')
    npar=4; % # of parameters defining each Lorentzian peak
    pfitvals=setLPeakBounds;
end

% If MT is to be fit to a super-Lorentzian, override MT fitting parameters
% accordingly
if pflgs.MTsuperLorentz
    pfitvals.MT=setSLPeakBounds;
    fixvals.MTsuperLorentz=true;
else
    fixvals.MTsuperLorentz=false;
end

% Set up variable fixvals to use for fixing/constraining values (NaN 
% for all values free to vary). If all Pseudo-Voigt peaks are to have 
% identical Lorentzian-Gaussian character (based upon 
% pflgs.PVcharconst), include a field fixvals.samePVchar=true (to be
% read in by ufzsMultiPeakFit.m)
for i = 1:numel(ppars.pools)
    results.peakfit.(ppars.pools{i})=zeros(size(results.zspec,1),1);
    results.peakLW_Hz.(ppars.pools{i})=zeros(size(results.zspec,1),1);
    fixvals.(ppars.pools{i})=NaN(npar,1);
end    
if strcmp(ppars.peaktype,'Pseudo-Voigt') && pflgs.PVcharconstr 
    %set flag to constrain peak character
    fixvals.samePVchar=true;
end

% If using offsets from specified power, move that power to the front
% of index list pvals
if pflgs.fix
    pvals = [ppars.fixind pvals(pvals ~= ppars.fixind)];
end

% Check for water, NOE, and MT pools specified, in order to fit first (if
% pflgs.water1st=true)
firstnames = ppars.pools(strcmp(ppars.pools,'water')|...
    strcmp(ppars.pools,'NOE')|strcmp(ppars.pools,'MT'));            

% Start loop for fitting z-spectra for each saturation power
ind_fits=cell(length(pvals),1);
for i = pvals

%         % If pflgs.phaseMTR=true, first add linear phase to the 
%         % spectrum such that the signal amplitudes at each end of the 
%         % spectrum are equal
%         if pflgs.phaseMTR
%             linphase=lsqnonlin(@(x) MTRphase(x,results.zspecppm,...
%                 1-squeeze(results.zspec(i,:))));
%             MTR=(1-squeeze(results.zspec(i,:))).*exp(-1i)
%         else
%         end

    % If pflgs.water1st=true, fit for water, NOE, and MT pools using 
    % negative ppm values only prior to fitting other peaks
    if pflgs.water1st
        EPfirst = ufzsMultiPeakFit(results.zspecppm, ...
            1-squeeze(results.zspec(i,:)),results.omega_0_MHz,...
            firstnames,pfitvals,fixvals);
        for j = 1:numel(firstnames)
            fixvals.(firstnames{j}) = EPfirst.(firstnames{j});
        end
    end
    
    % Then, fit the rest of the peaks
    [results.EstParams{i},~,~,~,ind_fits{i}] = ufzsMultiPeakFit(...
        results.zspecppm,1-squeeze(results.zspec(i,:)),results.omega_0_MHz,...
        ppars.pools,pfitvals,fixvals,ppars.ppmwt,true);
    title([ppars.peaktype ' fitting of z-spectrum, ' ...
        num2str(results.satT(i),'%1.2f') ' \muT (' ...
        num2str(results.satHz(i),'%3.1f') ' Hz)'])

    % Store peak amplitudes + LWs in results, and reset fixvals. If  
    % peak offsets will be fit only once, detect which power index they  
    % are fit for, then store peak offsets in fixvals. 
    for j = 1:numel(ppars.pools)
        name=ppars.pools{j};
        results.peakfit.(name)(i)=results.EstParams{i}.(name)(1);
%             if strcmp(pars.peaktype,'Pseudo-Voigt') %calculate FWHM from ind_fits
%                 PSshape_normd=ind_fits{i}.(name)./max(ind_fits{i}.(name));
%                 FWHM_ppmvals=results.zspecppm(abs(PSshape_normd-0.5)<1e-2);
%                 delta_ppm=abs(FWHM_ppmvals-circshift(FWHM_ppmvals,1));
%                 delta_ppm=abs(delta_ppm(2:end)); %remove 1st value, since it circles round
%                 maxind=find(delta_ppm==max(delta_ppm));
%                 ind_lr=floor(length(delta_ppm)/4);
%                 peakLW_ppm=sum(delta_ppm((maxind-ind_lr):(maxind+ind_lr)));
%             else
%                 peakLW_ppm=results.EstParams{i}.(name)(2);
%             end
%             results.peakLW_Hz.(name)(i)=peakLW_ppm*omega_0_MHz;
%                 %convert from ppm to Hz
        fixvals.(name)=NaN(npar,1);
        if pflgs.fix
            if strcmp(name,'MT') && pflgs.MTsuperLorentz
                fixvals.(name)(3)=results.EstParams{ppars.fixind}.(name)(3);
            else
                fixvals.(name)(end-1)=results.EstParams{ppars.fixind}.(name)(end-1); 
            end
        end
    end
end
results.indivFits=ind_fits;
set(0,'DefaultFigureWindowStyle','normal') %turn off figure docking
end
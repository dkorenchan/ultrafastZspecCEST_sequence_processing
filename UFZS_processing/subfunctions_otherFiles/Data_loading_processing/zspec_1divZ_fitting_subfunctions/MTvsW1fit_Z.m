function [FitResult,FitParam] = MTvsW1fit_Z(Offset,Saturation,satpwr_uT,FitParam)
% Fit the far-ppm values of multiple Z-spectra at different saturation 
% amplitudes using an MT super-Lorentzian model
%   Offset      -   Z-spectrum offsets, in ppm
%   Saturation  -   Z-spectrum amplitudes
%   FitParam    -   Structure containing fitting parameters

if size(Offset,1) == 1
    Offset = Offset';
end

if size(Saturation,1) ~= size(Offset,1)
    Saturation = Saturation';
end

if (Saturation(1,1)) > 10
    error('Z-spectrum have to be 0-1, can not be percentage')
end

if FitParam.MT.ppmExclude(1)>FitParam.MT.ppmExclude(2) %make sure 1st value is the smaller one!
    FitParam.MT.ppmExclude=FitParam.MT.ppmExclude([2,1]);
end

if FitParam.MT.ppmReinclude(1)>FitParam.MT.ppmReinclude(2) %make sure 1st value is the smaller one!
    FitParam.MT.ppmReinclude=FitParam.MT.ppmReinclude([2,1]);
end

if size(satpwr_uT,1)>size(satpwr_uT,2)
    satpwr_uT=satpwr_uT';
end

satpwr_radS=satpwr_uT*42.577*2*pi;  %convert saturation amplitudes to rad/s
OffsetHz=Offset*FitParam.Magfield;  %convert offsets to Hz


%% MT FITTING
% First, take only the part of the z-spectrum that is far from the regular
% CEST peaks, as indicated in FitParam.MTppmExclude
MTfitIdx=Offset<FitParam.MT.ppmExclude(1)|Offset>FitParam.MT.ppmExclude(2)...
    |(Offset>FitParam.MT.ppmReinclude(1) & Offset<FitParam.MT.ppmReinclude(2));
Saturation_MTfit=Saturation(MTfitIdx,:)';
Offset_MTfit=OffsetHz(MTfitIdx);

% If FitParam.MT.wtHigherSatAmpl==true, weight fitting towards the higher
% saturation amplitudes
if FitParam.MT.wtHigherSatAmpl
    disp('MT fitting: Focus on larger saturation amplitudes...')
    reswt=reshape(ones(size(Saturation_MTfit)).*(satpwr_uT)'./mean(satpwr_uT),...
        [],1);
else
    reswt=ones(numel(Saturation_MTfit),1);
end

% Set up fitting parameters
opts=fitoptions('Method','NonlinearLeastSquares','Display','iter',...
    'DiffMaxChange',1e-3,'DiffMinChange',1e-12,'MaxFunEvals',6000,'MaxIter',4000,...
    'TolFun',1e-15,'TolX',1e-15,'Weights',reswt);
opts.StartPoint =[0,    1,      .05,    40,     10e-6,    -0.5*FitParam.Magfield];
opts.Lower      =[0,    1,      .005,   1,      1e-6,     -1*FitParam.Magfield];
opts.Upper      =[0,    10,     .3,     200,    100e-6,   1*FitParam.Magfield];

% Fit MT
if isfield(FitParam,'MT') && isfield(FitParam.MT,'lineshape')
    lineshape = FitParam.MT.lineshape;
else
    lineshape = 'superlorentzian';
end
% lineshape is captured by the closure; fittype does not see it as a parameter
modelstr= @(t2a,rb,mb0,r,t2b,delta_MT,ra,delta,w1) ...
    mt_model_forCESTfit(t2a,ra,rb,mb0,r,t2b,1,delta_MT,delta,w1,lineshape);
ft=fittype(modelstr,'problem','ra','independent',{'delta','w1'},'dependent','z');
[delta_fit,w1_fit,Z_fit]=prepareSurfaceData(Offset_MTfit,satpwr_radS,...
    Saturation_MTfit); 
[FitResult.model,FitResult.resnorm,FitResult.fitoutput]=...
    fit([delta_fit,w1_fit],Z_fit,ft,opts,'problem',FitParam.R1);

% Generate MT-only curve by setting T2a=0 (which makes R_rfa=0)
cv=coeffvalues(FitResult.model);
FitResult.Z=mt_model_forCESTfit(0,FitParam.R1,cv(2),cv(3),cv(4),cv(5),1,cv(6),...
    OffsetHz,satpwr_radS,lineshape);
FitResult.coeffs=cv;


%% PLOTTING
if FitParam.MT.plot
    figure; scatter(Offset_MTfit,Saturation_MTfit); set(gca,'ColorOrderIndex',1);
    hold on; plot(OffsetHz,FitResult.Z); 
    title('MT Fitted Curve and Fit Points'); ylabel('Z'); xlabel('Offset (Hz)');
    leglbl=cell(numel(satpwr_uT),1);
    for ii=1:numel(satpwr_uT)
        leglbl{ii}=[num2str(satpwr_uT(ii),'%2.1f') ' \muT (' ...
            num2str(satpwr_radS(ii),'%3.0f') ' rad/s)'];
    end
    ylim([0 1])
    legend(leglbl);
    set(gca,'XDir','reverse');
end
end
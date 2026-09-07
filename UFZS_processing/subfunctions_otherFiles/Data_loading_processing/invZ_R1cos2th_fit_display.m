% invZ_R1cos2th_fit_display: Fits z-spectral peaks in the
% R1*cos^2(theta)*(1/Z-1) domain.  When the MT pool is selected, its
% contribution is first fitted in Z-space (MTvsW1fit_Z) and then used as a
% fixed background during the CEST peak fitting.  Fitted peak amplitudes
% are converted to equivalent Z-values (ZssPeak) and stored as
% (1 - ZssPeak) in results.peakfit, making them directly compatible with
% the existing QUESP_fit_display pipeline.
%
%   INPUTS:
%       results     -   Struct with processed z-spectral data, including:
%                         .zspec        (n_satpwr x n_ppm) Z-spectra
%                         .zspecppm     (1 x n_ppm) offset ppm axis
%                         .satT         saturation amplitudes in uT
%                         .omega_0_MHz  1H Larmor frequency in MHz
%       ppars       -   Processing parameters struct (pools, peaktype,
%                         PVcharconstr, fix, fixind, water1st, ppmwt, ...)
%       pflgs       -   Processing flags struct
%       timing      -   Timing struct (.tp saturation duration [s],
%                         .rd recovery delay [s])
%
%   OUTPUTS:
%       results     -   Updated struct with:
%                         .peakfit.(pool)  (1 - ZssPeak) per saturation power
%                         .T1w             water T1 [s] from user prompt
%                         .EstParams       fitted parameter cell array
%                         .indivFits       individual peak fit cell array
%       ppars       -   Updated params: MT removed from .pools when
%                         pre-fitted as background
%
function [results, ppars] = invZ_R1cos2th_fit_display(results, ppars, pflgs, timing)
set(0, 'DefaultFigureWindowStyle', 'docked')

n_satpwr    = size(results.zspec, 1);
omega_0_MHz = results.omega_0_MHz;
satpwr_uT   = results.satT;
satpwr_Hz   = results.satHz;
ppm         = results.zspecppm;
all_Zspec   = results.zspec;


%% PROMPT FOR T1
answer      = inputdlg({'Water 1H T1 (s):'}, '1/Z fitting: R1 parameter', 1, {'4.2'});
results.T1w = str2double(answer{1});
R1          = 1 / results.T1w;


%% SET UP PEAK FITTING PARAMETERS
if strcmp(ppars.peaktype, 'Pseudo-Voigt')
    disp('Performing Pseudo-Voigt peak fitting of R_1*cos^2(theta)*(1/Z-1) spectra...')
    npar     = 6;
    pfitvals = setPVPeakBounds_invZ;
elseif strcmp(ppars.peaktype, 'Lorentzian')
    disp('Performing Lorentzian peak fitting of R_1*cos^2(theta)*(1/Z-1) spectra...')
    npar     = 4;
    pfitvals = setLPeakBounds_invZ;
end

% MT is pre-fitted separately; only remaining pools enter the CEST fit
hasMT     = any(strcmp(ppars.pools, 'MT'));
cestPools = ppars.pools(~strcmp(ppars.pools, 'MT'));

% Initialise fixvals and peakfit storage for CEST pools
fixvals = struct;
for ii = 1:numel(cestPools)
    results.peakfit.(cestPools{ii})   = zeros(n_satpwr, 1);
    results.peakLW_Hz.(cestPools{ii}) = zeros(n_satpwr, 1);
    fixvals.(cestPools{ii})           = NaN(npar, 1);
end
if strcmp(ppars.peaktype, 'Pseudo-Voigt') && pflgs.PVcharconstr
    fixvals.samePVchar = true;
end
fixvals.MTsuperLorentz = false; % MT is handled via pre-fitting in this mode


%% STEP 1 - MT PRE-FITTING IN Z-SPACE (skipped when MT pool not selected)
if hasMT
    disp('***STEP 1: Fitting MT pool using high-ppm wings of the Z-spectrum...')
    dppmWings                   = 2;          % ppm range for each wing to fit for MT
    disp(['***Wing length on each side: ' num2str(dppmWings) ' ppm'])
    FitParam.MT.ppmExclude      = [round(min(ppm)+dppmWings), ...
                                   round(max(ppm)-dppmWings)];  % ppm range to exclude from MT fit
    if FitParam.MT.ppmExclude(1)>-4
        warning(['You may have rNOE contamination in MT fitting! ' ...
            'Consider a wider z-spectrum ppm display window'])
    end
    if FitParam.MT.ppmExclude(2)<6
        warning(['You may have CEST peak contamination in MT fitting! ' ...
            'Consider a wider z-spectrum ppm display window'])
    end
    FitParam.MT.ppmReinclude    = [0, 0];     % sub-range within exclusion to re-include
    FitParam.MT.wtHigherSatAmpl = false;
    FitParam.MT.fixUndercut     = false;
    FitParam.MT.plot            = true;
    FitParam.MT.lineshape       = ppars.MTlineshape;
    FitParam.R1                 = R1;
    FitParam.Magfield           = omega_0_MHz;
    [FitResult_MT, ~] = MTvsW1fit_Z(ppm, all_Zspec, satpwr_uT, FitParam);
    MT_Zfit = FitResult_MT.Z;
    % Ensure orientation matches all_Zspec (n_satpwr x n_ppm)
    if ~isequal(size(MT_Zfit), size(all_Zspec))
        MT_Zfit = MT_Zfit';
    end
    results.EstMT = FitResult_MT; % full MT fit result (model, coeffs, Z)
    disp('***STEP 1: COMPLETE!')
    % Remove MT from pool list since it is now encoded in the background
    ppars.pools = ppars.pools(~strcmp(ppars.pools, 'MT'));
else
    MT_Zfit = ones(size(all_Zspec));
    disp('No MT pool specified; using flat MT background (Z=1) for 1/Z fitting.')
end


%% STEP 2 - COMPUTE R1*(1/Z-1) AND FIT CEST PEAKS IN COS2TH-WEIGHTED SPACE
disp('***STEP 2: Computing R_1*(1/Z-1), then fitting CEST peaks...')

all_invZspec = R1 * (1 ./ all_Zspec - 1);
MT_invZ      = R1 * (1 ./ MT_Zfit  - 1);

% ppm range for CEST fitting (covers typical CEST/NOE peaks on both sides)
% ppmInclude = [-2, 12];
ppmInclude = [-2, 5];
CESTfitIdx = ppm > ppmInclude(1) & ppm < ppmInclude(2);
ppm_CEST   = ppm(CESTfitIdx);

% When fix mode is active, run the fix-index power first so its offsets
% are available for all subsequent powers
pvals = 1:n_satpwr;
if pflgs.fix
    pvals = [ppars.fixind, pvals(pvals ~= ppars.fixind)];
end

ind_fits  = cell(n_satpwr, 1);
EstParams = cell(n_satpwr, 1);

for ii = pvals
    w1_Hz      = satpwr_Hz(ii);
    cosSqTheta = (ppm_CEST .* omega_0_MHz).^2 ./ ...
                 ((ppm_CEST .* omega_0_MHz).^2 + w1_Hz.^2);

    fitspec    = all_invZspec(ii, CESTfitIdx) .* cosSqTheta;
    invZ_MT_bg = MT_invZ(ii, CESTfitIdx)      .* cosSqTheta;

    [EstParams{ii}, ~, ~, ~, ind_fits{ii}] = ...
        ufzsMultiPeakFit_1divZminus1xR1cos2th( ...
            ppm_CEST, w1_Hz, fitspec, invZ_MT_bg, R1, omega_0_MHz, ...
            cestPools, pfitvals, fixvals, ppars.ppmwt, true);
    title([ppars.peaktype ' fit of R_1cos^2\theta(1/Z-1)-spectrum, ' ...
        num2str(satpwr_uT(ii), '%1.2f') ' \muT (' ...
        num2str(w1_Hz * 2 * pi, '%3.1f') ' rad/s)'])

    % For each non-water pool, convert the fitted amplitude to an equivalent
    % Z-value (ZssPeak) using the QUESP steady-state expression, then store
    % (1 - ZssPeak) in results.peakfit for use by QUESP_fit_display.
    %
    % The cos^2(theta) terms cancel algebraically:
    %   ZssPeak = R1*cos2th / (p(1)*cos2th + R1*cos2th) = R1 / (p(1) + R1)
    % but we follow the explicit form from MT_CEST_fit_QUESP_dk.m for clarity.
    for jj = 1:numel(cestPools)
        name = cestPools{jj};
        if strcmp(name, 'water')
            continue
        end
        p           = EstParams{ii}.(name);
        RpeakOffset = p(end-1);   % ppm offset: index 3 (Lorentz) or 5 (PV)
        delta_Hz    = omega_0_MHz * RpeakOffset;
        cos2th_peak = delta_Hz.^2 ./ (w1_Hz.^2 + delta_Hz.^2);
        Rpeak       = p(1) * cos2th_peak;
        ZssPeak     = R1 * cos2th_peak ./ (Rpeak + R1 * cos2th_peak);
        results.peakfit.(name)(ii) = 1 - ZssPeak;
    end

    % Reset fixvals; re-apply offset constraints if pflgs.fix is active
    for jj = 1:numel(cestPools)
        name   = cestPools{jj};
        tmpfix = NaN(npar, 1);
        if pflgs.fix
            tmpfix(end-1) = EstParams{ppars.fixind}.(name)(end-1);
        end
        fixvals.(name) = tmpfix;
    end
    if strcmp(ppars.peaktype, 'Pseudo-Voigt') && pflgs.PVcharconstr
        fixvals.samePVchar = true;
    end
    fixvals.MTsuperLorentz = false;
end

results.EstParams = EstParams;
% Add MT model coefficients to each per-power EstParams entry
if hasMT
    for ii = 1:n_satpwr
        results.EstParams{ii}.MT = FitResult_MT.coeffs(:);
    end
end
results.indivFits = ind_fits;

%% RECONSTRUCT AND PLOT FITTED Z-SPECTRA (undocked)
set(0, 'DefaultFigureWindowStyle', 'normal')
satpwr_radS   = satpwr_uT * 42.577 * 2 * pi;
all_Zspec_fit = zeros(size(all_Zspec));
for ii = 1:n_satpwr
    w1_Hz_ii     = satpwr_Hz(ii);
    cos2th_full  = (ppm .* omega_0_MHz).^2 ./ ...
                   ((ppm .* omega_0_MHz).^2 + w1_Hz_ii.^2);
    R1cos2th_fit = MT_invZ(ii,:) .* cos2th_full;
    for jj = 1:numel(cestPools)
        pn = cestPools{jj};
        p  = EstParams{ii}.(pn);
        if strcmp(ppars.peaktype, 'Pseudo-Voigt')
            pool_fit = cos2th_full .* ufzsSinglePseudoVoigtModel(p, ppm);
        else
            pool_fit = cos2th_full .* ufzsSingleLorentzianModel(p, ppm);
        end
        R1cos2th_fit = R1cos2th_fit + pool_fit;
    end
    all_Zspec_fit(ii,:) = R1 * cos2th_full ./ (R1cos2th_fit + R1 * cos2th_full);
end
results.zspec_fit = all_Zspec_fit;
figure; hold on;
clrs   = lines(n_satpwr);
leglbl = cell(n_satpwr, 1);
for ii = 1:n_satpwr
    scatter(ppm, all_Zspec(ii,:), 15, clrs(ii,:), 'filled', 'HandleVisibility', 'off');
    plot(ppm, all_Zspec_fit(ii,:), 'Color', clrs(ii,:));
    leglbl{ii} = [num2str(satpwr_uT(ii),'%2.1f') ' \muT (' ...
                  num2str(satpwr_radS(ii),'%3.0f') ' rad/s)'];
end
legend(leglbl);
title('Fitted Z-spectra (from 1/Z fit)');
ylabel('Z'); xlabel('Offset (ppm)');
ylim([0 1]);
set(gca, 'XDir', 'reverse');

disp('***STEP 2: COMPLETE!')
end

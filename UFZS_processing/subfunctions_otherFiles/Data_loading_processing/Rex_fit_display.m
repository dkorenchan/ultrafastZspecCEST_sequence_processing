% Rex_fit_display: 2D joint fit of R1*cos^2(theta)*(1/Z-1) spectra using
%   the dominant first term of Zaiss & Bachert 2013 Eq. (23) with R2b=0:
%
%     Rex(Dw, B1) = fb * kb * dOmb^2/(B1^2+Dw^2) * B1^2/(B1^2+kb^2+Dwb^2)
%
%   where dOmb = chemical shift of pool b relative to water [Hz],
%   Dw = saturation offset from water [Hz], Dwb = Dw - dOmb [Hz].
%   Water is fitted with a 3-parameter Lorentzian.
%   MT (if selected) is pre-fitted in Z-space before Rex fitting.
%   Fitted parameters (fb, kb) are stored in results.Rex and the QUESP
%   step is skipped.
%
%   INPUTS / OUTPUTS: same as invZ_R1cos2th_fit_display
%
function [results, ppars] = Rex_fit_display(results, ppars, pflgs, timing)
set(0, 'DefaultFigureWindowStyle', 'docked')

n_satpwr    = size(results.zspec, 1);
omega_0_MHz = results.omega_0_MHz;
satpwr_uT   = results.satT;
satpwr_Hz   = satpwr_uT * 42.577;
ppm         = results.zspecppm;
all_Zspec   = results.zspec;


%% GET REX FITTING BOUNDS
pfitvals = setRexBounds;

% R1 is a free parameter in the Rex fit; use the start value for MT pre-fitting
R1_init = pfitvals.water.st(2);

% Determine pool lists
hasMT     = any(strcmp(ppars.pools, 'MT'));
cestPools = ppars.pools(~strcmp(ppars.pools, 'MT'));
rexPools  = cestPools(~strcmp(cestPools, 'water')); % non-water CEST pools


%% STEP 1 - MT PRE-FITTING (same as invZ mode)
if hasMT
    disp('***STEP 1: Fitting MT pool using high-ppm wings of the Z-spectrum...')
    dppmWings                   = 2;
    FitParam.MT.ppmExclude      = [round(min(ppm)+dppmWings), ...
                                   round(max(ppm)-dppmWings)];
    FitParam.MT.ppmReinclude    = [0, 0];
    FitParam.MT.wtHigherSatAmpl = false;
    FitParam.MT.fixUndercut     = false;
    FitParam.MT.plot            = true;
    FitParam.MT.lineshape       = ppars.MTlineshape;
    FitParam.R1                 = R1_init;  % initial R1 for MT pre-fit only
    FitParam.Magfield           = omega_0_MHz;
    [FitResult_MT, ~] = MTvsW1fit_Z(ppm, all_Zspec, satpwr_uT, FitParam);
    MT_Zfit = FitResult_MT.Z;
    if ~isequal(size(MT_Zfit), size(all_Zspec))
        MT_Zfit = MT_Zfit';
    end
    results.EstMT = FitResult_MT;
    disp('***STEP 1: COMPLETE!')
    ppars.pools = ppars.pools(~strcmp(ppars.pools, 'MT'));
else
    MT_Zfit = ones(size(all_Zspec));
    disp('No MT pool specified; using flat MT background for Rex fitting.')
end


%% STEP 2 - 2D JOINT REX FIT
disp('***STEP 2: 2D joint Rex fit across all offsets and B1 values...')

% ppm fitting range (same as invZ mode)
ppmInclude = [-2, 5];
CESTfitIdx = ppm > ppmInclude(1) & ppm < ppmInclude(2);
ppm_CEST   = ppm(CESTfitIdx);
ppm_Hz     = ppm_CEST * omega_0_MHz;

% cos^2(theta) matrix: (n_satpwr x n_ppm)
cos2th = zeros(n_satpwr, numel(ppm_CEST));
for ii = 1:n_satpwr
    cos2th(ii,:) = ppm_Hz.^2 ./ (ppm_Hz.^2 + satpwr_Hz(ii).^2);
end

% Pass raw Z-spectra sliced to fitting range so R1 can be a fitted parameter
Zspec_CEST   = all_Zspec(:, CESTfitIdx);
MT_Zfit_CEST = MT_Zfit(:,  CESTfitIdx);

% Build parameter vector: [R2w(1), R1(1), pool1(4), pool2(4), ...]
%   R2w: water T2 rate [s^-1]; R1: fitted water T1 rate [s^-1]
%   each Rex pool: [fb, kb_Hz, delta_ppm, R2b_Hz]  (Gamma approx: kb>>R2b)
p0 = pfitvals.water.st(:)';
lb = pfitvals.water.lb(:)';
ub = pfitvals.water.ub(:)';
for jj = 1:numel(rexPools)
    name = rexPools{jj};
    p0 = [p0, pfitvals.(name).st(:)']; %#ok<AGROW>
    lb = [lb, pfitvals.(name).lb(:)']; %#ok<AGROW>
    ub = [ub, pfitvals.(name).ub(:)']; %#ok<AGROW>
end

opts_fit = optimoptions('lsqnonlin','MaxFunctionEvaluations',6000,'MaxIterations',4000,...
    'FunctionTolerance',1e-12,'StepTolerance',1e-12,'Display','iter');

[p_opt, ~, residuals, ~, ~, ~, jacobian] = lsqnonlin(...
    @(p) rexResidual(p, Zspec_CEST, MT_Zfit_CEST, ppm_Hz, satpwr_Hz, cos2th, rexPools, omega_0_MHz),...
    p0, lb, ub, opts_fit);

CI_all    = nlparci(p_opt, residuals, 'jacobian', jacobian);
R1_fitted = p_opt(2);
all_invZspec = R1_fitted * (1 ./ all_Zspec - 1);
MT_invZ      = R1_fitted * (1 ./ MT_Zfit  - 1);
data_2D  = R1_fitted * (1./Zspec_CEST - 1./MT_Zfit_CEST) .* cos2th;
data_vec = data_2D(:);
SS_res = sum(residuals.^2);
SS_tot = sum((data_vec - mean(data_vec)).^2);
Rsq    = 1 - SS_res / SS_tot;

% Store global R²
results.Rex.Rsq = Rsq;
disp(['Rex fit overall R^2 = ' num2str(Rsq,'%0.4f')])

%% EXTRACT RESULTS PER POOL AND REPORT
% Water
R2w = p_opt(1);
results.Rex.water.R2w = R2w;
results.Rex.water.R1  = R1_fitted;
results.Rex.water.CI  = CI_all(1:2,:);
results.T1w = 1 / R1_fitted;
disp(['  Water: R2w = ' num2str(R2w,'%0.3f') ' s^-1,  R1 = ' ...
    num2str(R1_fitted,'%0.3f') ' s^-1  (T1 = ' num2str(1/R1_fitted,'%0.2f') ' s)'])

% CEST pools
for jj = 1:numel(rexPools)
    name   = rexPools{jj};
    idx    = 2 + 4*(jj-1) + 1;    % pool j at p(4j-1); p(1)=R2w, p(2)=R1
    fb     = p_opt(idx);
    kb     = p_opt(idx+1);
    dp     = p_opt(idx+2);
    R2b    = p_opt(idx+3);
    CI_p   = CI_all(idx:idx+3, :);

    results.Rex.(name).fb        = fb;
    results.Rex.(name).kb        = kb;
    results.Rex.(name).delta_ppm = dp;
    results.Rex.(name).R2b       = R2b;
    results.Rex.(name).CI        = CI_p;

    disp(['  ' name ': fb = ' num2str(fb,'%0.4e') ' [' ...
        num2str(CI_p(1,1),'%0.4e') ', ' num2str(CI_p(1,2),'%0.4e') ']'...
        ',  kb = ' num2str(kb,'%0.1f') ' Hz [' ...
        num2str(CI_p(2,1),'%0.1f') ', ' num2str(CI_p(2,2),'%0.1f') ']'...
        ',  delta = ' num2str(dp,'%0.3f') ' ppm'...
        ',  R2b = ' num2str(R2b,'%0.1f') ' s^-1'])

    % Derive per-power ZssPeak for downstream compatibility
    % At the peak centre (Dw_b = 0): Rex = fb*kb*dOmb^2/(B1^2+dOmb^2)
    results.peakfit.(name) = zeros(n_satpwr, 1);
    dOmb_Hz = dp * omega_0_MHz;
    for ii = 1:n_satpwr
        B1       = satpwr_Hz(ii);
        omega1   = 2*pi*B1;
        dOmb_rad = 2*pi*dOmb_Hz;
        bpeak_at_peak = omega1^2 / (omega1^2 + kb^2);
        term1 = fb * kb * dOmb_rad^2 / (omega1^2 + dOmb_rad^2) * bpeak_at_peak;
        term2 = fb * R2b * bpeak_at_peak;
        ZssPeak     = R1_fitted / (term1 + term2 + R1_fitted);
        results.peakfit.(name)(ii) = 1 - ZssPeak;
    end
end
disp('***STEP 2: COMPLETE!')


%% PLOT 1: Per-B1 docked figures in R1*cos2th*(1/Z-1) space
% (DefaultFigureWindowStyle is already 'docked' from function start)
satpwr_radS = satpwr_uT * 42.577 * 2 * pi;
plotvis     = {'g-','r-','m-','c-','y-','b-'};
R2w = p_opt(1);

for ii = 1:n_satpwr
    B1        = satpwr_Hz(ii);
    cos2th_ii = cos2th(ii,:);
    mt_ii     = MT_invZ(ii,CESTfitIdx) .* cos2th_ii;
    omega1_w  = 2*pi*B1;
    Dw_rad_w  = 2*pi*ppm_Hz;
    water_ii  = cos2th_ii .* (R2w * omega1_w^2 ./ (omega1_w^2 + Dw_rad_w.^2));

    pool_contribs = zeros(numel(rexPools), numel(ppm_CEST));
    for jj = 1:numel(rexPools)
        idx  = 2 + 4*(jj-1) + 1;
        fb   = p_opt(idx);  kb = p_opt(idx+1);  dp = p_opt(idx+2);  R2b = p_opt(idx+3);
        dOmb_rad = 2*pi * dp * omega_0_MHz;
        Dw_rad   = 2*pi * ppm_Hz;
        Dwb_rad  = Dw_rad - dOmb_rad;
        omega1   = 2*pi * B1;
        bpeak    = omega1.^2./(omega1.^2+kb.^2+Dwb_rad.^2);
        Rex  = fb.*bpeak.*(kb.*dOmb_rad.^2./(omega1.^2+Dw_rad.^2) + R2b);
        pool_contribs(jj,:) = cos2th_ii .* Rex;
    end
    total_ii = water_ii + sum(pool_contribs, 1);

    figure; h = axes;
    plot(ppm_CEST, all_invZspec(ii,CESTfitIdx).*cos2th_ii, 'Marker','o', 'Color',[0.5 0.5 0.5], ...
        'MarkerFaceColor',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','none')
    hold on;
    plot(ppm_CEST, mt_ii,    'b-',  'LineWidth',2)
    plot(ppm_CEST, total_ii, 'k-',  'LineWidth',2)
    plot(ppm_CEST, water_ii, 'b--', 'LineWidth',1.5)
    legstrs = {'Raw data','MT','Sum','Water'};
    for jj = 1:numel(rexPools)
        plot(ppm_CEST, pool_contribs(jj,:), ...
            plotvis{1+mod(jj,numel(plotvis))}, 'LineWidth',2)
        legstrs{end+1} = rexPools{jj}; %#ok<AGROW>
    end
    axis([min(ppm_CEST) max(ppm_CEST) 0 Inf]);
    set(h,'Xdir','reverse')
    xlabel('\Delta\omega [ppm]')
    ylabel('R_1cos^2\theta(1/Z-1)')
    legend(legstrs)
    grid on
    title(['Rex fit of R_1cos^2\theta(1/Z-1)-spectrum, ' ...
        num2str(satpwr_uT(ii),'%1.2f') ' \muT (' ...
        num2str(satpwr_radS(ii),'%3.1f') ' rad/s)'])
end

%% PLOT 2: Reconstructed Z-spectra (undocked)
set(0, 'DefaultFigureWindowStyle', 'normal')
ppm_Hz_full   = ppm * omega_0_MHz;
all_Zspec_fit = zeros(size(all_Zspec));
for ii = 1:n_satpwr
    B1           = satpwr_Hz(ii);
    cos2th_full  = ppm_Hz_full.^2 ./ (ppm_Hz_full.^2 + B1^2);
    omega1_f     = 2*pi*B1;
    Dw_rad_full  = 2*pi*ppm_Hz_full;
    water_B1     = R2w * omega1_f^2 ./ (omega1_f^2 + Dw_rad_full.^2);
    R1c2th_total = MT_invZ(ii,:).*cos2th_full + cos2th_full.*water_B1;
    for jj = 1:numel(rexPools)
        idx  = 2 + 4*(jj-1) + 1;
        fb   = p_opt(idx);  kb = p_opt(idx+1);  dp = p_opt(idx+2);  R2b = p_opt(idx+3);
        dOmb = dp * omega_0_MHz;
        Dwb_rad_full = 2*pi*(ppm_Hz_full - dOmb);
        omega1_sq_f  = (2*pi*B1)^2;
        bpeak_full   = omega1_sq_f./(omega1_sq_f+kb^2+Dwb_rad_full.^2);
        rex_full  = fb.*bpeak_full.*(kb.*dOmb^2./(B1^2+ppm_Hz_full.^2) + R2b);
        R1c2th_total = R1c2th_total + cos2th_full.*rex_full;
    end
    all_Zspec_fit(ii,:) = R1_fitted * cos2th_full ./ (R1c2th_total + R1_fitted * cos2th_full);
end
results.zspec_fit = all_Zspec_fit;

clrs = lines(n_satpwr);
figure; hold on;
leglbl = cell(n_satpwr, 1);
for ii = 1:n_satpwr
    scatter(ppm, all_Zspec(ii,:), 15, clrs(ii,:), 'filled', 'HandleVisibility','off');
    plot(ppm, all_Zspec_fit(ii,:), 'Color', clrs(ii,:));
    leglbl{ii} = [num2str(satpwr_uT(ii),'%2.1f') ' \muT (' ...
                  num2str(satpwr_radS(ii),'%3.0f') ' rad/s)'];
end
legend(leglbl);
title('Reconstructed fitted Z-spectra (Rex fit)');
ylabel('Z'); xlabel('Offset (ppm)');
ylim([0 1]);
set(gca,'XDir','reverse');
end


%% ---- LOCAL HELPER FUNCTIONS ----

% rexResidual: R1 is p(2); data is recomputed from raw Z-spectra each iteration
function res = rexResidual(p, Zspec_CEST, MT_Zfit_CEST, ppm_Hz, satpwr_Hz, cos2th, rexPools, omega_0_MHz)
R1 = p(2);
data_2D   = R1 * (1./Zspec_CEST - 1./MT_Zfit_CEST) .* cos2th;
model_vec = rexModel(p, ppm_Hz, satpwr_Hz, cos2th, rexPools, omega_0_MHz);
res       = model_vec(:) - data_2D(:);
end


% % rexCESTResidual/rexCESTModel: CEST-only versions (no water, idx starts at 1)
% function res = rexCESTResidual(p, data_vec, ppm_Hz, satpwr_Hz, cos2th, rexPools, omega_0_MHz)
% res = rexCESTModel(p, ppm_Hz, satpwr_Hz, cos2th, rexPools, omega_0_MHz);
% res = res(:) - data_vec;
% end

% function model = rexCESTModel(p, ppm_Hz, satpwr_Hz, cos2th, rexPools, omega_0_MHz)
% n_satpwr = numel(satpwr_Hz);
% model    = zeros(n_satpwr, numel(ppm_Hz));
% for jj = 1:numel(rexPools)
%     idx     = 4*(jj-1) + 1;   % pool j starts at index 4*(j-1)+1 (no water)
%     fb      = p(idx);  kb = p(idx+1);  dp = p(idx+2);  R2b = p(idx+3);
%     dOmb_Hz = dp * omega_0_MHz;
%     for ii = 1:n_satpwr
%         B1        = satpwr_Hz(ii);
%         omega1    = 2*pi*B1;
%         dOmb_rad  = 2*pi*dOmb_Hz;
%         Dw_rad    = 2*pi*ppm_Hz;
%         Dwb_rad   = Dw_rad - dOmb_rad;
%         bpeak = omega1.^2./(omega1.^2+kb.^2+Dwb_rad.^2);
%         Rex = fb.*bpeak.*(kb.*dOmb_rad.^2./(omega1.^2+Dw_rad.^2) + R2b);
%         model(ii,:) = model(ii,:) + cos2th(ii,:).*Rex;
%     end
% end
% end


% % waterRes: water-only Lorentzian residual for the water pre-fit step
% function res = waterRes(p, data_vec, ppm_CEST, cos2th)
% A_w    = p(1);  FWHM_w = p(2);  ofs_w = p(3);
% shape  = A_w * (FWHM_w/2)^2 ./ ((ppm_CEST - ofs_w).^2 + (FWHM_w/2)^2);
% model  = zeros(size(cos2th,1), numel(ppm_CEST));
% for ii = 1:size(cos2th,1)
%     model(ii,:) = cos2th(ii,:) .* shape;
% end
% res = model(:) - data_vec;
% end


% rexModel: evaluates the full Rex+water model, returns (n_satpwr x n_ppm) matrix
function model = rexModel(p, ppm_Hz, satpwr_Hz, cos2th, rexPools, omega_0_MHz)
n_satpwr = numel(satpwr_Hz);
n_ppm    = numel(ppm_Hz);
model    = zeros(n_satpwr, n_ppm);

% Water: R2w*sin^2(theta) model — B1-dependent, one parameter
R2w = p(1);
for ii = 1:n_satpwr
    omega1_w = 2*pi*satpwr_Hz(ii);
    Dw_rad_w = 2*pi*ppm_Hz;
    water_B1 = R2w * omega1_w^2 ./ (omega1_w^2 + Dw_rad_w.^2);
    model(ii,:) = model(ii,:) + cos2th(ii,:) .* water_B1;
end

% Rex contribution per CEST pool (4 params each; pool j starts at p(2+4*(j-1)+1))
for jj = 1:numel(rexPools)
    idx  = 2 + 4*(jj-1) + 1;
    fb   = p(idx);
    kb   = p(idx+1);
    dp   = p(idx+2);
    R2b  = p(idx+3);
    dOmb_Hz = dp * omega_0_MHz;
    Dwb_Hz  = ppm_Hz - dOmb_Hz;
    for ii = 1:n_satpwr
        B1       = satpwr_Hz(ii);
        omega1   = 2*pi*B1;
        Dw_rad   = 2*pi*ppm_Hz;
        Dwb_rad  = 2*pi*Dwb_Hz;
        dOmb_rad = 2*pi*dOmb_Hz;
        bpeak = omega1^2 ./ (omega1^2 + kb^2 + Dwb_rad.^2);
        Rex   = fb .* bpeak .* (kb * dOmb_rad^2 ./ (omega1^2 + Dw_rad.^2) + R2b);
        model(ii,:) = model(ii,:) + cos2th(ii,:) .* Rex;
    end
end
end

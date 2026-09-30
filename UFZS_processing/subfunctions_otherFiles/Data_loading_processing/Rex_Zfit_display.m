% Rex_Zfit_display: Fits CEST peaks using Rex (Zaiss & Bachert 2013, Eq. 23
%   terms 1+2, R2b included, kb>>R2b for Gamma) with residuals computed in
%   Z-spectrum space. The forward model sums all Rex pool contributions and
%   the water (R2w*sin^2(theta)) term in R1*(1/Z-1) space, then converts to
%   Z via Z = R1 / (R1/Z_MT + R2w*sin^2(theta) + sum(Rex)), and minimises
%   against the raw Z-spectral data. MT (if selected) is pre-fitted in
%   Z-space as in the standard 1/Z pipeline.
%
%   INPUTS / OUTPUTS: same as Rex_fit_display
%
function [results, ppars] = Rex_Zfit_display(results, ppars, pflgs, timing)
set(0, 'DefaultFigureWindowStyle', 'docked')

n_satpwr    = size(results.zspec, 1);
omega_0_MHz = results.omega_0_MHz;
satpwr_uT   = results.satT;
satpwr_Hz   = satpwr_uT * 42.577;
ppm         = results.zspecppm;
all_Zspec   = results.zspec;


%% PROMPT FOR T1
answer      = inputdlg({'Water 1H T1 (s):'}, 'Rex (Z-space) fitting: T1', 1, {'4.2'});
results.T1w = str2double(answer{1});
R1          = 1 / results.T1w;


%% GET REX FITTING BOUNDS
pfitvals = setRexBounds;

% Determine pool lists
hasMT     = any(strcmp(ppars.pools, 'MT'));
cestPools = ppars.pools(~strcmp(ppars.pools, 'MT'));
rexPools  = cestPools(~strcmp(cestPools, 'water'));


%% STEP 1 - MT PRE-FITTING (same as standard Rex pipeline)
if hasMT
    disp('***STEP 1: Fitting MT pool using high-ppm wings of the Z-spectrum...')
    dppmWings                   = 2;
    FitParam.MT.ppmExclude      = [round(min(ppm)+dppmWings), round(max(ppm)-dppmWings)];
    FitParam.MT.ppmReinclude    = [0, 0];
    FitParam.MT.wtHigherSatAmpl = false;
    FitParam.MT.fixUndercut     = false;
    FitParam.MT.plot            = true;
    FitParam.MT.lineshape       = ppars.MTlineshape;
    FitParam.R1                 = R1;
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
    disp('No MT pool specified; using flat MT background (Z=1).')
end


%% STEP 2 - 2D JOINT REX FIT (residuals in Z-space)
disp('***STEP 2: 2D joint Rex fit (Z-space residuals) across all offsets and B1 values...')

% ppm fitting range
ppmInclude = [-8, 8];
CESTfitIdx = ppm > ppmInclude(1) & ppm < ppmInclude(2);
ppm_CEST   = ppm(CESTfitIdx);
ppm_Hz     = ppm_CEST * omega_0_MHz;

% Slice raw Z-spectra to fitting range for Z-space residual computation
Zspec_CEST   = all_Zspec(:, CESTfitIdx);
MT_Zfit_CEST = MT_Zfit(:,  CESTfitIdx);

% Build parameter vector: [R2w(1), pool1(4), pool2(4), ...]
%   R2w: water T2 rate [s^-1]; each pool: [fb, kb, delta_ppm, R2b]
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
    @(p) rexZResidual(p, R1, MT_Zfit_CEST, Zspec_CEST, ppm_Hz, satpwr_Hz, rexPools, omega_0_MHz),...
    p0, lb, ub, opts_fit);

CI_all = nlparci(p_opt, residuals, 'jacobian', jacobian);
SS_res = sum(residuals.^2);
SS_tot = sum((Zspec_CEST(:) - mean(Zspec_CEST(:))).^2);
Rsq    = 1 - SS_res / SS_tot;

results.Rex.Rsq = Rsq;
disp(['Rex (Z-space) fit overall R^2 = ' num2str(Rsq,'%0.4f')])


%% EXTRACT RESULTS PER POOL AND REPORT
R2w = p_opt(1);
results.Rex.water.R2w = R2w;
results.Rex.water.CI  = CI_all(1,:);
disp(['  Water: R2w = ' num2str(R2w,'%0.3f') ' s^-1'])

for jj = 1:numel(rexPools)
    name   = rexPools{jj};
    idx    = 2 + 4*(jj-1) + 1;
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

    % ZssPeak for downstream QUESP compatibility
    results.peakfit.(name) = zeros(n_satpwr, 1);
    dOmb_Hz = dp * omega_0_MHz;
    for ii = 1:n_satpwr
        B1       = satpwr_Hz(ii);
        omega1   = 2*pi*B1;
        dOmb_rad = 2*pi*dOmb_Hz;
        bpeak_at_peak = omega1^2 / (omega1^2 + kb^2);
        term1 = fb * kb * dOmb_rad^2 / (omega1^2 + dOmb_rad^2) * bpeak_at_peak;
        term2 = fb * R2b * bpeak_at_peak;
        ZssPeak = R1 / (term1 + term2 + R1);
        results.peakfit.(name)(ii) = 1 - ZssPeak;
    end
end
disp('***STEP 2: COMPLETE!')


%% PLOT 1: Per-B1 docked figures — Z data vs fitted Z model
satpwr_radS = satpwr_uT * 42.577 * 2 * pi;
Z_model_CEST = rexZModel(p_opt, R1, MT_Zfit_CEST, ppm_Hz, satpwr_Hz, rexPools, omega_0_MHz);

for ii = 1:n_satpwr
    figure; h = axes;
    scatter(ppm_CEST, Zspec_CEST(ii,:), 15, [0.5 0.5 0.5], 'filled', 'DisplayName','Raw Z data')
    hold on;
    plot(ppm_CEST, Z_model_CEST(ii,:), 'k-', 'LineWidth', 2, 'DisplayName','Fitted Z')
    plot(ppm_CEST, MT_Zfit_CEST(ii,:), 'b--', 'LineWidth', 1.5, 'DisplayName','MT only')
    set(h,'Xdir','reverse')
    xlabel('\Delta\omega [ppm]'); ylabel('Z')
    legend; grid on; ylim([0 1]);
    title(['Rex (Z-space) fit,  ' num2str(satpwr_uT(ii),'%1.2f') ' \muT (' ...
        num2str(satpwr_radS(ii),'%3.1f') ' rad/s)'])
end


%% PLOT 2: Full reconstructed Z-spectra (undocked)
set(0, 'DefaultFigureWindowStyle', 'normal')
ppm_Hz_full = ppm * omega_0_MHz;
all_Zspec_fit = rexZModel(p_opt, R1, MT_Zfit, ppm_Hz_full, satpwr_Hz, rexPools, omega_0_MHz);
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
title('Reconstructed fitted Z-spectra (Rex Z-space fit)');
ylabel('Z'); xlabel('Offset (ppm)');
ylim([0 1]); set(gca,'XDir','reverse');
end


%% ---- LOCAL HELPER FUNCTIONS ----

% rexZResidual: residuals are Z_model - Z_measured (Z-space)
function res = rexZResidual(p, R1, MT_Zfit_CEST, Zspec_CEST, ppm_Hz, satpwr_Hz, rexPools, omega_0_MHz)
Z_model = rexZModel(p, R1, MT_Zfit_CEST, ppm_Hz, satpwr_Hz, rexPools, omega_0_MHz);
res     = Z_model(:) - Zspec_CEST(:);
end


% rexZModel: forward Z-spectrum model
%   Z = R1 / (R1/Z_MT + R2w*sin^2(theta) + sum_pools(Rex))
function Z = rexZModel(p, R1, MT_Zfit, ppm_Hz, satpwr_Hz, rexPools, omega_0_MHz)
n_satpwr = numel(satpwr_Hz);
n_ppm    = numel(ppm_Hz);
Z        = zeros(n_satpwr, n_ppm);
R2w      = p(1);

for ii = 1:n_satpwr
    B1      = satpwr_Hz(ii);
    omega1  = 2*pi*B1;
    Dw_rad  = 2*pi*ppm_Hz;
    sin2th  = omega1^2 ./ (omega1^2 + Dw_rad.^2);
    cos2th  = Dw_rad.^2  ./ (omega1^2 + Dw_rad.^2);

    % Water direct saturation in R1*(1/Z-1) space
    water_contrib = R2w * sin2th;
%     water_contrib = R2w * omega1^2 ./ (R2w^2 + Dw_rad.^2);  %using Mulkern & Williams, Med Phys 1993

    % Sum Rex pool contributions in R1*(1/Z-1) space
    Rex_total = zeros(1, n_ppm);
    for jj = 1:numel(rexPools)
        idx  = 2 + 4*(jj-1) + 1;
        fb   = p(idx);  kb = p(idx+1);  dp = p(idx+2);  R2b = p(idx+3);
        dOmb_rad = 2*pi * dp * omega_0_MHz;
        Dwb_rad  = Dw_rad - dOmb_rad;
        bpeak    = omega1^2 ./ (omega1^2 + kb^2 + Dwb_rad.^2);
        Rex_pool = fb .* bpeak .* (kb * dOmb_rad^2 ./ (omega1^2 + Dw_rad.^2) + R2b);
%         % to add in cross-term
%         Rex_pool = Rex_pool + fb * kb .* sin2th .* R2b .* (R2b + kb) ./ (omega1^2 + kb^2 + Dwb_rad.^2);
        Rex_total = Rex_total + Rex_pool;
    end

    % Z = R1*cos2th / (R1*cos2th/Z_MT + R2w*sin2th + Rex)  [spin-lock steady state]
    MT_Z_ii = MT_Zfit(ii,:);
    Z(ii,:) = (R1 .* cos2th) ./ (R1 .* cos2th ./ MT_Z_ii + water_contrib + Rex_total);
end
end

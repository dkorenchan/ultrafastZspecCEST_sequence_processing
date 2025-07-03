% zspecPlot: Plots raw spectra, z-spectra, and MTR asymmetry
%
%   INPUTS:
%       results     -   Struct containing processed data to plot along with
%                       other plotting information
%       convZflg    -   Logical: if true, then data are for conventional 
%                       Z-spectroscopy
%       xmax        -   (optional) Maximum ppm value for displaying the MTR
%                       asymmetry plot
%
%   OUTPUTS:    None - results in plotting windows generated
%
function zspecPlot(results,convZflg,xmax)
if nargin<2
    convZflg=false;
    xmax=8;
elseif nargin<3
%     xmax=Inf;
    xmax=8;
end

if convZflg
    xmax=Inf;
end

% Plot all raw spectra on one plot (only if ultrafast Z-spec)
%
if ~convZflg
    figure; hold on
    for i = 1:size(results.spec,1)
        plot(results.specppm,abs(results.spec(i,:)))
    end
    title('Raw spectra')
    ylabel('Signal (arb. units)')
    xlabel('Frequency (ppm)')
    legend(results.speclabels)
    set(gca,'Xdir','reverse')
    axis square
end

% Plot all z-spectra on one plot
%
figure; hold on
for i = 1:size(results.zspec,1)
    plot(results.zspecppm,abs(results.zspec(i,:)))
end
title('z-spectra')
ylabel('MTR')
xlabel('Frequency (ppm)')
legend(results.zspeclabels)
set(gca,'Xdir','reverse')
axis square

% Plot all MTR asymmetry curves on one plot
%
figure; hold on
for i = 1:size(results.zasym,1)
    plot(results.zasymppm,results.zasym(i,:))
end
xlim([0 xmax]);
title('MTR asymmetry')
ylabel('MTR asymmetry')
xlabel('Frequency (ppm)')
legend(results.zspeclabels,'Location','northwest')
set(gca,'Xdir','reverse')
axis square
end
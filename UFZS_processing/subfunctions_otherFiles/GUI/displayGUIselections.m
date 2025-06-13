% displayGUIselections: Generates a text report in the MATLAB command line
% indicating to the user the selections made with the GUI.
%
%   INPUTS:
%       procflgs    -   Struct containing values of logical flags for
%                       data processing options
%       params      -   Struct containing numerical/string processing 
%                       parameters
%
%   OUTPUTS:    NONE - Leads to text printed in the MATLAB command line
%
function displayGUIselections(procflgs,params)
if procflgs.override
    if ~isinf(params.pw90) && ~isnan(params.pw90)
        disp(['Parameter override: 90 pulse width set to ' ...
            num2str(params.pw90) ' us'])
    end
end
if procflgs.topproc 
    disp('Data will be loaded with processing already done in Topspin')
else
    disp('Raw FID data will be loaded and processed')
    if strcmp(params.filter,'exponential')
        disp(['--FID will be weighted with exponential filter, decay rate ' ...
            num2str(params.ap) ' Hz'])
    elseif strcmp(params.filter,'gaussian')
        disp(['--FID will be weighted with Gaussian filter, attenuation ' ...
            num2str(params.edge) ' dB at FID ends'])        
    end
end
if procflgs.norm
    disp('Normalization of raw spectra will be performed:')
    if strcmp(params.normtype,'nosatfit')
        disp(['--Raw spectra will be normalized by determining scaling factors '...
            'using the non-saturated spectra, clustering into two groups, and '...
            'adjusting spectra + calculating z-spectra accordingly'])
    elseif strcmp(params.normtype,'ppmval')
        disp(['--Raw spectra will be normalized to 1st B1 value, using signal at ' ...
            num2str(params.ppmnorm,'%2.2f') ' ppm'])
    end
else
    disp('No spectral normalization across B1 values will be performed')
end
disp(['All z-spectra will be displayed from -' num2str(params.ppmwdw) ' to '...
    num2str(params.ppmwdw) ' ppm'])
if ~procflgs.norm || (procflgs.norm && ~strcmp(params.normtype,'nosatfit'))
    disp('Z-spectra will be calculated using the largest non-saturated spectrum identified')
end
if procflgs.peakfit 
    disp('Peak fitting of z-spectra will be performed for QUESP analysis.')
    disp('--Peak fitting: selected peaks:')
    for i = 1:numel(params.pools)
        disp(params.pools{i});
    end
    if strcmp(params.peaktype,'Pseudo-Voigt')
        disp('--Peak fitting: will use Pseudo-Voigt lineshapes.')
        if procflgs.PVcharconstr
            disp('----Gaussian-Lorentzian peak character will be kept identical for all peaks except water (and MT).')
        end
    elseif strcmp(params.peaktype,'Lorentzian')
        disp('--Peak fitting: will use Lorentzian lineshapes.')
    end
    if procflgs.water1st
        disp(['--Peak fitting: water (+ NOE and/or MT) will be fit first for each power '...
            'using negative ppm values, then all other peaks will be fit'])
    else
        disp('--Peak fitting: all peaks will be fit simultaneously for each power')
    end    
    if procflgs.fix
        disp(['--Peak fitting: will use peak offsets from power index ' ...
            num2str(params.fixind) ' for all other fits'])
    else
        disp('--Peak fitting: peak offsets will be fit independently at each power')
    end
    if procflgs.ppmwt
        disp(['----Fitting will focus on the following ppm value(s): [' ...
            num2str(params.ppmwt) ']'])
    end
    if procflgs.MTsuperLorentz
        disp('--Peak fitting: MT pool will be fit to a super-Lorentzian lineshape')
    end    
else
    disp('MTR asymmetry at user-specified ppm value will be performed for QUESP analysis')
end
end
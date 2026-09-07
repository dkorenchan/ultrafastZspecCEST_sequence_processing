% QUESP_fit_display: Coordinating function to fit either peak-fitted
% amplitudes or MTR asymmetry values to the desired QUESP functions, then
% plot the results
%
%   INPUTS:
%       results     -   Struct containing processed z-spectral data as well
%                       as fitted peak amplitudes (if previously specified)
%                       and other information about the dataset
%       pars        -   Struct containing other processing parameters,
%                       including pool names for QUESP fitting
%       timing      -   Struct containing values of delays in the pulse
%                       sequence:
%                           .tp  - Saturation pulse duration (s)
%                           .rd  - Recovery delay after measurement (s)
%       opts        -   Struct containing user start, upper, and lower
%                       bounds for QUESP fitting
%
%   OUTPUTS:
%       results     -   Struct comprised of the input results but now also
%                       containing results from QUESP fitting
%
function results=QUESP_fit_display(results,ppars,timing,opts)
disp('Performing fitting for exchange and pool parameters using QUESP...')
fitFcns_UFZS={'Regular','Inverse','OmegaPlot'}; %fit to all 3 functions!
scrsz=get(0,'screensize');

% Use pars.pools to ID whether to fit amplitudes from peak fitting, or MTR
% asymmetry
if sum(strcmp(ppars.pools,'MTRasym'))>0
    MTRflg=true;
else
    MTRflg=false;
end

% Prompt for T1 only if not already provided by a prior fitting step
% (e.g. invZ_R1cos2th_fit_display stores it in results.T1w)
if ~isfield(results,'T1w') || isempty(results.T1w) || isnan(results.T1w)
    prompt = {'Water 1H T1 (s):'};
    dlg_title = 'QUESP fitting parameters';
    num_lines = 1;
    def = {'4.2'};
    if MTRflg %also ask for ppm value for MTRasym calculation
        prompt(2) = {'Chemical shift value for QUESP fitting (ppm):'};
        def(2) = {'3.5'};   
    end
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    results.T1w = str2double(answer{1});  
    if MTRflg %also ask for ppm value for MTRasym calculation
        results.QUESPppm = str2double(answer{2});  
        % Obtain z-spectral values corresponding to inputted ppm value
        ppm_error = (results.QUESPppm - results.zasymppm) ./ results.QUESPppm; 
            %use to find index for QUESP fitting
    else
        % Remove water from pools list prior to QUESP fitting
        ppars.pools = ppars.pools(~strcmp(ppars.pools,'water'));
    end
else
    disp(['Using T1w = ' num2str(results.T1w,'%1.3f') ' s from prior 1/Z fitting step.'])
    if MTRflg
        answer_ppm = inputdlg({'Chemical shift value for QUESP fitting (ppm):'}, ...
            'QUESP fitting parameters',1,{'3.5'});
        results.QUESPppm = str2double(answer_ppm{1});
        ppm_error = (results.QUESPppm - results.zasymppm) ./ results.QUESPppm;
    else
        ppars.pools = ppars.pools(~strcmp(ppars.pools,'water'));
    end
end

% Perform QUESP fitting for each pool, treating the 'MTRasym' pool
% differently by using the ppm value signal amplitudes
z_ref = ones(length(results.satT),1);
for i = 1:numel(ppars.pools)
    % Display data points, prompt user to select data points for
    % fitting
    a.(ppars.pools{i})=figure; 
    if strcmp(ppars.pools{i},'MTRasym')
        z_lab=1-results.zasym(:,find(abs(ppm_error)<1e-3));
        if size(z_lab,2) > 1 %catch error: took 2 ppm values
            z_lab = z_lab(:,1);
        end
        qptitle=['QUESP Fitting, ' num2str(results.QUESPppm) ' ppm'];
    else
        z_lab=1-results.peakfit.(ppars.pools{i});
        qptitle=['QUESP Fitting, ' ppars.pools{i} ' ' ppars.peaktype];
    end
    scatter(results.satT,1-z_lab);
    title(ppars.pools{i});
    xlabel('B_1 (\muT)'); ylabel('Peak fit amplitude')
    choices = cell(length(results.satT),1);
    for j = 1:length(results.satT)
        choices{j} = num2str(results.satT(j),'%3.1f');
    end
    prompt = 'Select all saturation powers (in uT) to use for processing:';
    selmod = 'multiple';
    [fitpts.(ppars.pools{i}),fitflg]=listdlg('ListString' , choices , ...
        'SelectionMode' , selmod , ...
        'ListSize' , [300 500] , ...
        'PromptString' , prompt);
    if fitflg
        close(a.(ppars.pools{i}));
        figure('Position',[1 1 scrsz(3) 550]); 
        sgtitle(qptitle,'FontSize',24,'FontWeight','bold')
        pause(0.2); %prevents plotting on top of another subplot!

        results.QUESP.(ppars.pools{i})=QUESPfitting(z_lab,z_ref,results.satT,...
            fitFcns_UFZS,results.T1w,timing,true,fitpts.(ppars.pools{i}),opts);

        % Display goodness-of-fit results:
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp(['Fitting results, ' ppars.pools{i} ':']) 
        disp(['QUESP fitting, goodness-of-fit: ' ...
            num2str(results.QUESP.(ppars.pools{i}).Regular.rsq,'%1.3f')])
        disp(['Inverse QUESP fitting, goodness-of-fit: ' ...
            num2str(results.QUESP.(ppars.pools{i}).Inverse.rsq,'%1.3f')]) 
        disp(['Omega plot fitting, goodness-of-fit: ' ...
            num2str(results.QUESP.(ppars.pools{i}).OmegaPlot.rsq,'%1.3f')]) 
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    else
        disp(['QUESP fitting cancelled by user for ' ppars.pools{i} ' pool.'])            
    end
end
end
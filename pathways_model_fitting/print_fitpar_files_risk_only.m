%% This script is meant to take both constrained/unconstrained fitpar files and print paramatric fit and non-paramatric summary to Excel
% it also create choice matrix for every subject
clearvars
close all

%Input
fitparwave = '07132021';
isconstrained = 0;

% exclude should match those in the fit_parameters.m script
exclude = []; 
% TEMPORARY: subjects incomplete data (that the script is not ready for)

%% folder and subjects
root = 'E:\Ruonan\Pathways2021\data';
data_path = root; 
[subjects, data_all] = getSubjectsData(fullfile(data_path, 'RA_pilot_data_preprocessed.csv'));
subjects = subjects(~ismember(subjects, exclude));

% subjects = [1210];
path = fullfile(root, 'Fitpar files', fitparwave, filesep);
cd(path)

% defining monetary values
valueP = [5 6 7 8 10 12 14 16 19 23 27 31 37 44 52 61 73 86 101 120];

    
    summary_file = ['param_nonparam_' fitparwave '.txt']; % parametric and nonparametric risk and ambig attitudes
%     choiceData_file = [path 'choice_data' outputwave '.xls']; % choice matrix

    % unconstrained
    if isconstrained == 0
        % results file
        fid = fopen([path summary_file],'w')
        fprintf(fid,'\talpha\tgamma\tr2\tLL\tAIC\tBIC\n')
    end


    % Fill in subject numbers separated by commas
    % subjects = {'87','88'};
for s = 1:length(subjects)

    subject = subjects(s); 

    % load gains file for subject and extract params & choice data
    load(['pathways_' num2str(subject) '.mat']);

    aP = model_fit.info.b(2);
    gP = model_fit.info.b(1);
    LLP = model_fit.info.LL;
    r2P = model_fit.info.r2;
    AICP = model_fit.info.AIC;
    BICP = model_fit.info.BIC;
    modelP = model_fit.info.model;
    exitFlagP = model_fit.info.exitflag;
    optimizerP = model_fit.info.optimizer;

    %write into param text file
    fprintf(fid,'%s\t%f\t%f\t%f\t%f\t%f\t%f\n',...
    num2str(subject), aP, gP, r2P, LLP, AICP, BICP);
        

end   
    fclose(fid);




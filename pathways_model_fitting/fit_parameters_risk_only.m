% NOTE: Requires MATLAB optim library
% notice to change the constrained/unconstrained function, and change the fitpat.mat file name to constrained/unconstrained

clearvars
close all
%% pool object 
% poolobj = parpool('local', 8);
%% Input set up
fitparwave = '07132021';
search = 'grid'; % which method for searching optimal parameters
model = 'risk'; % which utility function
isconstrained = 0; % if use constrained fitting. 0-unconstrained, 1-constrained, 2-both

%% Set up loading + subject selection
% TODO: Maybe grab & save condition somewhere?

root = 'E:\Ruonan\Pathways2021\data'; % Need to change if doing analysis in different folders
data_path = root; % root of folders is sufficient
fitpar_out_path = fullfile(root,'Fitpar files', fitparwave);
%graph_out_path  = fullfile(root, 'ChoiceGraphs/');

if exist(fitpar_out_path)==0
    mkdir(fullfile(root,'Fitpar files'),fitparwave)
end

[subjects, data_all] = getSubjectsData(fullfile(data_path, 'RA_pilot_data_preprocessed.csv'));

tic

for subj_idx = 1:length(subjects)

    subjectNum = subjects(subj_idx);
    
%     subjectNum = 3;

    % load single subject's data
    Data = data_all(data_all.id == subjectNum, :);
    
    %% Refine variables
       
    values = Data.val;
    probs  = Data.prob;
    choice = Data.choice;
    ambigs = zeros(size(probs));

    %% Prepare variables for model fitting

    fixed_valueP = 5; % Value of fixed reward
    fixed_prob = 1;   % prb of fixed reward 
    base = 0; % null outcome

    if strcmp(search, 'grid')
    % grid search
    % range of each parameter
        slopeRange = -4:0.2:1;
        aRange = 0.001:0.2:4;
       % three dimenstions
        [b1, b2] = ndgrid(slopeRange, aRange);
        % all posibile combinatinos of three parameters
        b0 = [b1(:) b2(:)];
    elseif strcmp(search,'single')
        % single search
        b0 = [-1 0.5]; % starting point of the search process, [gamma, beta, alpha]
    end


    refVal = fixed_valueP * ones(length(choice), 1);
    refProb = fixed_prob  * ones(length(choice), 1);        

    %% Fit model

    % Two versions of function, calculate both the unconstrained and constrained fittings:
    % fit_ambgiNrisk_model: unconstrained
    if isconstrained == 0 || isconstrained == 2
        [info_uncstr, p_uncstr] = fit_risk_model(choice', ...
            refVal', ...
            values', ...
            refProb', ...
            probs', ...
            ambigs', ...
            model, ...
            b0, ...
            base);

        slope_uncstr = info_uncstr.b(1);
        a_uncstr = info_uncstr.b(2);
        r2_uncstr = info_uncstr.r2;

        disp(['Subject ' num2str(subjectNum) ' unconstrained fitting completed'])

    end

    if isconstrained == 1 || isconstrained == 2
        % fit_ambigNrisk_model_Constrained: constrained on alpha and beta    
        [info_cstr, p_cstr] = fit_risk_model_Constrained(choice', ...
            refVal', ...
            values', ...
            refProb', ...
            probs', ...
            ambigs', ...
            model, ...
            b0, ...
            base);

        slope_cstr = info_cstr.b(1);
        a_cstr = info_cstr.b(2);
        r2_cstr = info_cstr.r2;

        disp(['Subject ' num2str(subjectNum) ' constrained fitting completed'])
    end

    %% Create choice matrices

    % One matrix per condition. Matrix values are binary (0 for sure
    % choice, 1 for lottery). Matrix dimensions are prob/ambig-level
    % x payoff values. Used for graphing and some Excel exports.

    % Outputs:
    %  riskyChoicesP

    % Create riskyChoicesP
    % Risk levels by payoff values
    valueP = unique(values(ambigs == 0 & values ~= 4));
    prob = unique(probs);
    
    riskyChoicesP = zeros(length(prob), length(valueP));
    for i = 1:length(prob)
        for j = 1:length(valueP)
            selection = find(probs == prob(i) & values == valueP(j) & ambigs == 0);
            if ~isempty(selection)
                riskyChoicesP(i, j) = choice(selection);
            else
                riskyChoicesP(i, j) = NaN;
            end
        end
    end
    
    % Creat risky/ambig choiecs by level (nonparametric), excluding the value 5
    riskyChoices_byLevel = zeros(1, length(prob));
    for i=1:length(prob)
        riskyChoices_byLevel(1,i) = nanmean(riskyChoicesP(i,2:length(riskyChoicesP)));
    end 

    %% Graph
%         colors =   [255 0 0;
%         180 0 0;
%         130 0 0;
%         52 181 233;
%         7 137 247;
%         3 85 155;
%         ]/255;
% 
%         figure
%         counter=5;
%         for i=1:3
%             subplot(3,2,counter)
%             plot(valueP,ambigChoicesP(i,:),'--*','Color',colors(3+i,:))
%             legend([num2str(ambig(i)) ' ambiguity'])
%             if counter==1
%                 title(['Beta = ' num2str(bP)])
%             end
%             ylabel('Chose Lottery')
%             if counter==5
%             xlabel('Lottery Value ($)')
%             end
%             counter=counter-2;
%         end
% 
%         counter=2;
%         for i=1:3
%             subplot(3,2,counter)
%             plot(valueP,riskyChoicesP(i,:),'--*','Color',colors(i,:))
%             legend([num2str(prob(i)) ' probability'])
%             if counter==2
%                 title(['Alpha = ' num2str(aP)])
%             end
%                 if counter==6
%             xlabel('Lottery Value ($)')
%                 end
%             counter=counter+2;
%         end
% 
%         set(gcf,'color','w');
%         figName=['RA_GAINS_' num2str(subjectNum) '_fitpar'];
%          exportfig(gcf,figName,'Format','eps','bounds','tight','color','rgb','LockAxes',1,'FontMode','scaled','FontSize',1,'Width',4,'Height',2,'Reference',gca);
% 
%        %% figure with fitted logistic line, only for gain
%         xP = 0:0.1:max(valueP);
%         uFP = fixed_prob * (fixed_valueP).^aP;
% 
%        figure
% 
%         % risk pos
%         for i = 1 :length(prob)
%             plot(valueP,riskyChoicesP(i,:),'o','MarkerSize',8,'MarkerEdgeColor',colors([1 1 1])...
%                 ,'MarkerFaceColor',colors(i,:),'Color',colors(i,:));
%               hold on
%             % logistic function
%             uA = prob(i) * xP.^aP;
%             p = 1 ./ (1 + exp(slopeP*(uA-uFP)));
% 
%             plot(xP,p,'-','LineWidth',4,'Color',colors(i,:));
%             axis([0 150 0 1])
%             set(gca, 'ytick', [0 0.5 1])
%             set(gca,'xtick', [0 20  40  60  80  100 120])
%             set(gca,'FontSize',25)
%             set(gca,'LineWidth',3)
%             set(gca, 'Box','off')
% 
% 
%         end
%         title(['  alpha gain = ' num2str(aP)]);
% 
%         figure
%         % ambig pos
%         for i = 1:length(ambig)
%             plot(valueP,ambigChoicesP(i,:),'o','MarkerSize',8,'MarkerEdgeColor',colors([1 1 1]),'MarkerFaceColor',colors(length(prob)+i,:));
%              hold on
%  
%             % logistic function
%             uA = (0.5 - bP.*ambig(i)./2) * xP.^aP;
%             p = 1 ./ (1 + exp(slopeP*(uA-uFP)));
% 
% 
%             plot(xP,p,'-','LineWidth',2,'Color',colors(length(prob)+i,:));
%             set(gca, 'ytick', [0 0.5 1])
%             set(gca,'xtick', [0 20  40  60  80  100 120])
%             set(gca,'FontSize',25)
%             set(gca,'LineWidth',3)
%             set(gca, 'Box','off')
% 
%         end
%         title([ '  beta gain = ' num2str(bP)]);
% 
%         %     % risk neg
%         subplot(2,2,2)
%             for i = 1:length(prob)
%             plot(valueN,riskyChoicesN(i,:),'o','MarkerSize',8,'MarkerEdgeColor',colors([1 1 1])...
%                 ,'MarkerFaceColor',colors(i,:),'Color',colors(i,:));
%             hold on
%     
%             % logistic function
%             uA = -prob(i) * (-xN).^aN;
%             p = 1 ./ (1 + exp(slopeN*(uA-uFN)));
%     
%             plot(xN,p,'-','LineWidth',2,'Color',colors(i,:));
%     
%             end
%         title([char(subject) '  alpha loss = ' num2str(aN)]);
%             
%       
%         % ambig neg
%         subplot(2,2,4)
%         for i = 1:length(ambig)
%             plot(valueN,ambigChoicesN(i,:),'o','MarkerSize',8,'MarkerEdgeColor',colors([1 1 1]),'MarkerFaceColor',colors(length(prob)+i,:));
%               hold on
%     
%             % logistic function
%             uA = -(0.5 - bN.*ambig(i)./2) * (-xN).^aN;
%             p = 1 ./ (1 + exp(slopeN*(uA-uFN)));
%     
%     
%             plot(xN,p,'-','LineWidth',2,'Color',colors(length(prob)+i,:));
%     
%         end
%         title([char(subject) '  beta loss = ' num2str(bN)]);       

    %% Save generated values
   model_fit = struct('id', subjectNum,...
        'info', info_uncstr,...
        'p', p_uncstr,...
        'riskyChoices', riskyChoicesP,...
        'search', search);
    
    save_mat(model_fit, subjectNum, fitpar_out_path)

end

toc

% delete(poolobj)
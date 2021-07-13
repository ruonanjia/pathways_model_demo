clearvars
% close all

%% generate parameters
alphas = [0.3:0.3:1.2]; % risk attitude
gamma = -1.5; % noise of choice

%% Problem set
SAFE_REWARD = 5;

PROBABILITIES = [25, 50, 75]./100;
RISKY_REWARDS = [5, 6, 7, 8, 10, 12, 14, 16, 19, 23, 27, 31, 37, 44, 52, 61, 73, 86, 101, 120];
[P, V] = meshgrid(PROBABILITIES,RISKY_REWARDS); % RISKY problem set
A = zeros(size(P));
n_trials = numel(P);

dimP = size(P);
probs = [reshape(P, dimP(1)*dimP(2), [])];
dimA = size(A);
ambigs = [zeros(dimA(1)*dimA(2), 1)];
values = [reshape(V, dimP(1)*dimP(2), [])];

refVal = 5*ones(size(values));
refProb = ones(size(probs));
model = 'risk';
% search_method = 'grid';
base = 0;

for i = 1:length(alphas)
    alpha = alphas(i);
    %% generate choice
    choice_prob = choice_prob_ambigNrisk(base,refVal,values,refProb,probs,ambigs,[gamma, alpha],model);

    choices = binornd(1,choice_prob);

    % choice matrix
    valueP = unique(values(ambigs == 0 & values ~= 4));
    prob = unique(probs);

    riskyChoicesP = zeros(length(prob), length(valueP));
    for i = 1:length(prob)
        for j = 1:length(valueP)
            selection = find(probs == prob(i) & values == valueP(j) & ambigs == 0);
            if ~isempty(selection)
                riskyChoicesP(i, j) = choices(selection);
            else
                riskyChoicesP(i, j) = NaN;
            end
        end
    end

    %% plot
    colors =   [255 0 0;
        180 0 0;
        130 0 0;
        ]/255;

    xP = 0:0.1:max(values); % x axis
    uFP = 1 * 5.^alpha; %SV of the safe choice

    figure

    for i = 1:length(prob)
         % logistic function
        uA = prob(i) * xP.^alpha;
        p = 1 ./ (1 + exp(gamma*(uA-uFP)));

        plot(xP,p,'-','LineWidth',4,'Color',colors(i,:));
        hold on

        axis([0 130 0 1])
        set(gca, 'ytick', [0 0.25 0.5 0.75 1])
        set(gca,'xtick', [0 20  40  60  80  100 120])
        set(gca,'FontSize',25)
        set(gca,'LineWidth',3)
        set(gca, 'Box','off')

    end

    for i = 1 :length(prob)
        plot(unique(values),riskyChoicesP(i,:),'o','MarkerSize',8,'MarkerEdgeColor',colors([1 1 1])...
            ,'MarkerFaceColor',colors(i,:),'Color',colors(i,:));
          hold on      
    end

    legend([num2str(prob(1)*100) '%'], [num2str(prob(2)*100) '%'], [num2str(prob(3)*100) '%'], 'Location','southeast')

    title(['  alpha = ' num2str(alpha)]);
end
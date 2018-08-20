%% Dissertation Chapter 4, Figure 1 - Final code used
% "Stage 46 tadpoles have more responding ROIs, and these ROIs respond more
% often but with similar peak dF/F0."
% "A. Proportion of cells that respond in an experiment. p(KS) = 0.0294,
% median(46)  = 0.538  median(49) = 0.317, N(46) = 9, N(49) = 18. B. 
% Proportion of trials with a response (peak dF/F0 > 0.1) for all ROIs that
% respond at least once in the experiment. p(KS) < 0.001, median(46)  = 0.061  
% median(49) = 0.0.042  N(46) = 9, n(46) = 475, N(49) = 18, n(49) = 527. 
% C. Average peak dF/F0 for all trials with a response (peak dF/F0 > 0.1) 
% over all cells, split by stage. Average peak dF/F0 in response to 
% Multisensory (MS) and Mechanosensory (M) distributions are not different 
% across stage but Visual (V) is. N(46) = 9, n(46) = 475, N(49) = 18, 
% n(49) = 527. MS: p(KS) = 0.89, median(46)  = 0.195  median(49) = 0.196;  
% V*: p(KS) = 0.025 median(46)  = 0.185  median(49) = 0.198; M: p(KS) = 0.12, 
% median(46)  = 0.198  median(49) = 0.196. Insets are the extra points 
% beyond 2 (the end of the X-axis)."

% This is the code to generate the figure from my disseration. It's
% copy/paste from the scripts actually used for organizational purposes.
% Begin with file [NAME] which contains the final tadpole{1,:} smoothed
% stimuli for all tads used. 

%% Figure 1A: Proportion of cells that respond in an exp 46v49



%% Figure 1B: Proportion of trials with a response across all respROIs
% from diss_fig3, 3B

prop_resp_allrespROI = [];
for t = 1:length(allData)
        stage = allData{1,t}.stage * ones(length(allData{1,t}.resp_ROIs), 1);
        data = allData{1,t}.sum_responses(allData{1,t}.resp_ROIs) / length(allData{1,t}.stimorder);
        prop_resp_allrespROI = [prop_resp_allrespROI; data, stage];
end

%plot(prop_resp_allrespROI(:,2), prop_resp_allrespROI(:,1), 'o')

x_vals = rand(1, 1064) *0.1;
s46_ids = find(prop_resp_allrespROI(:,2) == 46)
s49_ids = find(prop_resp_allrespROI(:,2) == 49)
s46_propresp = prop_resp_allrespROI(s46_ids, 1)
s49_propresp = prop_resp_allrespROI(s49_ids, 1)

% stat test
% Wilcoxen rank sum test
% (https://www.mathworks.com/help/stats/ranksum.html)

[p, h, states] = ranksum(s46_propresp, s49_propresp)
% p < 0.01

% error bars
s46_mean = mean(s46_propresp)
s49_mean = mean(s49_propresp)
s46_std = std(s46_propresp) / sqrt(length(s46_propresp))
s49_std = std(s49_propresp) / sqrt(length(s49_propresp))

% Make plot
figure;
hold on
plot(x_vals(1:475), s46_propresp, 'go')
plot(x_vals(1:527)+0.5, s49_propresp, 'mo')
errorbar([0.05, .55], [s46_mean, s49_mean], [s46_std, s49_std], 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
e.Color = 'k'
hold off 
xlim([-0.2, 0.8])
ylabel('proportion responses')
xlabel('stage')
ax = gca;
ax.XTick = [0.05 0.55];
ax.XTickLabel = [46 49];
set(gca, 'Fontsize', 20);
saveas(gcf, 'prop respond by respROI by stage', 'png')
saveas(gcf, 'prop respond by respROI by stage', 'epsc2')

%% Figure 1C, 
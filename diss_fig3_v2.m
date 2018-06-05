%% Dissertation figure 3: 46 vs 49 diffs
% this will actually be 3 figures: 
    % fig 2 = prop resp rois, prop responses, avg peak resp
    % fig 3 = avg onset time of response
    % fig 4 = multisensory index by peak and timing assessment (MSIndex
    % onset time, but don't call it that)
% using corrected for bad trials data 
% use xcorr_cleaned_20180313 file 
%load('D:\Torrey_calcium_imaging\compare_46-49\analysis_Feb 2018\corrected_for_badtrials\46_49_comparison\xcorr_corrected_20180313.mat')

% steal figures from compare_46_49_new, but make them pretty

%% 2A prop respROIs
for t = 1:length(allData)
    key(t) = allData{1,t}.stage;
end

st46 = find(key == 46);
st49 = find(key == 49);
% generate x values for plotting purposes 
st46_X = ones(length(st46),1);
st49_X = 1.5*ones(length(st49),1);

for t = 1:length(allData)
    ct_respROIs(t,1) = length(allData{1,t}.resp_ROIs);
    ct_respROIs(t,2) = size(allData{1,t}.ROIcenters,1);
    ct_respROIs(t,3) = ct_respROIs(t,1) / ct_respROIs(t,2);
end

% get mean and SD of proportion of respROIs
stats_propresp(1,1) = median(ct_respROIs(st46,3))
stats_propresp(2,1) = median(ct_respROIs(st49,3))
stats_propresp(1,2) = std(ct_respROIs(st46,3))
stats_propresp(2,2) = std(ct_respROIs(st49,3))
[h, p] = kstest2(ct_respROIs(st46,3), ct_respROIs(st49,3))
% h = 1, p = 0.0294

% scatterplot the values by stage - proportion of respROIs
figure;
hold on
scatter(st46_X, ct_respROIs(st46,3), 80, 'g', 'LineWidth', 2)
scatter(st49_X, ct_respROIs(st49,3), 80, 'm', 'LineWidth', 2)
plot([0.9 1.1], [stats_propresp(1,1), stats_propresp(1,1)], 'k', 'LineWidth', 6)
plot([1.4 1.6], [stats_propresp(2,1), stats_propresp(2,1)], 'k', 'LineWidth', 6)
hold off
xlim([0.5 2])
ylim([0 1.1])
ax = gca;
ax.XTick = [1 1.5];
ax.XTickLabel = [46 49];
%title('Proportion of Responding ROIs by Exp')
xlabel('Stage')
ylabel('Proportion of ROIs')
set(gca, 'FontSize', 30)
saveas(gcf, 'prop resROIs by stage', 'epsc2')

%% 2B: proportion responses in each respROI

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
[h p] = kstest2(s46_propresp, s49_propresp)
% p < 0.001

% error bars
s46_mean = median(s46_propresp)
s49_mean = median(s49_propresp)
s46_std = std(s46_propresp) / sqrt(length(s46_propresp))
s49_std = std(s49_propresp) / sqrt(length(s49_propresp))

% Make plot
figure;
hold on
plot(x_vals(1:475), s46_propresp, 'go')
plot(x_vals(1:527)+0.5, s49_propresp, 'mo')
plot([-0.05 0.15], [s46_mean s46_mean], 'k', 'LineWidth', 6)
plot([0.45 0.65], [s49_mean s49_mean], 'k', 'LineWidth', 6)
%errorbar([0.05, .55], [s46_mean, s49_mean], [s46_std, s49_std], 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
%e.Color = 'k'
hold off 
xlim([-0.1, 0.7])
ylabel('proportion responses')
xlabel('stage')
ax = gca;
ax.XTick = [0.05 0.55];
ax.XTickLabel = [46 49];
set(gca, 'Fontsize', 30);
%saveas(gcf, 'prop respond by respROI by stage', 'png')
%saveas(gcf, 'prop respond by respROI by stage', 'epsc2')


%% 2C Average peak response

% Add to allRespROIs with the data for these plots
idx = 1

for t = 1:length(allData)
    for roi = 1:length(allData{1,t}.resp_ROIs)
        r = allData{1,t}.resp_ROIs(roi);
        allRespROIs(idx, 14:15) = allData{1,t}.avg_peak2(2:3, r);
        allRespROIs(idx, 16:17) = allData{1,t}.avg_onsettime2(2:3, r);
        allRespROIs(idx, 18:19) = allData{1,t}.std_onsettime2(2:3, r);
        idx = idx +1;
    end
end

% identify rows with each stage
st46 = find(allRespROIs(:,3) == 46);
st49 = find(allRespROIs(:,3) == 49);

% MS avg peak response
% ECDF
[f1, x1] = ecdf(allRespROIs(st46, 8));
[f2, x2] = ecdf(allRespROIs(st49, 8));
figure;
hold on
plot(x1, f1, 'g', 'LineWidth', 1) %st 46 = green
plot(x2, f2, 'm', 'LineWidth', 1) %st 49 = purple
legend({'stage 46', 'stage 49'})
hold off
ylabel('ROIs')
xlabel('Mean Peak \DeltaF/F_{0}')
xlim([0 2])
set(gca,'FontSize',30)
saveas(gcf, 'ecdf of MS mean peak', 'png')
saveas(gcf, 'ecdf of MS mean peak', 'epsc2')

% 46 hist
figure;
histogram(allRespROIs(st46, 8), 40, 'FaceColor', 'g', 'EdgeColor', 'g', 'BinLimits', [0.1 2]) 
ylim([0 130])
ylabel('ROIs')
xlabel('Mean Peak \DeltaF/F_{0}')
%xlim([0 2])
set(gca,'FontSize',30)
saveas(gcf, 'hist of MS mean peak 46', 'png')
saveas(gcf, 'hist of MS mean peak 46', 'epsc2')

% 49 hist
figure;
histogram(allRespROIs(st49, 8), 40, 'FaceColor', 'm', 'EdgeColor', 'm', 'BinLimits', [0.1 2]) 
ylabel('ROIs')
xlabel('Mean Peak \DeltaF/F_{0}')
ylim([0 130])
%xlim([0 2])
set(gca,'FontSize',30)
saveas(gcf, 'hist of MS mean peak 49', 'png')
saveas(gcf, 'hist of MS mean peak 49', 'epsc2')

% V avg peak response
% ECDF
[f1, x1] = ecdf(allRespROIs(st46, 14));
[f2, x2] = ecdf(allRespROIs(st49, 14));
figure;
hold on
plot(x1, f1, 'g', 'LineWidth', 1) %st 46 = green
plot(x2, f2, 'm', 'LineWidth', 1) %st 49 = purple
legend({'stage 46', 'stage 49'})
hold off
ylabel('ROIs')
xlabel('Mean Peak \DeltaF/F_{0}')
xlim([0 2])
set(gca,'FontSize',30)
saveas(gcf, 'ecdf of V mean peak', 'png')
saveas(gcf, 'ecdf of V mean peak', 'epsc2')

% 46 hist
figure;
histogram(allRespROIs(st46, 14), 40, 'FaceColor', 'g', 'EdgeColor', 'g', 'BinLimits', [0.1 2]) 
ylim([0 130])
ylabel('ROIs')
xlabel('Mean Peak \DeltaF/F_{0}')
%xlim([0 2])
set(gca,'FontSize',30)
saveas(gcf, 'hist of V mean peak 46', 'png')
saveas(gcf, 'hist of V mean peak 46', 'epsc2')

% 49 hist
figure;
histogram(allRespROIs(st49, 14), 40, 'FaceColor', 'm', 'EdgeColor', 'm', 'BinLimits', [0.1 2]) 
ylabel('ROIs')
xlabel('Mean Peak \DeltaF/F_{0}')
ylim([0 130])
%xlim([0 2])
set(gca,'FontSize',30)
saveas(gcf, 'hist of V mean peak 49', 'png')
saveas(gcf, 'hist of V mean peak 49', 'epsc2')

% M avg peak response
% ECDF
[f1, x1] = ecdf(allRespROIs(st46, 15));
[f2, x2] = ecdf(allRespROIs(st49, 15));
figure;
hold on
plot(x1, f1, 'g', 'LineWidth', 1) %st 46 = green
plot(x2, f2, 'm', 'LineWidth', 1) %st 49 = purple
legend({'stage 46', 'stage 49'})
hold off
ylabel('ROIs')
xlabel('Mean Peak \DeltaF/F_{0}')
xlim([0 2])
set(gca,'FontSize',30)
saveas(gcf, 'ecdf of M mean peak', 'png')
saveas(gcf, 'ecdf of M mean peak', 'epsc2')

% 46 hist
figure;
histogram(allRespROIs(st46, 15), 40, 'FaceColor', 'g', 'EdgeColor', 'g', 'BinLimits', [0.1 2]) 
ylim([0 130])
ylabel('ROIs')
xlabel('Mean Peak \DeltaF/F_{0}')
%xlim([0 2])
set(gca,'FontSize',30)
saveas(gcf, 'hist of M mean peak 46', 'png')
saveas(gcf, 'hist of M mean peak 46', 'epsc2')

% 49 hist
figure;
histogram(allRespROIs(st49, 15), 40, 'FaceColor', 'm', 'EdgeColor', 'm', 'BinLimits', [0.1 2]) 
ylabel('ROIs')
xlabel('Mean Peak \DeltaF/F_{0}')
ylim([0 130])
%xlim([0 2])
set(gca,'FontSize',30)
saveas(gcf, 'hist of M mean peak 49', 'png')
saveas(gcf, 'hist of M mean peak 49', 'epsc2')

%stats
vars = [8 14 15]
for v = 1:length(vars)
    [h(v) p(v)] = kstest2(allRespROIs(st49, vars(v)), allRespROIs(st46, vars(v)))
    peak_median(v,1) = nanmedian(allRespROIs(st49, vars(v)));
    peak_median(v,2) = nanmedian(allRespROIs(st46, vars(v)));
end

% Find all values larger than 2
for v = 1:length(vars)
    tmp = allRespROIs(st46, vars(v)) 
    tmp2 = tmp > 2
    large_peaks{v,1} = tmp(tmp2)
    clear('tmp', 'tmp2')
    tmp = allRespROIs(st49, vars(v)) 
    tmp2 = tmp > 2
    large_peaks{v,2} = tmp(tmp2)    
    clear('tmp', 'tmp2')
end

% plot all values larger than 2 by stage for each stim
xvals = ones(1, 7)
for v = 1:length(vars)
    figure;
    hold on
    plot(xvals(1:length(large_peaks{v,1})), large_peaks{v,1}, 'go', 'LineWidth', 3)
    plot((xvals(1:length(large_peaks{v,2}))+0.5), large_peaks{v,2}, 'mo', 'LineWidth', 3)
    hold off
    ylabel('Mean Peak \DeltaF/F_{0}')
    xlim([0.5 2])
    %ylim([0 1.1])
    ax = gca;
    ax.XTick = [1 1.5];
    ax.XTickLabel = [46 49];
    title(sprintf('Var %d', vars(v)))
    set(gca, 'FontSize', 30)
end

%% Figure 3: avg onset time
% data stored in allRespROIs, MS=9 V=16 M=17

% MS avg onset time
% ECDF
[f1, x1] = ecdf(allRespROIs(st46, 9));
[f2, x2] = ecdf(allRespROIs(st49, 9));
figure;
hold on
plot(x1, f1, 'g', 'LineWidth', 1) %st 46 = green
plot(x2, f2, 'm', 'LineWidth', 1) %st 49 = purple
%legend({'stage 46', 'stage 49'})
hold off
ylabel('ROIs')
xlabel('Mean onset time (sec)')
xlim([0 7])
set(gca,'FontSize',30)
saveas(gcf, 'ecdf of MS mean onset time', 'png')
saveas(gcf, 'ecdf of MS mean onset time', 'epsc2')

% 46 hist
figure;
histogram(allRespROIs(st46, 9), 70, 'FaceColor', 'g', 'EdgeColor', 'g', 'BinLimits', [0.1 7]) 
%ylim([0 130])
ylabel('ROIs')
%xlabel('Mean Peak \DeltaF/F_{0}')
xlim([0 7])
set(gca,'FontSize',30)
saveas(gcf, 'hist of MS mean onset time 46', 'png')
saveas(gcf, 'hist of MS mean onset time 46', 'epsc2')

% 49 hist
figure;
histogram(allRespROIs(st49, 9), 70, 'FaceColor', 'm', 'EdgeColor', 'm', 'BinLimits', [0.1 7]) 
ylabel('ROIs')
xlabel('Mean Peak \DeltaF/F_{0}')
%ylim([0 220])
xlim([0 7])
set(gca,'FontSize',30)
saveas(gcf, 'hist of MS mean onset time 49', 'png')
saveas(gcf, 'hist of MS mean onset time 49', 'epsc2')

% V avg onset time
% ECDF
[f1, x1] = ecdf(allRespROIs(st46, 16));
[f2, x2] = ecdf(allRespROIs(st49, 16));
figure;
hold on
plot(x1, f1, 'g', 'LineWidth', 1) %st 46 = green
plot(x2, f2, 'm', 'LineWidth', 1) %st 49 = purple
%legend({'stage 46', 'stage 49'})
hold off
%ylabel('ROIs')
xlabel('Mean onset time (sec)')
xlim([0 7])
set(gca,'FontSize',30)
saveas(gcf, 'ecdf of V mean onset time', 'png')
saveas(gcf, 'ecdf of V mean onset time', 'epsc2')

% 46 hist
figure;
histogram(allRespROIs(st46, 16), 70, 'FaceColor', 'g', 'EdgeColor', 'g', 'BinLimits', [0.1 7]) 
%ylim([0 130])
%ylabel('ROIs')
%xlabel('Mean Peak \DeltaF/F_{0}')
xlim([0 7])
set(gca,'FontSize',30)
saveas(gcf, 'hist of V mean onset time 46', 'png')
saveas(gcf, 'hist of V mean onset time 46', 'epsc2')

% 49 hist
figure;
histogram(allRespROIs(st49, 16), 70, 'FaceColor', 'm', 'EdgeColor', 'm', 'BinLimits', [0.1 7]) 
%ylabel('ROIs')
xlabel('Mean onset time (sec)')
%ylim([0 130])
xlim([0 7])
set(gca,'FontSize',30)
saveas(gcf, 'hist of V mean onset time 49', 'png')
saveas(gcf, 'hist of V mean onset time 49', 'epsc2')

% M avg onset time
% ECDF
[f1, x1] = ecdf(allRespROIs(st46, 17));
[f2, x2] = ecdf(allRespROIs(st49, 17));
figure;
hold on
plot(x1, f1, 'g', 'LineWidth', 1) %st 46 = green
plot(x2, f2, 'm', 'LineWidth', 1) %st 49 = purple
%legend({'stage 46', 'stage 49'})
hold off
ylabel('ROIs')
xlabel('Mean onset time (sec)')
xlim([0 7])
set(gca,'FontSize',30)
saveas(gcf, 'ecdf of M mean onset time', 'png')
saveas(gcf, 'ecdf of M mean onset time', 'epsc2')

% 46 hist
figure;
histogram(allRespROIs(st46, 17), 70, 'FaceColor', 'g', 'EdgeColor', 'g', 'BinLimits', [0.1 7]) 
%ylim([0 130])
%ylabel('ROIs')
%xlabel('Mean Peak \DeltaF/F_{0}')
xlim([0 7])
set(gca,'FontSize',30)
saveas(gcf, 'hist of M mean onset time 46', 'png')
saveas(gcf, 'hist of M mean onset time 46', 'epsc2')

% 49 hist
figure;
histogram(allRespROIs(st49, 17), 70, 'FaceColor', 'm', 'EdgeColor', 'm', 'BinLimits', [0.1 7]) 
%ylabel('ROIs')
xlabel('Mean onset time (sec)')
%ylim([0 130])
xlim([0 7])
set(gca,'FontSize',30)
saveas(gcf, 'hist of M mean onset time 49', 'png')
saveas(gcf, 'hist of M mean onset time 49', 'epsc2')

%stats
vars = [9 16 17]
for v = 1:length(vars)
    [h(v) p(v)] = kstest2(allRespROIs(st49, vars(v)), allRespROIs(st46, vars(v)))
    onsettime_median(v,1) = nanmedian(allRespROIs(st49, vars(v)));
    onsettime_median(v,2) = nanmedian(allRespROIs(st46, vars(v)));
end
% all signficant, p < 0.001

% doesn't make sense to do a cut off for this one. 


%% Figure 4: Multisensory assessment

%% A: by avg peak (in allRespROIs 11)
% ECDF
[f1, x1] = ecdf(allRespROIs(st46, 11));
[f2, x2] = ecdf(allRespROIs(st49, 11));
figure;
hold on
plot(x1, f1, 'g', 'LineWidth', 1) %st 46 = green
plot(x2, f2, 'm', 'LineWidth', 1) %st 49 = purple
%legend({'stage 46', 'stage 49'})
hold off
ylabel('ROIs')
xlabel('MS Index')
xlim([-inf 2])
set(gca,'FontSize',30)
saveas(gcf, 'ecdf of MSind peak', 'png')
saveas(gcf, 'ecdf of MSInd peak', 'epsc2')

% 46 hist
figure;
histogram(allRespROIs(st46, 11), 60, 'FaceColor', 'g', 'EdgeColor', 'g' , 'BinLimits', [-1 2]) 
%ylim([0 130])
ylabel('ROIs')
xlabel('MS Index')
%xlim([0 7])
set(gca,'FontSize',30)
saveas(gcf, 'hist of MSInd peak 46', 'png')
saveas(gcf, 'hist of MSInd peak 46', 'epsc2')

% 49 hist
figure;
histogram(allRespROIs(st49, 11), 60, 'FaceColor', 'm', 'EdgeColor', 'm', 'BinLimits', [-1 2]) 
ylabel('ROIs')
xlabel('MS Index')
%ylim([0 130])
%xlim([0 7])
set(gca,'FontSize',30)
saveas(gcf, 'hist of MSInd peak 49', 'png')
saveas(gcf, 'hist of MSInd peak 49', 'epsc2')


%% B: Timing Enhancement

% ECDF
[f1, x1] = ecdf(allRespROIs(st46, 12));
[f2, x2] = ecdf(allRespROIs(st49, 12));
figure;
hold on
plot(x1, f1, 'g', 'LineWidth', 1) %st 46 = green
plot(x2, f2, 'm', 'LineWidth', 1) %st 49 = purple
%legend({'stage 46', 'stage 49'})
hold off
ylabel('ROIs')
xlabel('Ratio')
xlim([-inf 3])
set(gca,'FontSize',30)
saveas(gcf, 'ecdf of MSind onset time', 'png')
saveas(gcf, 'ecdf of MSInd onset time', 'epsc2')

% 46 hist
figure;
histogram(allRespROIs(st46, 12), 60, 'FaceColor', 'g', 'EdgeColor', 'g' , 'BinLimits', [-1 3]) 
%ylim([0 130])
ylabel('ROIs')
xlabel('Ratio')
%xlim([0 7])
set(gca,'FontSize',30)
saveas(gcf, 'hist of MSInd onset time 46', 'png')
saveas(gcf, 'hist of MSInd onset time 46', 'epsc2')

% 49 hist
figure;
histogram(allRespROIs(st49, 12), 60, 'FaceColor', 'm', 'EdgeColor', 'm', 'BinLimits', [-1 3]) 
ylabel('ROIs')
xlabel('Ratio')
%ylim([0 130])
%xlim([0 7])
set(gca,'FontSize',30)
saveas(gcf, 'hist of MSInd onset time  49', 'png')
saveas(gcf, 'hist of MSInd onset time  49', 'epsc2')



%% stats
vars = [11, 12]
for v = 1:length(vars)
    [h(v) p(v)] = kstest2(allRespROIs(st49, vars(v)), allRespROIs(st46, vars(v)))
    peak_median(v,1) = nanmedian(allRespROIs(st49, vars(v)));
    peak_median(v,2) = nanmedian(allRespROIs(st46, vars(v)));
end

%% Large Values
%Find all values larger than 3
vals = [2 3];
for v = 1:length(vars)
    tmp = allRespROIs(st46, vars(v)) 
    tmp2 = tmp > vals(v);
    large_peaks{v,1} = tmp(tmp2);
    clear('tmp', 'tmp2')
    tmp = allRespROIs(st49, vars(v)) ;
    tmp2 = tmp > vals(v);
    large_peaks{v,2} = tmp(tmp2); 
    clear('tmp', 'tmp2')
end

% plot all values larger than 2 by stage for each stim
xvals = ones(1, 29)
for v = 1:length(vars)
    figure;
    hold on
    plot(xvals(1:length(large_peaks{v,1})), large_peaks{v,1}, 'go', 'LineWidth', 3)
    plot((xvals(1:length(large_peaks{v,2}))+0.5), large_peaks{v,2}, 'mo', 'LineWidth', 3)
    hold off
    ylabel('MS Index')
    xlim([0.5 2])
    %ylim([0 1.1])
    ax = gca;
    ax.XTick = [1 1.5];
    ax.XTickLabel = [46 49];
    title(sprintf('Var %d', vars(v)))
    set(gca, 'FontSize', 30)
end

ylabel('Ratio')


%% C Multisensory Index by number of responses (response reliability)

% metric: Contrast Index = (MS - MaxUni) / (MS + MaxUni)

for t = 1:length(allData)
    MS_trials = find(allData{1,t}.stimorder == 1);
    length(MS_trials)
    V_trials = find(allData{1,t}.stimorder == 2);
    length(V_trials)
    M_trials = find(allData{1,t}.stimorder == 3);
    length(M_trials)
    for r = 1:size(allData{1,t}.boolean_response,1)
        allData{1,t}.CI_percentresp(r, 1) = nansum(allData{1,t}.boolean_response(r, MS_trials)) / length(MS_trials);
        allData{1,t}.CI_percentresp(r, 2) = nansum(allData{1,t}.boolean_response(r, V_trials)) / length(V_trials);
        allData{1,t}.CI_percentresp(r, 3) = nansum(allData{1,t}.boolean_response(r, M_trials)) / length(M_trials);
        unimax = nanmax([allData{1,t}.CI_percentresp(r, 2), allData{1,t}.CI_percentresp(r, 3)]);
        allData{1,t}.CI_percentresp(r, 4) = (allData{1,t}.CI_percentresp(r, 1) - unimax) / (allData{1,t}.CI_percentresp(r, 1) + unimax); % contrast index
        unimean = nanmean([allData{1,t}.CI_percentresp(r, 2), allData{1,t}.CI_percentresp(r, 3)]);
        allData{1,t}.CI_percentresp(r, 5) = allData{1,t}.CI_percentresp(r, 1) - unimean;
    end
    clear('MS_trials', 'V_trials', 'M_trials')
end

% add contrast index by proportion responses to allRespROIs
idx = 1
for t = 1:length(allData)
    for r = 1:length(allData{1,t}.resp_ROIs)
        allRespROIs(idx, 20) = allData{1,t}.CI_percentresp(allData{1,t}.resp_ROIs(r), 4);
        idx = idx+1;
    end
end

% change in success rate (b/c data is alerady normalized to 1)
idx = 1
for t = 1:length(allData)
    for r = 1:length(allData{1,t}.resp_ROIs)
        allRespROIs(idx, 21) = allData{1,t}.CI_percentresp(allData{1,t}.resp_ROIs(r), 5);
        idx = idx+1;
    end
end


%% Stat test 
vars = [20 21]
for v = 1:length(vars)
    [h(v) p(v)] = kstest2(allRespROIs(st49, vars(v)), allRespROIs(st46, vars(v)))
    peak_median(v,1) = nanmedian(allRespROIs(st49, vars(v)));
    peak_median(v,2) = nanmedian(allRespROIs(st46, vars(v)));
    peak_median(v,3) = nanstd(allRespROIs(st49, vars(v)));
    peak_median(v,4) = nanstd(allRespROIs(st46, vars(v)));
end


%% Plot change in success rate 

% ECDF
[f1, x1] = ecdf(allRespROIs(st46, 21));
[f2, x2] = ecdf(allRespROIs(st49, 21));
figure;
hold on
plot(x1, f1, 'g', 'LineWidth', 1) %st 46 = green
plot(x2, f2, 'm', 'LineWidth', 1) %st 49 = purple
%legend({'stage 46', 'stage 49'})
hold off
ylabel('ROIs')
xlabel('\Delta proportion responses')
%xlim([-inf 3])
set(gca,'FontSize',30)
saveas(gcf, 'ecdf of MSind change in response prop', 'png')
saveas(gcf, 'ecdf of MSInd change in response prop', 'epsc2')

% 46 hist
figure;
histogram(allRespROIs(st46, 21), 60, 'FaceColor', 'g', 'EdgeColor', 'g', 'BinLimits', [-0.4 0.5]) 
%ylim([0 130])
ylabel('ROIs')
xlabel('\Delta proportion responses')
%xlim([0 7])
set(gca,'FontSize',30)
saveas(gcf, 'hist of MSInd change in response prop 46', 'png')
saveas(gcf, 'hist of MSInd change in response prop 46', 'epsc2')

% 49 hist
figure;
histogram(allRespROIs(st49, 20), 60, 'FaceColor', 'm', 'EdgeColor', 'm', 'BinLimits', [-0.4 0.5]) 
ylabel('ROIs')
xlabel('\Delta proportion responses')
%ylim([0 130])
%xlim([0 7])
set(gca,'FontSize',30)
saveas(gcf, 'hist of MSInd change in response prop  49', 'png')
saveas(gcf, 'hist of MSInd change in response prop  49', 'epsc2')



%% Old stuff 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 46, MS, V, M
clear('MS', 'M', 'V')
MS = allRespROIs(st46, 7);
V = allRespROIs(st46, 8);
M = allRespROIs(st46, 9);
figure;
hold on
ecdf(MS) %blue
ecdf(V) %orange
ecdf(M) %yellow
xlim([-0.6 2])
h = get(gca, 'children')
set(h, 'LineWidth', 3)
ylabel('ROIs')
xlabel('Mean Peak \DeltaF/F_{0}')
set(gca,'FontSize',20)
saveas(gcf, 'st46 ecdf of MS, V, M mean peak', 'png')
saveas(gcf, 'st46 ecdf of MS, V, M mean peak', 'epsc2')

% 49, MS, V, M
clear('MS', 'M', 'V')
MS = allRespROIs(st49, 7);
V = allRespROIs(st49, 8);
M = allRespROIs(st49, 9);
figure;
hold on
ecdf(MS) %blue
ecdf(V) %orange
ecdf(M) %yellow
xlim([-0.6 2])
h = get(gca, 'children')
set(h, 'LineWidth', 3)
ylabel('ROIs')
xlabel('Mean Peak \DeltaF/F_{0}')
set(gca,'FontSize',20)
saveas(gcf, 'st49 ecdf of MS, V, M mean peak', 'png')
saveas(gcf, 'st49 ecdf of MS, V, M mean peak', 'epsc2')

%% 3C histograms of peaks by modality

% 49 MS
figure;
hist(allRespROIs(st49, 7),30)
ylabel('ROI count')
xlabel('Peak (\DeltaF/F_{0})')
set(gca,'FontSize',20)
xlim([-.5 2])
ylim([0 350])
saveas(gcf, 'hist peak 49 MS', 'png')
saveas(gcf, 'hist peak 49 MS', 'epsc2')
close;
% 49 V
figure;
hist(allRespROIs(st49, 8),30)
ylabel('ROI count')
xlabel('Peak (\DeltaF/F_{0})')
set(gca,'FontSize',20)
xlim([-.5 2])
ylim([0 350])
saveas(gcf, 'hist peak 49 V', 'png')
saveas(gcf, 'hist peak 49 V', 'epsc2')
close;
% 49 M
figure;
hist(allRespROIs(st49, 9),30)
ylabel('ROI count')
xlabel('Peak (\DeltaF/F_{0})')
set(gca,'FontSize',20)
xlim([-.5 2])
ylim([0 350])
saveas(gcf, 'hist peak 49 M', 'png')
saveas(gcf, 'hist peak 49 M', 'epsc2')
close;

% 46 MS
figure;
hist(allRespROIs(st46, 7),30)
ylabel('ROI count')
xlabel('Peak (\DeltaF/F_{0})')
set(gca,'FontSize',20)
xlim([-.5 2])
ylim([0 200])
saveas(gcf, 'hist peak 46 MS', 'png')
saveas(gcf, 'hist peak 46 MS', 'epsc2')
close;
% 46 V
figure;
hist(allRespROIs(st46, 8),30)
ylabel('ROI count')
xlabel('Peak (\DeltaF/F_{0})')
set(gca,'FontSize',20)
xlim([-.5 2])
ylim([0 200])
saveas(gcf, 'hist peak 46 V', 'png')
saveas(gcf, 'hist peak 46 V', 'epsc2')
close;
% 46 M
figure;
hist(allRespROIs(st46, 9),30)
ylabel('ROI count')
xlabel('Peak (\DeltaF/F_{0})')
set(gca,'FontSize',20)
xlim([-.5 2])
ylim([0 200])
saveas(gcf, 'hist peak 46 M', 'png')
saveas(gcf, 'hist peak 46 M', 'epsc2')
close;

%% 3D ECDF Onset times

% 46, MS, V, M
clear('MS', 'M', 'V')
MS = allRespROIs(st46, 18);
V = allRespROIs(st46, 19);
M = allRespROIs(st46, 20);
figure;
hold on
ecdf(MS) %blue
ecdf(V) %orange
ecdf(M) %yellow
xlim([0 7])
h = get(gca, 'children')
set(h, 'LineWidth', 3)
ylabel('ROIs')
xlabel('Mean Peak \DeltaF/F_{0}')
set(gca,'FontSize',20)
saveas(gcf, 'st46 ecdf of MS, V, M mean onset time', 'png')
saveas(gcf, 'st46 ecdf of MS, V, M mean onset time', 'epsc2')

% 49, MS, V, M
clear('MS', 'M', 'V')
MS = allRespROIs(st49, 18);
V = allRespROIs(st49, 19);
M = allRespROIs(st49, 20);
figure;
hold on
ecdf(MS) %blue
ecdf(V) %orange
ecdf(M) %yellow
xlim([0 7])
h = get(gca, 'children')
set(h, 'LineWidth', 3)
ylabel('ROIs')
xlabel('Mean Peak \DeltaF/F_{0}')
set(gca,'FontSize',20)
saveas(gcf, 'st49 ecdf of MS, V, M mean onset time', 'png')
saveas(gcf, 'st49 ecdf of MS, V, M mean onset time', 'epsc2')

%% 3D histograms of onset times

% 49 MS
figure;
hist(allRespROIs(st49, 18),30)
ylabel('Onset time (sec)')
xlabel('ROI count')
set(gca,'FontSize',20)
xlim([0 7])
ylim([0 120])
saveas(gcf, 'hist onsettime 49 MS', 'png')
saveas(gcf, 'hist onsettime 49 MS', 'epsc2')
close;
% 49 V
figure;
hist(allRespROIs(st49, 19),30)
ylabel('Onset time (sec)')
xlabel('ROI count')
set(gca,'FontSize',20)
xlim([0 7])
ylim([0 120])
saveas(gcf, 'hist onsettime 49 V', 'png')
saveas(gcf, 'hist onsettime 49 V', 'epsc2')
close;
% 49 M
figure;
hist(allRespROIs(st49, 20),30)
ylabel('Onset time (sec)')
xlabel('ROI count')
set(gca,'FontSize',20)
xlim([0 7])
ylim([0 120])
saveas(gcf, 'hist onsettime 49 M', 'png')
saveas(gcf, 'hist onsettime 49 M', 'epsc2')
close;
% 46 MS
figure;
hist(allRespROIs(st46, 18),30)
ylabel('Onset time (sec)')
xlabel('ROI count')
set(gca,'FontSize',20)
xlim([0 7])
ylim([0 120])
saveas(gcf, 'hist onsettime 46 MS', 'png')
saveas(gcf, 'hist onsettime 46 MS', 'epsc2')
close;
% 46 V
figure;
hist(allRespROIs(st46, 19),30)
ylabel('Onset time (sec)')
xlabel('ROI count')
set(gca,'FontSize',20)
xlim([0 7])
ylim([0 120])
saveas(gcf, 'hist onsettime 46 V', 'png')
saveas(gcf, 'hist onsettime 46 V', 'epsc2')
close;
% 46 M
figure;
hist(allRespROIs(st46, 20),30)
ylabel('Onset time (sec)')
xlabel('ROI count')
set(gca,'FontSize',20)
xlim([0 7])
ylim([0 120])
saveas(gcf, 'hist onsettime 46 M', 'png')
saveas(gcf, 'hist onsettime 46 M', 'epsc2')
close;

%% E. Multisensory Index by peak ECDF

[f1, x1] = ecdf(allRespROIs(st49, 11)); % MSInd_peak
%set(f1, 'LineColor', 'm')
[f2, x2] = ecdf(allRespROIs(st46, 11)) %, 'g') % MSInd_peak
figure;
hold on
plot(x1, f1, 'm', 'LineWidth', 3)
plot(x2, f2, 'g', 'LineWidth', 3)
hold off
ylabel('ROIs')
xlabel('MSIndex (peak \DeltaF/F_{0})')
xlim([-6 15])
set(gca,'FontSize',20)
saveas(gcf, 'ecdf MSIndex peak', 'png')
saveas(gcf, 'ecdf MSIndex peak', 'epsc2')


%% E.  Multisensory Index by peak histogram

figure;
hist(allRespROIs(st49, 11),100);
xlabel('MSIndex')
ylabel('ROI count')
xlim([-6 15])
ylim([0 300])
set(gca,'FontSize',20)
saveas(gcf, 'hist MSIndex peak 49', 'png')
saveas(gcf, 'hist MSIndex peak 49', 'epsc2')
figure;
hist(allRespROIs(st46, 11),40);
xlabel('MSIndex')
ylabel('ROI count')
xlim([-6 15])
ylim([0 300])
set(gca,'FontSize',20)
saveas(gcf, 'hist MSIndex peak 46', 'png')
saveas(gcf, 'hist MSIndex peak 46', 'epsc2')

%% E. Multisensory Index by onset time ECDF

[f1, x1] = ecdf(allRespROIs(st49, 13)); % MSInd_peak
%set(f1, 'LineColor', 'm')
[f2, x2] = ecdf(allRespROIs(st46, 13)) %, 'g') % MSInd_peak
figure;
hold on
plot(x1, f1, 'm', 'LineWidth', 3)
plot(x2, f2, 'g', 'LineWidth', 3)
hold off
ylabel('ROIs')
xlabel('MSIndex (peak \DeltaF/F_{0})')
xlim([-2 15])
set(gca,'FontSize',20)
saveas(gcf, 'ecdf MSIndex onset time', 'png')
saveas(gcf, 'ecdf MSIndex onset time', 'epsc2')


%% E. Multisensory Index by onset time histogram

figure;
hist(allRespROIs(st49, 13),40);
ylabel('MSIndex')
xlabel('ROI count')
set(gca,'FontSize',20)
xlim([-1 15])
ylim([0 180])
saveas(gcf, 'hist MSIndex onset time 49', 'png')
saveas(gcf, 'hist MSIndex onset time 49', 'epsc2')
figure;
hist(allRespROIs(st46, 13),100);
ylabel('MSIndex')
xlabel('ROI count')
set(gca,'FontSize',20)
xlim([-1 15])
ylim([0 180])
saveas(gcf, 'hist MSIndex onset time 46', 'png')
saveas(gcf, 'hist MSIndex onset time 46', 'epsc2')


%% NOtes

    s46 = allRespROIs(st46, i);
    s49 = allRespROIs(st49, i);
    figure;
    hold on
    ecdf(s46)

    %set(h, 'Color', 'g')
    ecdf(s49)
    h = get(gca, 'children')
    set(h, 'LineWidth', 3)
    %j = get(gca, 'children')
    %set(j, 'LineWidth', 3)
    %set(j, 'Color', [0.5 0 0.5])
    hold off
    title(sprintf('%s ECDF by stage', labels{i-2}))
    xlabel(labels{i-2})
    ylabel('ROI count')
    fig_filename = sprintf('st 46 vs 49 ECDF of %s', labels{i-2})
    saveas(gcf, fig_filename, 'png')
    close;
%% Dissertation figure 3: 46 vs 49 diffs
% using corrected for bad trials data 
% use xcorr_cleaned_20180313 file 
%load('D:\Torrey_calcium_imaging\compare_46-49\analysis_Feb 2018\corrected_for_badtrials\46_49_comparison\xcorr_corrected_20180313.mat')

% steal figures from compare_46_49_new, but make them pretty


%% 3B: prop resp in respROI

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

%% 3C: peaks EDCF

st46 = find(allRespROIs(:,28) == 46);
st49 = find(allRespROIs(:,28) == 49);

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
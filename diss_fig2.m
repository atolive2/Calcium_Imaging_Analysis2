%% characterise all tads
% dissertation figure 2 

% use xcorr_cleaned_20180313 file 
%load('D:\Torrey_calcium_imaging\compare_46-49\analysis_Feb 2018\corrected_for_badtrials\46_49_comparison\xcorr_corrected_20180313.mat')


%% 2B - prop ROIs respond by tad

% scatterplot the values by stage - proportion of respROIs
% taken from compare_46_49_new
figure;
hold on
scatter(st46_X, ct_respROIs(s46_tads,3), 40, 'g', 'filled')
scatter(st49_X, ct_respROIs(s49_tads,3), 40, [0.5 0 0.5], 'filled')
hold off
xlim([0.5 2.5])
ax = gca;
ax.XTick = [1 2];
ax.XTickLabel = [46 49];
%title('Proportion of Responding ROIs by Exp')
xlabel('Stage')
ylabel('Proportion of ROIs that respond')
set(gca,'FontSize',20)
saveas(gcf, 'prop resROIs by stage', 'png')
saveas(gcf, 'prop resROIs by stage', 'epsc2')

%% 2C avearge peaks by modality
% 46 and 49 combined

% multisensory
figure;
hist(allRespROIs(:, 7), 40)
xlim([-0.5 2])
ylim([0 500])
xlabel('Mean Peak \DeltaF/F_{0}')
ylabel('ROI count')
set(gca,'FontSize',20)
fig_filename = 'peaks by mod respROIs MS'
saveas(gcf, fig_filename, 'png')
saveas(gcf, fig_filename, 'epsc2')
%visual
figure;
hist(allRespROIs(:, 8), 40)
xlim([-0.5 2])
ylim([0 500])
xlabel('Mean Peak \DeltaF/F_{0}')
ylabel('ROI count')
set(gca,'FontSize',20)
fig_filename = 'peaks by mod respROIs V'
saveas(gcf, fig_filename, 'png')
saveas(gcf, fig_filename, 'epsc2')
%mechanosensory
figure;
hist(allRespROIs(:, 9), 40)
xlim([-0.5 2])
ylim([0 500])
xlabel('Mean Peak \DeltaF/F_{0}')
ylabel('ROI count')
set(gca,'FontSize',20)
fig_filename = 'peaks by mod respROIs M'
saveas(gcf, fig_filename, 'png')
saveas(gcf, fig_filename, 'epsc2')


%% 2D proportion of responses respROIs

% Get the number of responses / total number of good trial presentations
for t = 1:length(allData)
    for i = 1:size(allData{1,t}.include, 1) %over all stimtypes
        stim_resp = allData{1,t}.boolean_response & squeeze(allData{1,t}.include(i,:,:));
        stim_resp_ct = sum(stim_resp,2);
        stim_resp_prop(i,:) = stim_resp_ct ./ sum(squeeze(allData{1,t}.include(i,:,:)),2);
    end
    allData{1,t}.resp_prop_bystim = stim_resp_prop;
    clear('stim_resp', 'stim_resp_ct', 'stim_resp_prop')
end

% Put the prop responses values for respROIs into a matrix
% only use 1-3 (high/high)
prop_resp_bystim = [];
for t = 1:length(allData)
    prop_resp_bystim = [prop_resp_bystim, allData{1,t}.resp_prop_bystim(1:3, allData{1,t}.resp_ROIs)];
end

%make histograms from prop_resp_bystim

%multisensory
figure;
hist(prop_resp_bystim(1, :), 40)
xlim([0 1])
ylim([0 500])
xlabel('Proportion of trials with response')
ylabel('ROI count')
set(gca,'FontSize',20)
fig_filename = 'prop responses respROIs MS'
saveas(gcf, fig_filename, 'png')
saveas(gcf, fig_filename, 'epsc2')

%visual
figure;
hist(prop_resp_bystim(2, :), 40)
xlim([0 1])
ylim([0 500])
xlabel('Proportion of trials with response')
ylabel('ROI count')
set(gca,'FontSize',20)
fig_filename = 'prop responses respROIs V'
saveas(gcf, fig_filename, 'png')
saveas(gcf, fig_filename, 'epsc2')

%mech
figure;
hist(prop_resp_bystim(3, :), 40)
xlim([0 1])
ylim([0 500])
xlabel('Proportion of trials with response')
ylabel('ROI count')
set(gca,'FontSize',20)
fig_filename = 'prop responses respROIs M'
saveas(gcf, fig_filename, 'png')
saveas(gcf, fig_filename, 'epsc2')













%% Notes

set(findall(gcf,'-property','FontSize'),'FontSize',30)


% apply to all
xlim([-0.5 2])
ylim([0 500])
xlabel('Mean Peak \DeltaF/F_{0}')
ylabel('ROI count')
set(gca,'FontSize',20)
fig_filename = 'peaks by mod respROIs MS'
saveas(gcf, fig_filename, 'png')
saveas(gcf, fig_filename, 'epsc2')

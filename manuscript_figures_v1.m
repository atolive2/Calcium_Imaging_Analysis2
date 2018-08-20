%% Ca img manuscript figures 
% This is the updated figures based on changes necessary from my
% dissertation. 

%% Began with load('D:\Torrey_calcium_imaging\compare_46-49\analysis_Feb 2018\corrected_for_badtrials\46_49_comparison\AllData_20180320.mat')

% HOWEVER, all generated plots will be saved in a separate folder, 
folder = 'D:\Torrey_calcium_imaging\manuscript_figures\Version_1'
% use to save files to the right place: saveas(gcf, fullfile(folder, filename), 'png')

%% Prepare workspace:

% get list of 46 and 49 tads
st49_tads = [];
st46_tads = [];
for t = 1:length(allData)
    if allData{1,t}.stage == 49
    st49_tads = [st49_tads, t];
    elseif allData{1,t}.stage == 46
    st46_tads = [st46_tads, t];
    end
end

st_ID = [46*ones(1,length(st46_tads)), 49*ones(1,length(st49_tads))]
%% Figure 1 part 1: St 46 tads have more generally responsive cells 

% Proportion of ROIs that respond by tadpole
for t = 1:length(st49_tads)
    st49_proprespROIs(t) = length(allData{1, st49_tads(t)}.resp_ROIs) / size(allData{1, st49_tads(t)}.ROIcenters, 1);
end

for t = 1:length(st46_tads)
    st46_proprespROIs(t) = length(allData{1, st46_tads(t)}.resp_ROIs) / size(allData{1, st46_tads(t)}.ROIcenters, 1);
end

%make plot
% first put all values in 1 vector and create a grouping variable
st_ID = [46*ones(1,length(st46_tads)), 49*ones(1,length(st49_tads))]
all_proprespROIs = [st46_proprespROIs, st49_proprespROIs]

figure;
hold on
boxplot(all_proprespROIs, st_ID, 'Colors', 'gm')
plot(ones(length(st46_tads)), st46_proprespROIs, 'go')
plot(2*ones(length(st49_tads)), st49_proprespROIs, 'mo')
hold off 
xlabel('Stage')
ylabel('Proportion Responding ROIs')
set(gca,'FontSize',20)
set(findobj(gca,'type','line'),'linew',2)
fig_filename = 'f1p1_proprespROIsbytad'
saveas(gcf, fullfile(folder, fig_filename), 'png')
saveas(gcf, fullfile(folder, fig_filename), 'epsc2')


% % proportion responses for ROIs by tadpole
for t = 1:length(allData)
    tmp = allData{1,t}.sum_responses(allData{1,t}.resp_ROIs) / length(allData{1,t}.stimorder)
    prop_resp_trials(t) = mean(tmp)
end

%make plot
sorted_propresptrials = [prop_resp_trials(st46_tads), prop_resp_trials(st49_tads)];
figure;
hold on
boxplot(sorted_propresptrials, st_ID, 'Colors', 'gm')
plot(ones(length(st46_tads)), prop_resp_trials(st46_tads), 'go')
plot(2*ones(length(st49_tads)), prop_resp_trials(st49_tads), 'mo')
hold off 
xlabel('Stage')
ylabel('Proportion trials with response')
set(gca,'FontSize',20)
set(findobj(gca,'type','line'),'linew',2)
fig_filename = 'f1p1_propresptrialssbyROIbytad'
saveas(gcf, fullfile(folder, fig_filename), 'png')
saveas(gcf, fullfile(folder, fig_filename), 'epsc2')

% Response size for respROIs by tadpole
% % peak size for trials with resp for each resp_roi
% %     if boolean response, include peakbytrial2 in avg
% % avg peaks = per roi
% % avg rois = per tad

for t = 1:length(allData)
    avg_peaksize(t) = mean(mean(allData{1,t}.avg_peak(:,allData{1,t}.resp_ROIs)))
end

%make plot
sorted_avg_peaksize = [avg_peaksize(st46_tads), avg_peaksize(st49_tads)];
figure;
hold on
boxplot(sorted_avg_peaksize, st_ID, 'Colors', 'gm')
plot(ones(length(st46_tads)), avg_peaksize(st46_tads), 'go')
plot(2*ones(length(st49_tads)), avg_peaksize(st49_tads), 'mo')
hold off 
xlabel('Stage')
ylabel('Peak Size')
set(gca,'FontSize',20)
set(findobj(gca,'type','line'),'linew',2)
fig_filename = 'f1p1_avg_peaksize_bytad'
saveas(gcf, fullfile(folder, fig_filename), 'png')
saveas(gcf, fullfile(folder, fig_filename), 'epsc2')

%% Figure 1 part 2: 46 respond more quickly

% Peak response size distribution (avg each tad's distribution)


% Latency of response distribution (avg each tad's distribution)





%% Figure 2 Multisensory integration, norm by tad
% Standard MSI (peak response)
for t = 1:length(allData)
    avg_MSInd_peak(t) = mean(mean(allData{1,t}.MSInd_peak(:,allData{1,t}.resp_ROIs)));
end

%make plot
sorted_avg_MSInd_peak = [avg_MSInd_peak(st46_tads), avg_MSInd_peak(st49_tads)];
figure;
hold on
boxplot(sorted_avg_MSInd_peak, st_ID, 'Colors', 'gm')
plot(ones(length(st46_tads)), avg_MSInd_peak(st46_tads), 'go')
plot(2*ones(length(st49_tads)), avg_MSInd_peak(st49_tads), 'mo')
hold off 
xlabel('Stage')
ylabel('MSInd')
set(gca,'FontSize',20)
set(findobj(gca,'type','line'),'linew',2)
fig_filename = 'f2p1_avg_MSInd_peak_bytad'
saveas(gcf, fullfile(folder, fig_filename), 'png')
saveas(gcf, fullfile(folder, fig_filename), 'epsc2')

% MS Timing
for t = 1:length(allData)
    avg_MSInd_onsettime(t) = mean(mean(allData{1,t}.MSInd_onsettime(:,allData{1,t}.resp_ROIs)));
end

%make plot
sorted_avg_MSInd_onsettime = [avg_MSInd_onsettime(st46_tads), avg_MSInd_onsettime(st49_tads)];
figure;
hold on
boxplot(sorted_avg_MSInd_onsettime, st_ID, 'Colors', 'gm')
plot(ones(length(st46_tads)), avg_MSInd_onsettime(st46_tads), 'go')
plot(2*ones(length(st49_tads)), avg_MSInd_onsettime(st49_tads), 'mo')
hold off 
xlabel('Stage')
ylabel('MSInd')
set(gca,'FontSize',20)
set(findobj(gca,'type','line'),'linew',2)
fig_filename = 'f2p1_avg_MSInd_onsettime_bytad'
saveas(gcf, fullfile(folder, fig_filename), 'png')
saveas(gcf, fullfile(folder, fig_filename), 'epsc2')


% Change in response reliability



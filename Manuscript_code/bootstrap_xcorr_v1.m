%% Shuffle and bootstrap xcorr values

% start with allData containing respROIs_ST, smoothed_trunc, and the list
% of trials in order. 
% folder: D:\Torrey_calcium_imaging\manuscript_figures\Version_2
load('allData_workspace_20181015.mat')

for k = 1:length(allData)
    bootData(k).smoothed_trunc = allData(k).smoothed_trunc;
end
for k = 1:length(allData)
    bootData(k).respROIs_ST = allData(k).respROIs_ST;
    bootData(k).incl_trials = allData(k).incl_trials;
end
%% Shuffle the trials, calculate xcorr and report maxR, lagmaxR and 0lagR into corresponding fields 

its = 100; % how many bootstrap iterations?
stims = [1 2 3 4] % MS, V, M, NS
% allData(k).incl_trials contains the trial numbers for each ROI pair for
% each stim type
tic
for k = 1:length(allData)
    if length(bootData(k).respROIs_ST) > 5
        set_lag = length(bootData(k).smoothed_trunc{1,1});
        for s = 1:length(stims) %all 4 stims
            for r1 = 1:length(bootData(k).respROIs_ST) %all responding ROIs
                for r2 = 1:length(bootData(k).respROIs_ST) %all responding ROIs 
                    all_tr = bootData(k).incl_trials{s, r1, r2};   
                    if ~isempty(all_tr)
                        for b = 1:its % number of iterations
                            r1_tr = all_tr(randperm(length(all_tr)));
                            r2_tr = all_tr(randperm(length(all_tr)));
                            tmp1 = [];
                            for tr = 1:length(r1_tr)
                                tmp1 = [tmp1; bootData(k).smoothed_trunc{bootData(k).respROIs_ST(r1), r1_tr(tr)}];
                            end
                            tmp2 = [];
                            for t2 = 1:length(r2_tr)
                                tmp2 = [tmp2; bootData(k).smoothed_trunc{bootData(k).respROIs_ST(r2), r2_tr(t2)}];
                            end
                            [boot_xcorr(k).R{s, r1, r2, b}, boot_xcorr(k).lag{s, r1, r2, b}] = xcorr(tmp1, tmp2, set_lag, 'coeff');
                            [boot_xcorr(k).maxR{s, r1, r2}(b), boot_xcorr(k).lagmaxR{s, r1, r2}(b)] = max(boot_xcorr(k).R{s, r1, r2, b}(:,:));
                            boot_xcorr(k).lag0R{s, r1, r2}(b) = boot_xcorr(k).R{s, r1, r2, b}(set_lag+1);
                        end
                        clear('r1_tr', 'r2_tr')
                    end 
                    clear('all_tr')
                end
            end
        end
        toc
    end
end

                    
%% Generate 95% confidence intervals from the bootstrap values for each ROI pair
tic
for k = 1:length(allData)
    if length(allData(k).respROIs_ST) > 5
        for s = 1:length(stims) %all 4 stims
            for r1 = 1:length(allData(t).respROIs_ST) %all responding ROIs
                for r2 = 1:length(allData(t).respROIs_ST) %all responding ROIs 
                    % maxR calculations
                    boot_xcorr(k).maxR_info{s,r1,r2}(1) = mean(boot_xcorr(k).maxR{s,r1,r2}(:)); %mean
                    boot_xcorr(k).maxR_info{s,r1,r2}(2) = std(boot_xcorr(k).maxR{s,r1,r2}(:)); %std
                    boot_xcorr(k).maxR_info{s,r1,r2}(3) = boot_xcorr(k).maxR_info{s,r1,r2}(1) + (2*boot_xcorr(k).maxR_info{s,r1,r2}(2)); %upper 95% CI
                    boot_xcorr(k).maxR_info{s,r1,r2}(4) = boot_xcorr(k).maxR_info{s,r1,r2}(1) - (2*boot_xcorr(k).maxR_info{s,r1,r2}(2)); %lower 95% CI
                    % lag of maxR calculations
                    boot_xcorr(k).lagmaxR_info{s,r1,r2}(1) = mean(boot_xcorr(k).lagmaxR{s,r1,r2}(:)); %mean
                    boot_xcorr(k).lagmaxR_info{s,r1,r2}(2) = std(boot_xcorr(k).lagmaxR{s,r1,r2}(:)); %std
                    boot_xcorr(k).lagmaxR_info{s,r1,r2}(3) = boot_xcorr(k).lagmaxR_info{s,r1,r2}(1) + (2*boot_xcorr(k).lagmaxR_info{s,r1,r2}(2)); %upper 95% CI
                    boot_xcorr(k).lagmaxR_info{s,r1,r2}(4) = boot_xcorr(k).lagmaxR_info{s,r1,r2}(1) - (2*boot_xcorr(k).lagmaxR_info{s,r1,r2}(2)); %lower 95% CI
                    % lag 0 R calculations
                    boot_xcorr(k).lag0R_info{s,r1,r2}(1) = mean(boot_xcorr(k).lag0R{s,r1,r2}(:)); %mean
                    boot_xcorr(k).lag0R_info{s,r1,r2}(2) = std(boot_xcorr(k).lag0R{s,r1,r2}(:)); %std
                    boot_xcorr(k).lag0R_info{s,r1,r2}(3) = boot_xcorr(k).lag0R_info{s,r1,r2}(1) + (2*boot_xcorr(k).lag0R_info{s,r1,r2}(2)); %upper 95% CI
                    boot_xcorr(k).lag0R_info{s,r1,r2}(4) = boot_xcorr(k).lag0R_info{s,r1,r2}(1) - (2*boot_xcorr(k).lag0R_info{s,r1,r2}(2)); %lower 95% CI
                end
            end
        end
        toc
    end
end
%%%%%%%%%%%%%% ABOVE THIS POINT: ADD CHANGES FROM CCV

%% Compare actual data point to 95% confidence interval for each ROI pair 
% have allData and bootData_CI open

for k = 1:length(boot_xcorr)
    if length(allData(k).respROIs_ST) > 10
        for s = 1:length(stims) %all 4 stims
            for r1 = 1:length(allData(k).respROIs_ST) %all responding ROIs
                for r2 = 1:length(allData(k).respROIs_ST) %all responding ROIs 
                    
                    if allData(k).respROIdff0_maxR_sq(s, r1, r2) < boot_xcorr(k).maxR_info(s, r1, r2, 4) % A<L smaller than L CI => outside CI
                        boot_xcorr(k).comparetoCI(s, r1, r2, 1) = 1; % actual data is outside CI
                        boot_xcorr(k).comparetoCI(s, r1, r2, 2) = boot_xcorr(k).maxR_info(s, r1, r2, 4) - allData(k).respROIdff0_maxR_sq(s, r1, r2); % LCI - A to get +diff
                        boot_xcorr(k).comparetoCI(s, r1, r2, 3) = -1; %ID those smaller than lower CI bound with -1
                    else %larger than L CI
                        if allData(k).respROIdff0_maxR_sq(s, r1, r2) < boot_xcorr(k).maxR_info(s, r1, r2, 3) % A<U %inside CI
                            boot_xcorr(k).comparetoCI(s, r1, r2, 1) = 0; % actual data is inside CI
                            tmp(1) = boot_xcorr(k).maxR_info(s, r1, r2, 3) - allData(k).respROIdff0_maxR_sq(s, r1, r2); %U - A (positive val)
                            tmp(2) = allData(k).respROIdff0_maxR_sq(s, r1, r2) - boot_xcorr(k).maxR_info(s, r1, r2, 4); %A - L (positive val)
                            boot_xcorr(k).comparetoCI(s, r1, r2, 2) = -min(tmp);
                            boot_xcorr(k).comparetoCI(s, r1, r2, 3) = 0; %ID those inside CI with 0
                        else %A>U larger than U CI => outside CI
                            boot_xcorr(k).comparetoCI(s, r1, r2, 1) = 1; % actual data is outside CI
                            boot_xcorr(k).comparetoCI(s, r1, r2, 2) = allData(k).respROIdff0_maxR_sq(s, r1, r2) - boot_xcorr(k).maxR_info(s, r1, r2, 3); %A-U
                            boot_xcorr(k).comparetoCI(s, r1, r2,3) = 1; % ID those larger than CI with +1
                        end
                    end
                end
            end
        end
    end
end
toc

%% Plot difference between edge of 95% confidence interval and the actual datapoint for each tad
% subplot each stim type
% histogram
% postive = outside (significantly diff from noise)
% negative = inside distribution (not sig diff from noise)
k = 22
% for k = 1:length(allData)
%     if length(allData(k).respROIs_ST) > 5
         for s = 1 %MS only now %:length(stims) %all 4 stims
             for r1 = 1:size(boot_xcorr(k).maxR_info, 2) 
                 for r2 = 1:size(boot_xcorr(k).maxR_info, 3) 
                     figure;
                     hold on
                     histogram(boot_xcorr(k).maxR(s, r1, r2, :),'Normalization', 'probability')
                     plot([boot_xcorr(k).maxR_info(s, r1, r2, 4), boot_xcorr(k).maxR_info(s, r1, r2, 3)], [0.1 0.1], '*k')
                     plot([allData(k).respROIdff0_maxR_sq(s, r1, r2)], 0.1, 'dm')
                     hold off
                     xlabel('maxR value')
                     ylabel('counts')
                     title(sprintf('Tad %d (k=%d), s=%d r1=%d, r2=%d, bootstrap CI hist with actual', allData(k).expnum, k, s, r1, r2))
                     fig_filename = sprintf('Tad %d (k=%d), s=%d r1=%d, r2=%d,bootstrap CI hist with actual', allData(k).expnum, k, s, r1, r2)
                     saveas(gcf, fig_filename, 'png')
                     saveas(gcf, fig_filename, 'epsc2')
                     close;
                 end
             end
         end
         

%% Plot proportion of pairs outside 95% confidence interval

for k = 1:length(boot_xcorr)
    if ~isempty(boot_xcorr(k).comparetoCI)
    for s = length(stims)
        bar_data(k,1) = length(find(boot_xcorr(k).comparetoCI(s, :, :, 3) == 1));
        bar_data(k, 2) = length(find(boot_xcorr(k).comparetoCI(s, :, :, 3) == 0));
        bar_data(k, 3) = length(find(boot_xcorr(k).comparetoCI(s, :, :, 3) == -1));
    end
    end
end
bar_sum = sum(bar_data,2)
for m = 1:size(bar_data, 2)
    prop_data(:,m) = bar_data(:, m) ./ bar_sum
end
sum(prop_data, 2)

% sort by stage
st46_tads = find(stage == 46)
st49_tads = find(stage == 49)

% eliminate bad tads
for k = 1:length(allData)
    keepers(k) = (length(allData(k).respROIs_ST) > 10);
end
list = find(keepers)
st46_keep = intersect(st46_tads, list)
st49_keep = intersect(st49_tads, list)

plot_data = [prop_data(st46_keep, :); prop_data(st49_keep, :)]
tad_colors = [repmat([0 1 0], [length(st46_keep), 1]); repmat([1 0 1], [length(st49_keep), 1])];
figure;
bar(plot_data, 'stacked', 'FaceColor', 'flat')
annotation('line', [.498 .498], [0 .805], 'LineWidth', 2)
xlabel('stage 46                    stage 49')
ylabel('Proportion ROI pairs')

legend('Larger', 'Inside', 'Smaller', 'Location', 'northoutside', 'Orientation', 'horizontal')
title('All ROI pairs actual maxR vs BS CI')
set(gca, 'FontSize', 25)
fig_filename = 'All ROI pairs actual maxR vs BS CI'
saveas(gcf, fig_filename, 'png')
saveas(gcf, fig_filename, 'epsc2')
close;

%% Plot histogram of values for difference between CI and actual
stim_ID = {'MS', 'V', 'M', 'N'};
stim_color = [255/256,215/256,0; 1 0 0; 0 0 1; 0.25 0.25 0.25];
for k = 1:length(boot_xcorr)
    if ~isempty(boot_xcorr(k).comparetoCI)
        figure;
        for s = 1:length(stims)
            subplot(2, 2, s)
            histogram(boot_xcorr(k).comparetoCI(s, :, :, 2), 50, 'FaceColor', stim_color(s, :), 'normalization', 'probability')
            axis tight
            title(sprintf('%s', stim_ID{s}))
        end
        suptitle(sprintf('Tad %d (k=%d) All ROI pairs diff from CI maxR', allData(k).expnum, k))
        fig_filename = sprintf('Tad %d (k=%d) All ROI pairs diff from CI maxR', allData(k).expnum, k);
        saveas(gcf, fig_filename, 'png')
    saveas(gcf, fig_filename, 'epsc2')
    close;
    end
end

%% Calculate how many ROI pairs are outside the CI for each ROI
% use this information to develop HighCorr definition

for k = 1:length(boot_xcorr)
    if ~isempty(boot_xcorr(k).comparetoCI)
        for s = 1:length(stims)
            for r = 1:size(boot_xcorr(k).comparetoCI, 2)
                boot_xcorr(k).ROI_pair_info(s,r) = length(find(boot_xcorr(k).comparetoCI(s,r,:,1)));
            end
        end
    end
end

% histogram each tad
for k = 1:length(boot_xcorr)
    if ~isempty(boot_xcorr(k).comparetoCI)
        figure;
        for s = 1:length(stims)
            subplot(2, 2, s)
            histogram((boot_xcorr(k).ROI_pair_info(s,:) ./ size(boot_xcorr(k).ROI_pair_info, 2)), 'FaceColor', stim_color(s, :), 'normalization', 'probability')
            axis tight
            title(sprintf('%s', stim_ID{s}))
        end
        suptitle(sprintf('Tad %d (k=%d) All ROIs how many pairs outside CI maxR', allData(k).expnum, k))
        fig_filename = sprintf('Tad %d (k=%d) All ROI how many pairs outside CI maxR', allData(k).expnum, k);
        saveas(gcf, fig_filename, 'png')
    saveas(gcf, fig_filename, 'epsc2')
    close;
    end
end

%% How many ROI pairs are outside the CI for each ROI grouped plots

% make histogram objects that are all the same range and number of bins
% min = -0.5, max = 1. 10 bins per 0.1 = 150 bins
num_bins = 150;
for k = 1:length(boot_xcorr)
    if ~isempty(boot_xcorr(k).comparetoCI)
        for s = 1:length(stims)
            ROI_CIdata(k, s, :) = histcounts(boot_xcorr(k).comparetoCI(s, :, :, 2), num_bins, 'BinLimits', [-0.5 1],   'normalization', 'probability');
            %ROI_CIdata(k, s, :) = h.Values
        end
    end
end

% get median and CI for each stage's histogram
% nonparametric CI can be approximated using this algorithm: https://onlinecourses.science.psu.edu/stat414/node/316/
% for N=9, the 96% CI is between the second and 9th values. For N=10, the
% 97% CI is between the second and 10th values 
for s = 1:length(stims)
    CI_median(s,:,1) = median(squeeze(ROI_CIdata(st46_keep, s, :))); %46 median
    CI_median(s,:,2) = median(squeeze(ROI_CIdata(st49_keep, s, :))); %49 median
    tmp1 = sort(squeeze(ROI_CIdata(st46_keep, s, :))); %sort st 46 vals
    CI_median(s, :, 3) = CI_median(s,:, 1) - tmp1(2, :); % st46 lower CI
    CI_median(s, :, 4) = tmp1(9, :) - CI_median(s,:, 1); %st46 upper CI
    tmp2 = sort(squeeze(ROI_CIdata(st49_keep, s, :))); %sort st 49 vals
    CI_median(s, :, 5) = CI_median(s, :, 2) - tmp2(2, :); % st49 lower CI
    CI_median(s, :, 6) = tmp2(10, :) - CI_median(s, :, 2); %st49 upper CI
end

% create figure 
X_val = -0.49:0.01:1;
figure;
for s = 1:length(stims)
    subplot(2, 2, s)
    hold on
%     for t = 1:length(st46_keep)
%         plot(X_val, squeeze(ROI_CIdata(st46_keep(t), s, :)), 'g')
%     end
%     for t = 1:length(st49_keep)
%         plot(X_val, squeeze(ROI_CIdata(st49_keep(t), s, :)), 'm')
%     end
    boundedline(X_val, CI_median(s,:,1), squeeze(CI_median(s,:,3:4)), 'cmap', [0 (100/256) 0], 'alpha', 'transparency', 0.3) %dark green is more visible
    boundedline(X_val, CI_median(s,:,2), squeeze(CI_median(s,:,5:6)), 'm', 'alpha', 'transparency', 0.1)
    
    hold off
    axis tight
    title(sprintf('%s', stim_ID{s}))
end

    suptitle('All Tads All ROI pairs diff from CI maxR')
    fig_filename = 'All Tads All ROI pairs diff from CI maxR';
    saveas(gcf, fig_filename, 'png')
    saveas(gcf, fig_filename, 'epsc2')
    %close;

%% How many ROIs have at least 20% of ROI pairs with actual maxR outside of CI?
threshold = 0.3;
for k = 1:length(boot_xcorr)
    if ~isempty(boot_xcorr(k).comparetoCI)
        numROIs = size(boot_xcorr(k).comparetoCI, 2);
        for s = 1:length(stims)
            for r = 1:numROIs
                boot_xcorr(k).highcorr(s,r) = length(find(boot_xcorr(k).comparetoCI(s, r, :, 1)));
            end
        end
        boot_xcorr(k).highcorrPCT = boot_xcorr(k).highcorr ./ numROIs;
        for s = 1:length(stims)
            boot_xcorr(k).highcorrROIs{s} = find(boot_xcorr(k).highcorrPCT(s, :) > threshold);
        end
    end
end

% calculate the number of highcorr ROIs for each stage and stimtype
for k = 1:length(boot_xcorr)
    for s = 1:length(boot_xcorr(k).highcorrROIs)
    numHCROIs(k, s) = length(boot_xcorr(k).highcorrROIs{s});
    end
end

%proportion highcorr ROIs
for k = 1:length(allData)
    numrespROIs(k) = length(allData(k).respROIs_ST);
end

for s = 1:size(numHCROIs,2)
    propHCROIs(:, s) = numHCROIs(:,s) ./ numrespROIs(1:28)'
end

% plot the proportion of ROIs that are HC for each stage and stim type
% get data for boxplot into one vector of labels and 1 vector of data
box_ids = [ones(1, 9), 2*ones(1, 9), 3*ones(1, 9), 4*ones(1, 10), 5*ones(1, 10), 6*ones(1, 10)]
box_data = [propHCROIs(st46_keep, 1)', propHCROIs(st46_keep, 2)', propHCROIs(st46_keep, 3)', propHCROIs(st49_keep, 1)', propHCROIs(st49_keep, 2)', propHCROIs(st49_keep, 3)']
%box_colors = [0 0.39 0; 0 0.39 0; 0 0.39 0; 1 0 1; 1 0 1; 1 0 1]
box_colors = [255/256,215/256,0; 1 0 0; 0 0 1; 255/256,215/256,0; 1 0 0; 0 0 1]
box_labels = {'MS', 'V', 'M', 'MS', 'V', 'M'}
plot_xval46 = [ones(9, 1), 2*ones(9,1), 3*ones(9,1)]
plot_xval49 = [4*ones(10, 1), 5*ones(10, 1), 6*ones(10, 1)]
figure;
hold on
plot(plot_xval46', [propHCROIs(st46_keep, 1), propHCROIs(st46_keep, 2), propHCROIs(st46_keep, 3)]', 'o', 'Color', [0 0.6 0], 'LineWidth', 0.1)
plot(plot_xval46', [propHCROIs(st46_keep, 1), propHCROIs(st46_keep, 2), propHCROIs(st46_keep, 3)]', 'Color', [0 0.6 0], 'LineWidth', 0.1)
plot(plot_xval49', [propHCROIs(st49_keep, 1), propHCROIs(st49_keep, 2), propHCROIs(st49_keep, 3)]', 'o', 'Color', [238/256, 130/256, 238/256], 'LineWidth', 0.1)
plot(plot_xval49', [propHCROIs(st49_keep, 1), propHCROIs(st49_keep, 2), propHCROIs(st49_keep, 3)]', 'Color', [238/256, 130/256, 238/256], 'LineWidth', 0.1)
h = boxplot(box_data, box_ids, 'BoxStyle', 'outline', 'Colors', box_colors, 'Labels', box_labels, 'PlotStyle', 'traditional')
set(h,{'linew'},{2})
hold off
xlabel('stage 46                    stage 49')
ylabel('prop HC ROIs')
ylim([0 1])
title('Prop HC ROIs by stage and stimtype')
set(gca, 'FontSize', 20)
fig_filename = 'Prop HC ROIs by stage and stimtype 30pct thresh'
saveas(gcf, fig_filename, 'png')
saveas(gcf, fig_filename, 'epsc2')
            

           % light magenta [238/256, 130/256, 238/256]
           % light darkgreen [152/256, 251/256, 152/256]
%% How many ROIs are HC in multiple modalities?

for k = 1:length(boot_xcorr)
    if ~isempty(boot_xcorr(k).highcorrROIs)
        HCROI_overlap(k,1) = length(intersect(boot_xcorr(k).highcorrROIs{1}, boot_xcorr(k).highcorrROIs{2})); %MS and V
        HCROI_overlap(k,2) = length(intersect(boot_xcorr(k).highcorrROIs{1}, boot_xcorr(k).highcorrROIs{3})); %MS and M
        HCROI_overlap(k,3) = length(intersect(boot_xcorr(k).highcorrROIs{2}, boot_xcorr(k).highcorrROIs{3})); %M and V
        HCROI_overlap(k,4) = length(intersect(intersect(boot_xcorr(k).highcorrROIs{1}, boot_xcorr(k).highcorrROIs{2}), boot_xcorr(k).highcorrROIs{3})); %MS, M and V
    end
end

% values for venn diagram
for k = 1:length(boot_xcorr)
    if ~isempty(boot_xcorr(k).highcorrROIs)
        for t = 1:length(boot_xcorr(k).highcorrROIs)
            total(t) = length(boot_xcorr(k).highcorrROIs{t});
        end
        Venn_results(k,1) = total(1) - (boot_xcorr(k).HCROI_overlap(1) + boot_xcorr(k).HCROI_overlap(2) - boot_xcorr(k).HCROI_overlap(4)); % MS only
        Venn_results(k,2) = total(2) - (boot_xcorr(k).HCROI_overlap(1) + boot_xcorr(k).HCROI_overlap(3) - boot_xcorr(k).HCROI_overlap(4)); % V only
        Venn_results(k,3) = total(3) - (boot_xcorr(k).HCROI_overlap(2) + boot_xcorr(k).HCROI_overlap(3) - boot_xcorr(k).HCROI_overlap(4)); % M only
        Venn_results(k,4) = boot_xcorr(k).HCROI_overlap(1) - boot_xcorr(k).HCROI_overlap(4); % MS and V only
        Venn_results(k,5) = boot_xcorr(k).HCROI_overlap(2) - boot_xcorr(k).HCROI_overlap(4); % MS and M only
        Venn_results(k,6) = boot_xcorr(k).HCROI_overlap(3) - boot_xcorr(k).HCROI_overlap(4); % M and V only
        Venn_results(k,7) = boot_xcorr(k).HCROI_overlap(4); % HC in all stimtypes
        Venn_results(k,8) = length(allData(k).respROIs_ST); % total responding ROIs
        Venn_results(k,9) = sum(Venn_results(k,1:7)); % total number of HC ROIS
        clear('total')
    end
end

% proportional HC ROIs for Venn Diagram
st46_HCvenn = Venn_results(st46_keep, 1:7) ./ repmat(Venn_results(st46_keep, 8), [1, 7]) 
st49_HCvenn = Venn_results(st49_keep, 1:7) ./ repmat(Venn_results(st49_keep, 8), [1, 7]) 
% mean for each category
mean_HCvenn(1, :) = mean(st46_HCvenn)
mean_HCvenn(2, :) = mean(st49_HCvenn)

% stat test it (Carlos recommends ANOVA, used anovan (for unbalanced design
% and 2 ways
stage = [repmat(46, [1 7*length(st46_HCvenn)]), repmat(49, [ 1, 7*length(st49_HCvenn)])]
stim_list = {'MS', 'V', 'M', 'MS_V', 'MS_M', 'M_V', 'MS_V_M'}
stimtype = [repmat(stim_list, [1 length(st46_HCvenn)]), repmat(stim_list, [1 length(st49_HCvenn)])];
Venn_data = [reshape(st46_HCvenn', [], 1); reshape(st49_HCvenn', [], 1)]
[p tbl stats] = anovan(Venn_data, {stage, stimtype}, 'varnames', {'stage', 'stimtype'})
multcompare(stats, 'Dimension', 2)

% not significant, p = 0.6 for stage and p = 0.056 for stimtype

% reduce to only the mains - MS, V and M
sm_stage = [repmat(46, [1 3*length(st46_HCvenn)]), repmat(49, [ 1, 3*length(st49_HCvenn)])]
sm_stimlist = {'MS', 'V', 'M'}
sm_stimtype = [repmat(sm_stimlist, [1 length(st46_HCvenn)]), repmat(sm_stimlist, [1 length(st49_HCvenn)])];
sm_Venn_data = [reshape(st46_HCvenn(:, 1:3)', [], 1); reshape(st49_HCvenn(:, 1:3)', [], 1)]
[p_sm, tbl_sm, stats_sm] = anovan(sm_Venn_data, {sm_stage, sm_stimtype}, 'varnames', {'stage', 'stimtype'})

% still not significant p = 0.66 for stage and 0.32 for stimtype


%% Anything different between HC and NHC? 
% by tad
% peak size, number of responses, response reliability, MS index, MS
% timing, MS by response reliability

% list of NHC but responding ROIs
for k = 1:length(boot_xcorr)
    if ~isempty(boot_xcorr(k).highcorrROIs)
        for t = 1:length(boot_xcorr(k).highcorrROIs)
            boot_xcorr(k).notHCROIs{t} = setdiff(allData(k).respROIs_ST, allData(k).respROIs_ST(boot_xcorr(k).highcorrROIs{t}));
        end
    end
end

% get basic properties
for k = 1:length(allData)
    stimmask = get_stimmask(allData(k).stimorder);
    allData(k).avg_meanpeak = mean_by_stimtype2(allData(k).meanpeak_bytrial_ST, stimmask);
    allData(k).avg_meanpeakloc = mean_by_stimtype2(allData(k).peakloc_bytrial_ST, stimmask);
    allData(k).MSenh_peak = calc_MSenhancement(allData(k).avg_meanpeak);
    allData(k).MSenh_peakloc = calc_MSenhancement(allData(k).avg_meanpeakloc);
end

%% stat test properties 
for k = 1:length(boot_xcorr)
    if ~isempty(boot_xcorr(k).highcorrROIs)
        for s = 1:4
            if ~isempty(boot_xcorr(k).highcorrROIs{s})
            roiH = allData(k).respROIs_ST(boot_xcorr(k).highcorrROIs{s});
            roiN = boot_xcorr(k).notHCROIs{s};
            [compareHC_P(k, s, 1), compareHC_H(k, s, 1), compareHC_stats{k, s, 1}] = ranksum(allData(k).avg_meanpeak(s, roiH), allData(k).avg_meanpeak(s, roiN)); %peak
            [compareHC_P(k, s, 2), compareHC_H(k, s, 2), compareHC_stats{k, s, 2}] = ranksum(allData(k).avg_meanpeakloc(s, roiH), allData(k).avg_meanpeakloc(s, roiN)); %peak location
            [compareHC_P(k, s, 3), compareHC_H(k, s, 3), compareHC_stats{k, s, 3}] = ranksum(allData(k).MSenh_peak(:, roiH), allData(k).MSenh_peak(:, roiN)); %MS index - peak
            [compareHC_P(k, s, 4), compareHC_H(k, s, 4), compareHC_stats{k, s, 4}] = ranksum(allData(k).MSenh_peakloc(:, roiH), allData(k).MSenh_peakloc(:, roiN)); % MS index - peak location
                % convert actual ROIs to index into respROIs for MSI_resprel
                [C roiNI] = intersect(allData(k).respROIs_ST, roiN);
                [C roiHI] = intersect(allData(k).respROIs_ST, roiH);
            [compareHC_P(k, s, 5), compareHC_H(k, s, 5), compareHC_stats{k, s, 5}] = ranksum(allData(k).MSI_resprel(roiHI), allData(k).MSI_resprel(roiNI)); % MS index - resp reliability
            [compareHC_P(k, s, 6), compareHC_H(k, s, 6), compareHC_stats{k, s, 6}] = ranksum(allData(k).sum_responses_ST(roiH), allData(k).sum_responses_ST(roiN)); % response reliability
             clear('roiH', 'roiN', 'roiHI', 'roiNI') 
            end
        end
    end
end

% how many sigs for each property?
total_sig = squeeze(sum(compareHC_H))
% LOL it's maximum 6 H =1, so nope, no differences. 

%% Do a PCA on the ROIs to see how HC and NHC pan out 

% first need to re-calculate a few missing properties
% peak variability, onset time with stim onset taken into account, primary
% modality, response reliability

for k = 1:length(allData) 
    if ~isempty(allData(k).respROIs_ST)
        

% get all HC and NHC into a matrix together with the relevant data points
idx = 1
for k = 1:length(boot_xcorr) %removes k = 29, which has no HC rois
    if ~isempty(allData(k).respROIs_ST)
        for r = 1:length(allData(k).respROIs_ST)
        allRespROIs(idx, 1) = allData(k).expnum;
        allRespROIs(idx, 2) = allData(k).respROIs_ST(r);
        allRespROIs(idx, 3) = allData(k).stage;
        allRespROIs(idx, 4).dff0 = allData(k).smoothed_trunc(allData(k).respROIs_ST(r), :);
        allRespROIs(idx, 5:7) = allData(k).avg_meanpeak(allData(k).respROIs_ST(r), 1:3);
        allRespROIs(idx, 8:10) = avg_onsettime
        allRespROIs(idx, 8:10) = std_onsettime
        allRespROIs(idx, 11) = allData(k).MSenh_peak(allData(k).respROIs_ST(r));
        allRespROIs(idx, 12) = allData(k).MSenh_peakloc(allData(k).respROIs_ST(r));
        
        allRespROIs(idx, 7) = ismember(boot_xcorr(k).highcorrROIs{1,1}); %MS HC
        allRespROIs(idx, 7) = ismember(boot_xcorr(k).highcorrROIs{1,2}); %V HC
        allRespROIs(idx, 7) = ismember(boot_xcorr(k).highcorrROIs{1,3}); %M HC
        %allRespROIs(idx, ).euclid_dist = allData(k).ROIdists(allData(k).respROIs_ST(r), :);
        %allRespROIs(idx, ).roiloc = allData(k).somaticROICenters{1,r}(1).Centroid;
        %allRespROIs(idx, ).stimorder = allData(k).stimorder;
        
        idx = idx + 1
    end
end
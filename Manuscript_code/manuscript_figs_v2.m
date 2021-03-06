%% Calcium Imaging Paper Figures Version 2
% started September 2018
% Start with new allData, saved as a struct. 
% D:\Torrey_calcium_imaging\2018August
% load('allData_workspace_7.mat')
% figures and workspaces are saved into D:\Torrey_calcium_imaging\manuscript_figures\Version_2

% downloaded boundedline, an open source matlab function for creating
% shaded error bars.
% https://github.com/kakearney/boundedline-pkg

%% Dataset details
% 29 Tadpoles in dataset. 10 stage 46 (9 keepers) and 18 stage 49 (10
% keepers)
% keepers = at least 25% of ROIs respond at least 1x in a stimulus trial. 
    % response = peak > 0.1, area >0 

folder = 'D:\Torrey_calcium_imaging\manuscript_figures\Version_2'    
    
%% Generate Data Values using smoothed and truncated data, for only good trials. 
% allData(k).smoothed_trunc is the df/f to use and allData(k).include is the good trials based on older analysis. 

% All fields generated from this smoothed and truncated data are designated
% _ST in the field name.

% Generate new peak, area and response values
for k = 1:length(allData)
    % Extract parameters for each trial
    [ allData(k).area_bytrial_ST ] = calc_area( allData(k).smoothed_trunc, 43 );
    [ allData(k).meanpeak_bytrial_ST, allData(k).peakloc_bytrial_ST] = calc_peak2( allData(k).smoothed_trunc, 5, 5);

    % Define response/no response
    [ allData(k).boolean_response_ST, allData(k).sum_responses_ST ] = get_respondingROIs3( allData(k).area_bytrial_ST, allData(k).meanpeak_bytrial_ST, allData(k).peakloc_bytrial_ST );
end

% recalculate trials to include
for k = 1:length(allData)
    allData(k).include2 = zeros(size(allData(k).meanpeak_bytrial_ST))
    for t = 1:size(allData(k).meanpeak_bytrial_ST,1)
        for r = 1:size(allData(k).meanpeak_bytrial_ST,2)
            if allData(k).meanpeak_bytrial_ST(t,r) < 5
                if allData(k).meanpeak_bytrial_ST(t,r) > -0.5
                    allData(k).include2(t,r) = 1;
                end
            end
        end
    end
end

% recalculate responding ROIs
for k = 1:length(allData)
    trs = find(allData(k).stimorder ~= 4);
    length(trs)
    tots = sum(allData(k).boolean_response_ST(:, trs),2);
    length(tots)
    allData(k).respROIs_ST = find(tots);
end

% are these roughly the same rOIs?
for k = 1:length(allData)
    ROI_cts(k, 1) = size(allData(k).resp_ROIs, 1)
    ROI_cts(k, 2) = length(allData(k).respROIs_ST)
    if ROI_cts(k, 1) > 0 && ROI_cts(k,2) > 0
    ROI_cts(k, 3) = length(intersect(allData(k).resp_ROIs(:, 1), allData(k).respROIs_ST))
    end
end
% mostly, yes. Although there's some added ROIs in exps that didn't have
% any in the previous iteration. 


%% PRepare workspace:
% get list of 46 and 49 tads
st49_tads = [];
st46_tads = [];
for t = 1:length(allData)
    if allData(t).stage == 49
    st49_tads = [st49_tads, t];
    elseif allData(t).stage == 46
    st46_tads = [st46_tads, t];
    end
end

st_ID = [46*ones(1,length(st46_tads)), 49*ones(1,length(st49_tads))]

% find good tadpoles
good_tads = zeros(1, length(allData))
for k = 1:length(allData)
    tmp = length(allData(k).respROIs_ST)
    tmp2 = size(allData(k).sum_resp_all_df_f0, 1)
    if tmp > 10
    %if tmp > (tmp2 * 0.2)
        good_tads(k) = 1
    end
end
keepers = find(good_tads)
st49_keep = intersect(st49_tads, keepers) % 10 total
st46_keep = intersect(st46_tads, keepers) % 9 total

%% Figure 1 Stage 46 tads have more generally responsive cells
% This code is copied from manuscript_figures_v1 but plots only have the 19
% keeper tadpoles included. Also moved data formats from allData{1,t} cell array of structs to 
% allData(k) struct 

%% Figure 1A St 46 tads have more generally responsive cells 
% changes from manuscript_figures_v1: 
    %st49_tads -> st49_keep 
    %allData{1,t} -> allData(t)
        %resp_ROIs -> respROIs_ST
        %ROIcenters -> ROICenters

% Proportion of ROIs that respond by tadpole
for t = 1:length(st49_keep)
    st49_proprespROIs(t) = length(allData(st49_keep(t)).respROIs_ST) / size(allData(st49_keep(t)).ROICenters, 1);
end

for t = 1:length(st46_keep)
    st46_proprespROIs(t) = length(allData(st46_keep(t)).respROIs_ST) / size(allData(st46_keep(t)).ROICenters, 1);
end

%make plot
% first put all values in 1 vector and create a grouping variable
st_ID = [46*ones(1,length(st46_keep)), 49*ones(1,length(st49_keep))]
all_proprespROIs = [st46_proprespROIs, st49_proprespROIs]

figure;
hold on
boxplot(all_proprespROIs, st_ID, 'Colors', 'gm')
plot(ones(length(st46_keep), 1), st46_proprespROIs, 'go')
plot(2*ones(length(st49_keep), 1), st49_proprespROIs, 'mo')
hold off 
xlabel('Stage')
ylabel('Proportion Responding ROIs')
set(gca,'FontSize',20)
set(findobj(gca,'type','line'),'linew',2)
fig_filename = '\figure_1\f1p1_proprespROIsbytad'
saveas(gcf, fullfile(folder, fig_filename), 'png')
saveas(gcf, fullfile(folder, fig_filename), 'epsc2')

% stat test (KS)
[h p stat] = kstest2(st46_proprespROIs, st49_proprespROIs)
% p = 0.0506


% % proportion responses for ROIs by tadpole
for t = 1:length(allData)
    tmp = allData(t).sum_responses_ST(allData(t).respROIs_ST) / length(allData(t).stimorder)
    prop_resp_trials(t) = mean(tmp)
end

%make plot
sorted_propresptrials = [prop_resp_trials(st46_keep), prop_resp_trials(st49_keep)];
figure;
hold on
boxplot(sorted_propresptrials, st_ID, 'Colors', 'gm')
plot(ones(length(st46_keep)), prop_resp_trials(st46_keep), 'go')
plot(2*ones(length(st49_keep)), prop_resp_trials(st49_keep), 'mo')
hold off 
xlabel('Stage')
ylabel('Proportion trials with response')
set(gca,'FontSize',20)
set(findobj(gca,'type','line'),'linew',2)
fig_filename = '\figure_1\f1p1_propresptrialssbyROIbytad'
saveas(gcf, fullfile(folder, fig_filename), 'png')
saveas(gcf, fullfile(folder, fig_filename), 'epsc2')

% stat test (KS)
[h p stat] = kstest2(prop_resp_trials(st46_keep), prop_resp_trials(st49_keep))
% p = 0.7475


% Response size for respROIs by tadpole
% % peak size for trials with resp for each resp_roi
% %     if boolean response, include peakbytrial2 in avg
% % avg peaks = per roi
% % avg rois = per tad


for t = 1:length(allData)
    avg_peaksize(t) = mean(mean(allData(t).avg_meanpeak(:,allData(t).respROIs_ST)))
end

%make plot
sorted_avg_peaksize = [avg_peaksize(st46_keep), avg_peaksize(st49_keep)];
figure;
hold on
boxplot(sorted_avg_peaksize, st_ID, 'Colors', 'gm')
plot(ones(length(st46_keep)), avg_peaksize(st46_keep), 'go')
plot(2*ones(length(st49_keep)), avg_peaksize(st49_keep), 'mo')
hold off 
xlabel('Stage')
ylabel('Peak Size')
set(gca,'FontSize',20)
set(findobj(gca,'type','line'),'linew',2)
fig_filename = '\figure_1\f1p1_avg_peaksize_bytad'
saveas(gcf, fullfile(folder, fig_filename), 'png')
saveas(gcf, fullfile(folder, fig_filename), 'epsc2')

% stat test (KS)
[h p stat] = kstest2(avg_peaksize(st46_keep), avg_peaksize(st49_keep))
% p = 0. 1273

%% Figure 2 Stage 49 tads have more MSI 
% MSIndex, MS timing, MS Resp Reliability 

% 2A Standard MSI (peak response)
for t = 1:length(allData)
    avg_MSInd_peak(t) = mean(mean(allData(t).MSenh_peak(:,allData(t).respROIs_ST)));
end

%make plot
sorted_avg_MSInd_peak = [avg_MSInd_peak(st46_keep), avg_MSInd_peak(st49_keep)];
figure;
hold on
boxplot(sorted_avg_MSInd_peak, st_ID, 'Colors', 'gm')
plot(ones(length(st46_keep)), avg_MSInd_peak(st46_keep), 'go')
plot(2*ones(length(st49_keep)), avg_MSInd_peak(st49_keep), 'mo')
hold off 
xlabel('Stage')
ylabel('MSInd')
set(gca,'FontSize',20)
set(findobj(gca,'type','line'),'linew',2)
fig_filename = '\figure_2\f2p1_avg_MSInd_peak_bytad'
saveas(gcf, fullfile(folder, fig_filename), 'png')
saveas(gcf, fullfile(folder, fig_filename), 'epsc2')

% stat test (KS)
[h p stat] = kstest2(avg_MSInd_peak(st46_keep), avg_MSInd_peak(st49_keep))
% p = 0.9485

% 2B MS Timing
for t = 1:length(allData)
    avg_MSInd_onsettime(t) = mean(mean(allData(t).MSenh_peakloc(:,allData(t).respROIs_ST)));
end

%make plot
sorted_avg_MSInd_onsettime = [avg_MSInd_onsettime(st46_keep), avg_MSInd_onsettime(st49_keep)];
figure;
hold on
boxplot(sorted_avg_MSInd_onsettime, st_ID, 'Colors', 'gm')
plot(ones(length(st46_keep)), avg_MSInd_onsettime(st46_keep), 'go')
plot(2*ones(length(st49_keep)), avg_MSInd_onsettime(st49_keep), 'mo')
hold off 
xlabel('Stage')
ylabel('MSInd')
set(gca,'FontSize',20)
set(findobj(gca,'type','line'),'linew',2)
fig_filename = '\figure_2\f2p1_avg_MSInd_onsettime_bytad'
saveas(gcf, fullfile(folder, fig_filename), 'png')
saveas(gcf, fullfile(folder, fig_filename), 'epsc2')

% stat test (KS)
[h p stat] = kstest2(avg_MSInd_onsettime(st46_keep), avg_MSInd_onsettime(st49_keep))
% p = 0.8646

% 2C MS reliability

for k = 1:length(allData)
    MS = find(allData(k).stimorder == 1);
    US = [find(allData(k).stimorder == 2), find(allData(k).stimorder == 3)];
    for r = 1:length(allData(k).respROIs_ST)
        Uval = sum(allData(k).boolean_response_ST(allData(k).respROIs_ST(r), US))
        Mval = sum(allData(k).boolean_response_ST(allData(k).respROIs_ST(r), MS))
        if Mval + Uval > 0
            allData(k).MSI_resprel(r) = Uval / (Mval + Uval) ; 
        else
            allData(k).MSI_resprel(r) = nan
        end
    end
    clear('MS', 'US')
end

% plot mean for each tad
for t = 1:length(allData)
    avg_MSInd(t) = nanmean(allData(t).MSI_resprel);
end

%make plot
sorted_avg_MSInd = [avg_MSInd(st46_keep), avg_MSInd(st49_keep)];
figure;
hold on
boxplot(sorted_avg_MSInd, st_ID, 'Colors', 'gm')
plot(ones(length(st46_keep)), avg_MSInd(st46_keep), 'go')
plot(2*ones(length(st49_keep)), avg_MSInd(st49_keep), 'mo')
hold off 
xlabel('Stage')
ylabel('MSInd')
set(gca,'FontSize',20)
set(findobj(gca,'type','line'),'linew',2)
fig_filename = '\figure_2\f2p3_avg_MSInd_responsereliability_bytad'
saveas(gcf, fullfile(folder, fig_filename), 'png')
saveas(gcf, fullfile(folder, fig_filename), 'epsc2')

% statistically different? 
% https://www.mathworks.com/help/stats/ranksum.html
[h, p, stat] = kstest2(avg_MSInd(st46_keep), avg_MSInd(st49_keep))
% not significant (p = 0.7044)

%% Figure 1 and 2 histograms to display with boxplots

% use histcounts to get data, then boundedline to plot mean and 95% CI
% (2*SD)

% figure 1, get data, (peak and latency for MS, V, M)
stims = [1 2 3 4]
for k = 1:length(allData)
    if ~isempty(allData(k).respROIs_ST)
        for s = 1:length(stims)
            
        end
    end
end

%% Figure 3: Correlation Methods; max correlation occurs at 0ms lag
% recalculate xcorr with truncated df/f
% replot xcorr by tad
% calculate xcorr with shuffled df/f
% plot shuffled xcorr by tad
% plot avg R vs response reliability
% plot distance vs R val to rule out spillover 

%% Recalculate xcorr using same code as in dissertation
% C:\Users\Torrey\Documents\GitHub\Calcium_Imaging_Analysis2\correlations\xcorr_cleaneddata
% contains code I worked from (needed to by converted to struct from cell
% array of structs

%%  First, find trials that meet "bad" criterion
for t = 1:length(allData)
    mask_peakU = allData(t).meanpeak_bytrial_ST < 5;
    mask_peakL = allData(t).meanpeak_bytrial_ST > -1;
    %mask_area = allData(t).area_bytrial_ST > -1;
    allData(t).trial_mask = mask_peakU & mask_peakL % mask_area;
    allData(t).ROI_oktrial_ct = sum(allData(t).trial_mask,2);
    clear('mask_peak', 'mask_area')
end

%how many trials are ok?
for t = 1:length(allData)
    tot = length(allData(t).stimorder)
    trialct(2,t) = mean(allData(t).ROI_oktrial_ct(allData(t).respROIs_ST) / tot)
end

%% re-run xcorr with the mask
% eliminate trials based on the mask, and concatenate only trials where we
% have both as good trials for the ROI pair

%% calculate correlation coefficient for responding ROIs using all trials
% maxlag defaults to 2N-1, with N = greater of the lengths of x and y
% therefore, to prevent crashing Matlab, set maxlag. I will set max lag to
% the length of one trial (so the maximum tested lag is that response 1
% starts at the beginning of the trial and response y starts at the end of
% the trial. 
% chose 'coeff' as the normalization option because it normalizes the
% sequence so that the autocorrelations at 0 lag equal 1.

% assign stim numbers to groups 
% 1 = MS, 2 = V, 3 = M, 4 = NS
stims = {[1, 5, 8, 9, 10, 11], [3, 7, 12], [2, 6], [4]}

for t = 1:length(allData)
    set_lag = length(allData(t).smoothed_trunc{1,1});
    for s = 1:length(stims)
    all_trials = ismember(allData(t).stimorder, stims{s});
        for r1 = 1:length(allData(t).respROIs_ST)
            for r2 = 1:length(allData(t).respROIs_ST)
                good_trials = allData(t).trial_mask(allData(t).respROIs_ST(r1), :) & allData(t).trial_mask(allData(t).respROIs_ST(r2), :);
                trials = intersect(find(good_trials), find(all_trials));
                allData(t).incl_trials{s, r1, r2} = trials;
                num_tr = length(trials);
                if ~isempty(trials)
                    tmp_data = [];
                    for tr = 1:length(trials)
                        tmp_data = [tmp_data; allData(t).smoothed_trunc{[allData(t).respROIs_ST(r1), allData(t).respROIs_ST(r2)], trials(tr)}];
                    end
                    allData(t).OKtrialpairs{s,r1,r2} = tmp_data;
                    clear('include')
                    [allData(t).respROIdff0_R{s, r1, r2}, allData(t).respROIdff0_lag_{s, r1, r2}] = xcorr(allData(t).OKtrialpairs{s, r1,r2}(:,1), allData(t).OKtrialpairs{s, r1,r2}(:,2), set_lag, 'coeff');
                else
                    fprintf('no data in t=%d, r1=%d r2=%d \n', t, r1, r2)
                end
            end
        end
        clear('trials')
    end
end

%% Find maxR and lag maxR for each xcorr
for t = 1:length(allData)
    for s = 1:length(stims)
            for r1 = 1:size(allData(t).respROIdff0_R, 2)
                for r2 = 1:size(allData(t).respROIdff0_R, 3)
                    if ~isempty(allData(t).respROIdff0_R{s, r1, r2})
                    [allData(t).respROIdff0_maxR_sq(s, r1, r2), allData(t).respROIdff0_maxRlag_sq(s,r1, r2)] = max(allData(t).respROIdff0_R{s,r1, r2}(:,:));
                    else
                       allData(t).respROIdff0_maxR_sq(s,r1, r2) = 0;
                       allData(t).respROIdff0_maxRlag_sq(s,r1, r2) = nan;
                    end
                end
            end
    end
end

%% Find R at 0lag
for t = 1:length(allData)
    loc = length(allData(t).smoothed_trunc{1,1}) + 1;
    for s = 1:length(stims)
            for r1 = 1:size(allData(t).respROIdff0_R, 2)
                for r2 = 1:size(allData(t).respROIdff0_R, 3)
                    if ~isempty(allData(t).respROIdff0_R{s, r1, r2})
                        allData(t).respROIdff0_0lagR_sq(s, r1, r2) = max(allData(t).respROIdff0_R{s,r1, r2}(loc,:));
                    else
                        allData(t).respROIdff0_0lagR_sq(s,r1, r2) = 0;
                       
                    end
                end
            end
    end
end

%% Plot xcorr [maxR, lag of maxR, 0lagR]
% sort from high to low maxR in the MS case, then use that sort for
% everything

% sort values
for t = 1:length(allData)
    if ~isempty(allData(t).respROIs_ST)
        for s = 1:length(stims)
            avg = mean(allData(t).respROIdff0_maxR_sq(s, :, :));
            [B, I] = sort(avg, 'descend');
            allData(t).respROI_maxR_I{s} = I;
            allData(t).maxR_sorted(s, :,:) = sort_ROIs(I, squeeze(allData(t).respROIdff0_maxR_sq(s, :,:)));
            allData(t).maxRlag_sorted(s, :,:) = sort_ROIs(I, squeeze(allData(t).respROIdff0_maxRlag_sq(s, :,:)));
            allData(t).lag0R_sorted(s, :,:) = sort_ROIs(I, squeeze(allData(t).respROIdff0_0lagR_sq(s, :,:)));
        end
    end
end

% make plots maxR
for t = 1:length(allData)
    if ~isempty(allData(t).respROIs_ST)
        for s = 1:length(stims)
            plot_xcorr(squeeze(allData(t).maxR_sorted(s, :, :)), 'hot')
            title(sprintf('Tad%d(t=%d) sorted maxR sm_trunc stim %d', allData(t).expnum, t, s))
            fig_filename = sprintf('Tad%d(t=%d) sorted maxR sm_trunc stim %d', allData(t).expnum, t, s)
            %saveas(gcf, fig_filename, 'epsc2')
            saveas(gcf, fig_filename, 'png')
            close;
        end
    end
end

% make plots 0lag R
for t = 1:length(allData)
    if ~isempty(allData(t).respROIs_ST)
        for s = 1:length(stims)
            plot_xcorr(squeeze(allData(t).lag0R_sorted(s, :, :)), 'hot')
            title(sprintf('Tad%d(t=%d) sorted lag 0 R sm_trunc stim %d', allData(t).expnum, t, s))
            fig_filename = sprintf('Tad%d(t=%d) sorted lag 0 R sm_trunc stim %d', allData(t).expnum, t, s)
            %saveas(gcf, fig_filename, 'epsc2')
            saveas(gcf, fig_filename, 'png')
            close;
        end
    end
end


% make plots maxR lag
for t = 1:length(allData)
    if ~isempty(allData(t).respROIs_ST)
        for s = 1:length(stims)
            plot_xcorr(squeeze(allData(t).maxRlag_sorted(s, :, :)), 'jet')
            title(sprintf('Tad%d(t=%d) sorted maxRlag sm_trunc stim %d', allData(t).expnum, t, s))
            fig_filename = sprintf('Tad%d(t=%d) sorted maxRlag sm_trunc stim %d', allData(t).expnum, t, s)
            %saveas(gcf, fig_filename, 'epsc2')
            saveas(gcf, fig_filename, 'png')
            close;
        end
    end
end

%% Figure 3A: put diss figure 4 in as background 

%% Figure 3B
    
%% Shuffle the xcorr trials to bootstrap 
% Shuffle 1000 times. 
% Take the maxR from each and make a distribution. 
% Get confidence intervals of each ROI. 
% Calculate how many ROIs' actual maxR is outside confidence interval

% Create a new struct with the bootstrap, xcorr_boot

%%%%%%%%%%%%%%%%%% COME BACK TO THIS %%%%%%%%%%%%%%%%%%%%%


%% Average maxR vs distance
% plot by tadpole and get line of best fit 

% get distances between ROIs
for k = 1:length(allData)
    if ~isempty(allData(k).ROICenters)
        for r1 = 1:length(allData(k).respROIs_ST)
            rr1 = allData(k).respROIs_ST(r1);
            for r2 = 1:length(allData(k).respROIs_ST)
                rr2 = allData(k).respROIs_ST(r2);
                allData(k).ROIdists(r1, r2) = sqrt((allData(k).ROICenters(rr1,1) - allData(k).ROICenters(rr2,1))^2 + (allData(k).ROICenters(rr1,2) - allData(k).ROICenters(rr2,2))^2);
            end
        end
    end
end

%% plot distance vs maxR
stim_IDs = {'MS', 'V', 'M', 'N'};
%stim_colors = [
for k = 1:length(allData)
    if length(allData(k).respROIs_ST) > 1
        if ~isempty(allData(k).ROIdists)
    X_vals_long = reshape(allData(k).ROIdists, 1, []);
    X_vals = X_vals_long(X_vals_long ~= 0);
        for s = 1:length(stims)
            Y_vals_all = reshape(allData(k).respROIdff0_maxR_sq(s, :, :), 1, []);
            Y_vals = Y_vals_all(Y_vals_all < 0.9999);
            [p,err] = polyfit(X_vals,Y_vals,1);   % First order polynomial
            y_fit = polyval(p,X_vals,err);   % Values on a line
            y_dif = Y_vals - y_fit;          % y value difference (residuals)
            SSdif = sum(y_dif.^2);      % Sum square of difference
            SStot = (length(Y_vals)-1)*var(Y_vals);   % Sum square of y taken from variance
            rsq = 1-SSdif/SStot;        % Correlation 'r' value. If 1.0 the correlelation is perfect
            rsq_dist_maxR(k) = rsq;
            R_val_text = sprintf('R^2 = %.4f', rsq)
            figure;
            hold on
            scatter(X_vals, Y_vals)
            plot(X_vals, y_fit, 'LineWidth', 3)
            title(sprintf('Tad%d (t=%d) Euclid Dist vs maxR stim %s', allData(k).expnum, k, stim_IDs{s}))
            xlabel('Distance (pixels)')
            ylabel('max R val (norm)')
            annotation('textbox', [0.5 0.8 0.2 0.1], 'String', R_val_text)
            set(gca,'fontsize',20)
            fig_filename = sprintf('Tad%d (t=%d) Euclid Dist vs maxR stim %s', allData(k).expnum, k, stim_IDs{s});
            saveas(gcf, fig_filename, 'epsc2')
            saveas(gcf, fig_filename, 'png')
            clear('fig_filename', 'p', 'err', 'y_fit', 'Y_vals', 'SSdif', 'SStot', 'rsq', 'R_val_text')
            close;
        end
        end
    end
end

% Histogram of all R^2 values
histogram(rsq_dist_maxR, 10)
title('All R squared for Euclid Dist vs maxR')
xlabel('R^2')
ylabel('Counts')
set(gca,'fontsize',20)
fig_filename = 'All R squared for Euclid Dist vs maxR';
saveas(gcf, fig_filename, 'epsc2')
saveas(gcf, fig_filename, 'png')

%% Average maxR vs response reliability

% response reliability is sum_responses_ST / length(stimorder)
% avg maxR is the mean of the row for the ROI 

for k = 1:length(allData)
    if length(allData(k).respROIs_ST) > 1
        X_vals = allData(k).sum_responses_ST(allData(k).respROIs_ST) ./ length(allData(k).stimorder);
        for s = 1:length(stims)
            Y_vals = squeeze(mean(allData(k).respROIdff0_maxR_sq(s, :, :)));
            [p,err] = polyfit(X_vals,Y_vals,1);   % First order polynomial
            y_fit = polyval(p,X_vals,err);   % Values on a line
            y_dif = Y_vals - y_fit;          % y value difference (residuals)
            SSdif = sum(y_dif.^2);      % Sum square of difference
            SStot = (length(Y_vals)-1)*var(Y_vals);   % Sum square of y taken from variance
            rsq = 1-SSdif/SStot;        % Correlation 'r' value. If 1.0 the correlelation is perfect
            rsq_RespRel_maxR(k) = rsq;
            R_val_text = sprintf('R^2 = %.4f', rsq)
            figure;
            hold on
            scatter(X_vals, Y_vals)
            plot(X_vals, y_fit, 'LineWidth', 3)
            title(sprintf('Tad%d (t=%d) Response reliability vs mean maxR stim %s', allData(k).expnum, k, stim_IDs{s}))
            xlabel('Response Reliability')
            ylabel('max R val (norm)')
            annotation('textbox', [0.5 0.8 0.2 0.1], 'String', R_val_text)
            set(gca,'fontsize',20)
            fig_filename = sprintf('Tad%d (t=%d) Response REliability vs maxR stim %s', allData(k).expnum, k, stim_IDs{s});
            saveas(gcf, fig_filename, 'epsc2')
            saveas(gcf, fig_filename, 'png')
            clear('fig_filename', 'p', 'err', 'y_fit', 'Y_vals', 'SSdif', 'SStot', 'rsq', 'R_val_text')
            close;
        end
        clear('X_vals')
    end
end

% Histogram of all R^2 values
histogram(rsq_RespRel_maxR, 10)
title('All R squared for Response Reliability vs maxR')
xlabel('R^2')
ylabel('Counts')
set(gca,'fontsize',20)
fig_filename = 'All R squared for Response Reliability vs maxR';
saveas(gcf, fig_filename, 'epsc2')
saveas(gcf, fig_filename, 'png')

%% Get highcorr ROIs 

for k = 1:length(allData)
    if length(allData(k).respROIs_ST) > 1
        for s = 1:length(stims)
            for r = 1:size(allData(k).respROIdff0_maxR_sq, 2)
                for c = 1:size(allData(k).respROIdff0_maxR_sq, 3)
                    if allData(k).respROIdff0_maxR_sq(s, r, c) > 0.5
                        highcorr(r, c) = 1;
                    else
                        highcorr(r, c) = 0; 
                    end
                end
            end
            allData(k).respROIdff0_HC(s, :, :) = logical(highcorr);
            clear('highcorr')
        end
    end
end

%% get indexes of ROIs that are highly correlated, by ROI

for k = 1:length(allData)
    if length(allData(k).respROIs_ST) > 1
        for s = 1:length(stims)
            for r = 1:size(allData(k).respROIdff0_HC, 2)
                allData(k).respROIdff0_correlROIs{s, r} = find(allData(k).respROIdff0_HC(s, r, :));
            end
        end
    end
end

% determine the overlap in which other ROIs are correlated with that ROI
for k = 1:length(allData)
    if length(allData(k).respROIs_ST) > 1
        for s = 1:length(stims)
            for r = 1:size(allData(k).respROIdff0_correlROIs, 2)
                for c = 1:size(allData(k).respROIdff0_correlROIs, 2)
                    allData(k).correlROIs_int{s, r, c} = intersect(allData(k).respROIdff0_correlROIs{s, r}, allData(k).respROIdff0_correlROIs{s, c});
                end
            end
        end
    end
end
                
%% Determine the ROIs that are the same across all cells with significant
% overlap with a given ROI
% significan overlap = at least 1/6*total num ROIS in correlated_ROIs_alldff0_int

clear('lens', 'int', 'roi_count')
for k = 1:length(allData)
    if length(allData(k).respROIs_ST) > 10
        roi_count = (1/6) * length(allData(k).respROIs_ST);
        for s = 1:length(stims)
            for row = 1:size(allData(k).correlROIs_int, 2)
                for ct = 1:size(allData(k).correlROIs_int, 3)                
                    lens(ct) = length(allData(k).correlROIs_int{s, row, ct});
                end
                first_roi = find((lens > roi_count), 1);
                if isempty(first_roi) 
                    continue
                else
                    int = allData(k).correlROIs_int{s, row, first_roi};
                    for col = first_roi:size(allData(k).correlROIs_int, 3)
                        if lens(col) > roi_count
                            int = intersect(int, allData(k).correlROIs_int{s, row, col});
                        else
                            continue
                        end
                    end
                    allData(k).correlROIs_common{s,row} = int;
                end
                clear('lens', 'int')
            end
            
        end
    end
    clear('roi_count')
end

% Index ROI numbers to the actual ROIs 
for k = 1:length(allData)
    if length(allData(k).respROIs_ST) > 10
        roi_list = allData(k).respROIs_ST;
        for s = 1:length(stims)
            if size(allData(k).correlROIs_common, 1) >= s
                for i = 1:size(allData(k).correlROIs_common,2)
                    if ~isempty(allData(k).correlROIs_common{s,i})
                        allData(k).correlROIs_common_AROI{s,i} = roi_list(allData(k).correlROIs_common{s,i});
                    end
                end
%             else
%                 allData(k).correlROIs_common_AROI{s} = [];
            end
        end
        clear('roi_list')
    end
end

%% How many HighCorr ROIs do I have? 

% Get unique HC count for each tad
for k = 1:length(allData)
    if ~isempty(allData(k).correlROIs_common_AROI)
        for s = 1:length(stims)
            tmp = [];
            for q = 1:size(allData(k).correlROIs_common_AROI, 2)
                tmp = [tmp; allData(k).correlROIs_common_AROI{s, q}];
            end
            allData(k).uniqueHC_AROI{s} = unique(tmp);
        end
    end
end
                
% Get HC and NHC count, and percentage in a separate variable
for k = 1:length(allData)
    if ~isempty(allData(k).correlROIs_common_AROI)
        for s = 1:length(stims)
            highcorr_counts(k, s, 1) = length(allData(k).correlROIs_common_AROI{s}); %HC ROIs
            highcorr_counts(k, s, 2) = length(allData(k).respROIs_ST); %all ROIs
            highcorr_counts(k, s, 3) = highcorr_counts(k, s, 1) / highcorr_counts(k, s, 2); % proportion HC
        end
    else
        highcorr_counts(k, s, 1) = NaN;
    end
end

                
% Sort by Stage
for k = 1:length(allData)
    stage(k) = allData(k).stage;
end

propHC46_5 = [];
propHC49_5 = [];
for k = 1:length(stage)
    if length(allData(k).respROIs_ST) > 10
        if stage(k) == 46
            propHC46_5 = [propHC46_5; highcorr_counts(k, :, 3)];
        elseif stage(k) == 49
            propHC49_5 = [propHC49_5; highcorr_counts(k, :, 3)];
        else
            continue
        end
    end
end

                
%% Figure 4: Paired stats on xcorr

% calculate mean maxR by stage and stimtype

for k = 1:length(st46_keep)
    mean_maxR(k, 1:3) = squeeze(mean(mean(allData(st46_keep(k)).maxR_sorted(1:3, :, :), 2), 3));
end

st = length(st46_keep)
for k = 1:length(st49_keep)
    mean_maxR(k+st, 1:3) = squeeze(mean(mean(allData(st49_keep(k)).maxR_sorted(1:3, :, :), 2), 3));
end

% stat test
data = reshape(mean_maxR, 1, []) % preserves column ordering
stage = repmat([repmat(46, length(st46_keep), 1); repmat(49, length(st49_keep), 1)], 3, 1);
id = [1*ones(size(mean_maxR, 1), 1); 2*ones(size(mean_maxR, 1), 1); 3*ones(size(mean_maxR, 1), 1)];
[p tbl stats] = anovan(data, {stage, id})                
%p = 0.3191 for stage, 0.5512 for id                
                
                
% get data for boxplot into one vector of labels and 1 vector of data
box_ids = [ones(1, 9), 2*ones(1, 9), 3*ones(1, 9), 4*ones(1, 10), 5*ones(1, 10), 6*ones(1, 10)]
box_data = [mean_maxR(1:9, 1)', mean_maxR(1:9, 2)', mean_maxR(1:9, 3)', mean_maxR(10:19, 1)', mean_maxR(10:19, 2)', mean_maxR(10:19, 3)']
%box_colors = [0 0.39 0; 0 0.39 0; 0 0.39 0; 1 0 1; 1 0 1; 1 0 1]
box_colors = [255/256,215/256,0; 1 0 0; 0 0 1; 255/256,215/256,0; 1 0 0; 0 0 1]
box_labels = {'MS', 'V', 'M', 'MS', 'V', 'M'}
plot_xval46 = [ones(9, 1), 2*ones(9,1), 3*ones(9,1)]
plot_xval49 = [4*ones(10, 1), 5*ones(10, 1), 6*ones(10, 1)]
figure;
hold on
plot(plot_xval46(:,1:3), [mean_maxR(1:9, 1), mean_maxR(1:9, 2), mean_maxR(1:9, 3)], 'o', 'Color', [0 0.6 0], 'LineWidth', 0.1)
plot(plot_xval46(:, 1:3)', [mean_maxR(1:9, 1), mean_maxR(1:9, 2), mean_maxR(1:9, 3)]', 'Color', [0 0.6 0], 'LineWidth', 0.1)
plot(plot_xval49(:, 1:3), [mean_maxR(10:19, 1), mean_maxR(10:19, 2), mean_maxR(10:19, 3)], 'o', 'Color', [238/256, 130/256, 238/256], 'LineWidth', 0.1)
plot(plot_xval49(:, 1:3)', [mean_maxR(10:19, 1), mean_maxR(10:19, 2), mean_maxR(10:19, 3)]', 'Color', [238/256, 130/256, 238/256], 'LineWidth', 0.1)
h = boxplot(box_data, box_ids, 'BoxStyle', 'outline', 'Colors', box_colors, 'Labels', box_labels, 'PlotStyle', 'traditional')
set(h,{'linew'},{2})
hold off
xlabel('stage 46                    stage 49')
ylabel('mean correlation')
%ylim([0 1])
title('Mean corr by stage and stimtype')
set(gca, 'FontSize', 20)
fig_filename = 'Mean corr by stage and stimtype'
saveas(gcf, fig_filename, 'png')
saveas(gcf, fig_filename, 'epsc2')                
                
                
                
                
                
                
                
                
%% Implement data cleaning 

%remove trials with problems
%peak > 5, area < 0 
%if > 75% of trials are bad for a single ROI, eliminate ROI

% Using allData in D:\Torrey_calcium_imaging\compare_46-49\analysis_Feb 2018
load('xcorr_analysis_20180302')
%%%% This is SMOOTHED df/f0!! 

%% Find trials that meet "bad" criterion

% some of the areas are still type cell...
% some of the tadpoles have both regular and smoothed df/f0...
% why the fuck am I so inconsistent?????
% I want to use ONLY SMOOOTHED DATA FOR THIS ANALYSIS. 

%% Ok, let's start over. Only include smoothed data in the new structure. 
% original, inconsistent allData is renamed allData_old
% New, consistent allData is allData
% take only the smoothed df_f0 and associated files. Give them all the same
% names. 

for t = 1:length(allData_old)
    %fields that are the same 
    allData{1,t}.stimorder = allData_old{1,t}.stimorder;
    allData{1,t}.expnum = allData_old{1,t}.expnum;
    allData{1,t}.stage = allData_old{1,t}.stage;
    allData{1,t}.bathTBX = allData_old{1,t}.bathRBX; % TBX = tubacurarine, previous code put an R by mistake
    allData{1,t}.onset_time = allData_old{1,t}.onset_time;
    allData{1,t}.resp_ROIs = allData_old{1,t}.resp_ROIs;
    allData{1,t}.ROIcenters = allData_old{1,t}.ROIcenters;
    % fields that have different names in the tads with "smoothed"
    if  isfield(allData_old, 'smoothed')
        allData{1,t}.df_f0 = allData_old{1,t}.smoothed;
        if iscell(allData_old{1,t}.area_bytrial_sm)
            allData{1,t}.area_bytrial = cell2mat(allData_old{1,t}.area_bytrial_sm);
        else
            allData{1,t}.area_bytrial = allData_old{1,t}.area_bytrial_sm;
        end
        allData{1,t}.peak_bytrial = allData_old{1,t}.meanpeak_bytrial_sm;
    else
        allData{1,t}.df_f0 = allData_old{1,t}.df_f0;
        if iscell(allData_old{1,t}.area_bytrial)
            allData{1,t}.area_bytrial = cell2mat(allData_old{1,t}.area_bytrial);
        else
            allData{1,t}.area_bytrial = allData_old{1,t}.area_bytrial;
        end
        allData{1,t}.peak_bytrial = allData_old{1,t}.peak_bytrial;
    end
end

%% Make sure all df/f0 traces are the same dimensions (1 row, N cols)
% this will fix matrix concatenation issues with getting the data for xcorr
% later. 

for t= 1:length(allData)
    if size(allData{1,t}.df_f0{1,1}, 1) == 1
        for r = 1:size(allData{1,t}.df_f0,1)
            for tr = 1:size(allData{1,t}.df_f0,2)
                allData{1,t}.df_f0{r, tr} = allData{1,t}.df_f0{r, tr}';
            end
        end
    end
end


%% Now we can actually find trials that meet "bad" criterion
for t = 1:length(allData)
    mask_peak = allData{1,t}.peak_bytrial2 < 5;
    if iscell(allData{1,t}.area_bytrial)
        mask_area = cell2mat(allData{1,t}.area_bytrial) > 0;
    else
        mask_area = allData{1,t}.area_bytrial > 0;
    end
    allData{1,t}.trial_mask = mask_peak; %& mask_area;
    allData{1,t}.ROI_oktrial_ct = sum(allData{1,t}.trial_mask,2);
    clear('mask_peak', 'mask_area')
end

%how many trials are ok?
for t = 1:length(allData)
    tot = length(allData{1,t}.stimorder)
    trialct(2,t) = mean(allData{1,t}.ROI_oktrial_ct(allData{1,t}.resp_ROIs) / tot)
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

for t = 1:length(allData)
    set_lag = length(allData{1,t}.df_f0{1,1});
    for r1 = 1:length(allData{1,t}.resp_ROIs)
        for r2 = 1:length(allData{1,t}.resp_ROIs)
            include = allData{1,t}.trial_mask(allData{1,t}.resp_ROIs(r1), :) & allData{1,t}.trial_mask(allData{1,t}.resp_ROIs(r2), :);
            trials = find(include);
            if ~isempty(trials)
                tmp_data = [];
                for tr = 1:length(trials)
                    tmp_data = [tmp_data; allData{1,t}.df_f0{[allData{1,t}.resp_ROIs(r1), allData{1,t}.resp_ROIs(r2)], trials(tr)}];
                end
                allData{1,t}.OKtrialpairs_all{r1,r2} = tmp_data;
                clear('trials', 'include')
                [allData{1,t}.respROIdff0_R_all{r1, r2}, allData{1,t}.respROIdff0_lag_all{r1, r2}] = xcorr(allData{1,t}.OKtrialpairs_all{r1,r2}(:,1), allData{1,t}.OKtrialpairs_all{r1,r2}(:,2), set_lag, 'coeff');
            end
            end
    end
    clear('set_lag')
end

%% calculate correlation coefficient for responding ROIs using each modality separately


    %1 = multisensory high M / high V
    %2 = visual high/crash
    %3 = mechanosensory high
    %4 = no stimulus
    
    %5 = multisensory high M / low V
    %6 = visual low/scrambled crash
    %7 = mechanosensory low
    
    %8 = multisensory low M / high V
    %9 = multisensory low M / low V 
    
    %10 = multisensory med M / high V
    %11 = multisensory med M / low V
    %12 = mechanosensoy med
    MS_stim = [1, 5, 8, 9, 10, 11]
    M_stim = [3, 7, 12]
    V_stim = [2, 6]
    N_stim = [4]
    
% mutlisensory
for t = 1:length(allData)
    set_lag = length(allData{1,t}.df_f0{1,1});
    MS_trials = ismember(allData{1,t}.stimorder, MS_stim);
    for r1 = 1:length(allData{1,t}.resp_ROIs)
        for r2 = 1:length(allData{1,t}.resp_ROIs)
            good_trials = allData{1,t}.trial_mask(allData{1,t}.resp_ROIs(r1), :) & allData{1,t}.trial_mask(allData{1,t}.resp_ROIs(r2), :);
            include = good_trials & MS_trials;
            trials = find(include);
            num_tr = length(trials)
            if ~isempty(trials)
                tmp_data = [];
                for tr = 1:length(trials)
                    tmp_data = [tmp_data; allData{1,t}.df_f0{[allData{1,t}.resp_ROIs(r1), allData{1,t}.resp_ROIs(r2)], trials(tr)}];
                end
                allData{1,t}.OKtrialpairs_MS{r1,r2} = tmp_data;
                clear('trials', 'include')
                [allData{1,t}.respROIdff0_R_MS{r1, r2}, allData{1,t}.respROIdff0_lag_MS{r1, r2}] = xcorr(allData{1,t}.OKtrialpairs_MS{r1,r2}(:,1), allData{1,t}.OKtrialpairs_MS{r1,r2}(:,2), set_lag, 'coeff');
            end
        end
    end
end

% visual
for t = 1:length(allData)
    set_lag = length(allData{1,t}.df_f0{1,1});
    V_trials = ismember(allData{1,t}.stimorder, V_stim);
    for r1 = 1:length(allData{1,t}.resp_ROIs)
        for r2 = 1:length(allData{1,t}.resp_ROIs)
            good_trials = allData{1,t}.trial_mask(allData{1,t}.resp_ROIs(r1), :) & allData{1,t}.trial_mask(allData{1,t}.resp_ROIs(r2), :);
            include = good_trials & V_trials;
            trials = find(include);
            num_tr = length(trials)
            if ~isempty(trials)
                tmp_data = [];
                for tr = 1:length(trials)
                    tmp_data = [tmp_data; allData{1,t}.df_f0{[allData{1,t}.resp_ROIs(r1), allData{1,t}.resp_ROIs(r2)], trials(tr)}];
                end
                allData{1,t}.OKtrialpairs_V{r1,r2} = tmp_data;
                clear('trials', 'include')
                [allData{1,t}.respROIdff0_R_V{r1, r2}, allData{1,t}.respROIdff0_lag_V{r1, r2}] = xcorr(allData{1,t}.OKtrialpairs_V{r1,r2}(:,1), allData{1,t}.OKtrialpairs_V{r1,r2}(:,2), set_lag, 'coeff');
            end
        end
    end
end


% mechanosensory
for t = 1:length(allData)
    set_lag = length(allData{1,t}.df_f0{1,1});
    M_trials = ismember(allData{1,t}.stimorder, M_stim);
    for r1 = 1:length(allData{1,t}.resp_ROIs)
        for r2 = 1:length(allData{1,t}.resp_ROIs)
            good_trials = allData{1,t}.trial_mask(allData{1,t}.resp_ROIs(r1), :) & allData{1,t}.trial_mask(allData{1,t}.resp_ROIs(r2), :);
            include = good_trials & M_trials;
            trials = find(include);
            num_tr = length(trials);
            if ~isempty(trials)
                tmp_data = [];
                for tr = 1:length(trials)
                    tmp_data = [tmp_data; allData{1,t}.df_f0{[allData{1,t}.resp_ROIs(r1), allData{1,t}.resp_ROIs(r2)], trials(tr)}];
                end
                allData{1,t}.OKtrialpairs_M{r1,r2} = tmp_data;
                clear('trials', 'include')
                [allData{1,t}.respROIdff0_R_M{r1, r2}, allData{1,t}.respROIdff0_lag_M{r1, r2}] = xcorr(allData{1,t}.OKtrialpairs_M{r1,r2}(:,1), allData{1,t}.OKtrialpairs_M{r1,r2}(:,2), set_lag, 'coeff');
            end
        end
    end
%     %%%% If the last ROI has no data,OKtrialpairs ends up 1 unit short.
%     %%%% This code adds an empty row.
%     if isfield(allData{1,t}, 'OKtrialpairs_M')
%         if size(allData{1,t}.OKtrialpairs_M, 1) ~= length(allData{1,t}.resp_ROIs)
%             allData{1,t}.OKtrialpairs_M{length(allData{1,t}.resp_ROIs), 1} = []
%         end
%         if size(allData{1,t}.OKtrialpairs_M, 2) ~= length(allData{1,t}.resp_ROIs)
%             allData{1,t}.OKtrialpairs_M{1, length(allData{1,t}.resp_ROIs)} = []
%         end
%         for r1 = 1:length(allData{1,t}.resp_ROIs)
%             for r2 = 1:length(allData{1,t}.resp_ROIs)
%                 if ~isempty(allData{1,t}.OKtrialpairs_M{r1,r2})
%                 [allData{1,t}.respROIdff0_R_M{r1, r2}, allData{1,t}.respROIdff0_lag_M{r1, r2}] = xcorr(allData{1,t}.OKtrialpairs_M{r1,r2}(:,1), allData{1,t}.OKtrialpairs_M{r1,r2}(:,2), set_lag, 'coeff');
%                 end
%             end
%         end
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%
end

% no stim 
for t = 1:length(allData)
    set_lag = length(allData{1,t}.df_f0{1,1});
    N_trials = ismember(allData{1,t}.stimorder, N_stim);
    for r1 = 1:length(allData{1,t}.resp_ROIs)
        for r2 = 1:length(allData{1,t}.resp_ROIs)
            good_trials = allData{1,t}.trial_mask(allData{1,t}.resp_ROIs(r1), :) & allData{1,t}.trial_mask(allData{1,t}.resp_ROIs(r2), :);
            include = good_trials & N_trials;
            trials = find(include);
            num_tr = length(trials)
            if ~isempty(trials)
                tmp_data = [];
                for tr = 1:length(trials)
                    tmp_data = [tmp_data; allData{1,t}.df_f0{[allData{1,t}.resp_ROIs(r1), allData{1,t}.resp_ROIs(r2)], trials(tr)}];
                end
                allData{1,t}.OKtrialpairs_N{r1,r2} = tmp_data;
                clear('trials', 'include')
                [allData{1,t}.respROIdff0_R_N{r1, r2}, allData{1,t}.respROIdff0_lag_N{r1, r2}] = xcorr(allData{1,t}.OKtrialpairs_N{r1,r2}(:,1), allData{1,t}.OKtrialpairs_N{r1,r2}(:,2), set_lag, 'coeff');
            end
        end
    end
end


%% Generate data for plots of each xcorr ([maxR, lag of maxR, 0 lag R] for [all, MS, V, M, N])

%% Find maxR and lag maxR for each xcorr

% all
for t = 1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_R_all')
        for r1 = 1:size(allData{1,t}.respROIdff0_R_all, 1)
            for r2 = 1:size(allData{1,t}.respROIdff0_R_all, 2)
                if ~isempty(allData{1,t}.respROIdff0_R_all{r1, r2})
                [allData{1,t}.respROIdff0_maxR_sq_all(r1, r2), allData{1,t}.respROIdff0_maxRlag_sq_all(r1, r2)] = max(allData{1,t}.respROIdff0_R_all{r1, r2}(:,:));
                else
                   allData{1,t}.respROIdff0_maxR_sq_all(r1, r2) = 0;
                   allData{1,t}.respROIdff0_maxRlag_sq_all(r1, r2) = nan;
                end
            end
        end
    end
end

% multisensory
for t = 1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_R_MS')
        for r1 = 1:size(allData{1,t}.respROIdff0_R_MS, 1)
            for r2 = 1:size(allData{1,t}.respROIdff0_R_MS, 2)
                if ~isempty(allData{1,t}.respROIdff0_R_MS{r1, r2})
                [allData{1,t}.respROIdff0_maxR_sq_MS(r1, r2), allData{1,t}.respROIdff0_maxRlag_sq_MS(r1, r2)] = max(allData{1,t}.respROIdff0_R_MS{r1, r2}(:,:));
                else
                   allData{1,t}.respROIdff0_maxR_sq_MS(r1, r2) = 0;
                   allData{1,t}.respROIdff0_maxRlag_sq_MS(r1, r2) = nan;
                end
            end
        end
    end
end

% visual
for t = 1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_R_V')
        for r1 = 1:size(allData{1,t}.respROIdff0_R_V, 1)
            for r2 = 1:size(allData{1,t}.respROIdff0_R_V, 2)
                if ~isempty(allData{1,t}.respROIdff0_R_V{r1, r2})
                [allData{1,t}.respROIdff0_maxR_sq_V(r1, r2), allData{1,t}.respROIdff0_maxRlag_sq_V(r1, r2)] = max(allData{1,t}.respROIdff0_R_V{r1, r2}(:,:));
                else
                   allData{1,t}.respROIdff0_maxR_sq_V(r1, r2) = 0;
                   allData{1,t}.respROIdff0_maxRlag_sq_V(r1, r2) = nan;
                end
            end
        end
    end
end

% mechanosensory
for t = 1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_R_M')
        for r1 = 1:size(allData{1,t}.respROIdff0_R_M, 1)
            for r2 = 1:size(allData{1,t}.respROIdff0_R_M, 2)
                if ~isempty(allData{1,t}.respROIdff0_R_M{r1, r2})
                [allData{1,t}.respROIdff0_maxR_sq_M(r1, r2), allData{1,t}.respROIdff0_maxRlag_sq_M(r1, r2)] = max(allData{1,t}.respROIdff0_R_M{r1, r2}(:,:));
                else
                   allData{1,t}.respROIdff0_maxR_sq_M(r1, r2) = 0;
                   allData{1,t}.respROIdff0_maxRlag_sq_M(r1, r2) = nan;
                end
            end
        end
    end
end

% no stim
for t = 1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_R_N')
        for r1 = 1:size(allData{1,t}.respROIdff0_R_N, 1)
            for r2 = 1:size(allData{1,t}.respROIdff0_R_N, 2)
                if ~isempty(allData{1,t}.respROIdff0_R_N{r1, r2})
                [allData{1,t}.respROIdff0_maxR_sq_N(r1, r2), allData{1,t}.respROIdff0_maxRlag_sq_N(r1, r2)] = max(allData{1,t}.respROIdff0_R_N{r1, r2}(:,:));
                else
                   allData{1,t}.respROIdff0_maxR_sq_N(r1, r2) = 0;
                   allData{1,t}.respROIdff0_maxRlag_sq_N(r1, r2) = nan;
                end
            end
        end
    end
end


%% Find xcorr R at 0 lag

% all
for t = 1:length(allData)
    loc = length(allData{1,t}.df_f0{1,1}) + 1;
    if isfield(allData{1,t}, 'respROIdff0_R_all')
        for r1 = 1:size(allData{1,t}.respROIdff0_R_all, 1)
            for r2 = 1:size(allData{1,t}.respROIdff0_R_all, 2)
                if ~isempty(allData{1,t}.respROIdff0_R_all{r1, r2})
                    allData{1,t}.respROIdff0_0lagR_sq_all(r1, r2) = allData{1,t}.respROIdff0_R_all{r1, r2}(loc,:);
                else
                    allData{1,t}.respROIdff0_0lagR_sq_all(r1, r2) = 0;
                end
            end
        end
    end
end

%multisensory
for t = 1:length(allData)
    loc = length(allData{1,t}.df_f0{1,1}) + 1;
    if isfield(allData{1,t}, 'respROIdff0_R_MS')
        for r1 = 1:size(allData{1,t}.respROIdff0_R_MS, 1)
            for r2 = 1:size(allData{1,t}.respROIdff0_R_MS, 2)
                if ~isempty(allData{1,t}.respROIdff0_R_MS{r1, r2})
                    allData{1,t}.respROIdff0_0lagR_sq_MS(r1, r2) = allData{1,t}.respROIdff0_R_MS{r1, r2}(loc,:);
                else
                    allData{1,t}.respROIdff0_0lagR_sq_MS(r1, r2) = 0;
                end
            end
        end
    end
end

%visual
for t = 1:length(allData)
    loc = length(allData{1,t}.df_f0{1,1}) + 1;
    if isfield(allData{1,t}, 'respROIdff0_R_V')
        for r1 = 1:size(allData{1,t}.respROIdff0_R_V, 1)
            for r2 = 1:size(allData{1,t}.respROIdff0_R_V, 2)
                if ~isempty(allData{1,t}.respROIdff0_R_V{r1, r2})
                    allData{1,t}.respROIdff0_0lagR_sq_V(r1, r2) = allData{1,t}.respROIdff0_R_V{r1, r2}(loc,:);
                else
                    allData{1,t}.respROIdff0_0lagR_sq_V(r1, r2) = 0;
                end
            end
        end
    end
end

%mechanosensory
for t = 1:length(allData)
    loc = length(allData{1,t}.df_f0{1,1}) + 1;
    if isfield(allData{1,t}, 'respROIdff0_R_M')
        for r1 = 1:size(allData{1,t}.respROIdff0_R_M, 1)
            for r2 = 1:size(allData{1,t}.respROIdff0_R_M, 2)
                if ~isempty(allData{1,t}.respROIdff0_R_M{r1, r2})
                    allData{1,t}.respROIdff0_0lagR_sq_M(r1, r2) = allData{1,t}.respROIdff0_R_M{r1, r2}(loc,:);
                else
                    allData{1,t}.respROIdff0_0lagR_sq_M(r1, r2) = 0;
                end
            end
        end
    end
end

%no stim
for t = 1:length(allData)
    loc = length(allData{1,t}.df_f0{1,1}) + 1;
    if isfield(allData{1,t}, 'respROIdff0_R_N')
        for r1 = 1:size(allData{1,t}.respROIdff0_R_N, 1)
            for r2 = 1:size(allData{1,t}.respROIdff0_R_N, 2)
                if ~isempty(allData{1,t}.respROIdff0_R_N{r1, r2})
                    allData{1,t}.respROIdff0_0lagR_sq_N(r1, r2) = allData{1,t}.respROIdff0_R_N{r1, r2}(loc,:);
                else
                    allData{1,t}.respROIdff0_0lagR_sq_N(r1, r2) = 0;
                end
            end
        end
    end
end

%% Get highcorr cells by locating all maxR > 0.5

% based on all df/f0
for t=1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_R_all')
        if length(allData{1,t}.resp_ROIs) > 0
            for r = 1:size(allData{1,t}.respROIdff0_maxR_sq_all, 1)
                for c = 1:size(allData{1,t}.respROIdff0_maxR_sq_all, 2)
                    if allData{1,t}.respROIdff0_maxR_sq_all(r,c) > 0.5
                        highcorr(r, c) = 1;
                    else
                        highcorr(r, c) = 0;
                    end
                end
            end
            allData{1,t}.respROIdff0_highcorr = logical(highcorr)
            clear('highcorr')
        end
    end
end

% based on multi df/f0
for t=1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_R_MS')
        if length(allData{1,t}.resp_ROIs) > 0
            for r = 1:size(allData{1,t}.respROIdff0_maxR_sq_MS, 1)
                for c = 1:size(allData{1,t}.respROIdff0_maxR_sq_MS, 2)
                    if allData{1,t}.respROIdff0_maxR_sq_MS(r,c) > 0.5
                        highcorr(r, c) = 1;
                    else
                        highcorr(r, c) = 0;
                    end
                end
            end
            allData{1,t}.respROIdff0_highcorr_MS = logical(highcorr)
            clear('highcorr')
        end
    end
end            
            
% based on vis df/f0
for t=1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_R_V')
        if length(allData{1,t}.resp_ROIs) > 0
            for r = 1:size(allData{1,t}.respROIdff0_maxR_sq_V, 1)
                for c = 1:size(allData{1,t}.respROIdff0_maxR_sq_V, 2)
                    if allData{1,t}.respROIdff0_maxR_sq_V(r,c) > 0.5
                        highcorr(r, c) = 1;
                    else
                        highcorr(r, c) = 0;
                    end
                end
            end
            allData{1,t}.respROIdff0_highcorr_V = logical(highcorr)
            clear('highcorr')
        end
    end
end 
                
% based on mech df/f0
for t=1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_R_M')
        if length(allData{1,t}.resp_ROIs) > 0
            for r = 1:size(allData{1,t}.respROIdff0_maxR_sq_M, 1)
                for c = 1:size(allData{1,t}.respROIdff0_maxR_sq_M, 2)
                    if allData{1,t}.respROIdff0_maxR_sq_M(r,c) > 0.5
                        highcorr(r, c) = 1;
                    else
                        highcorr(r, c) = 0;
                    end
                end
            end
            allData{1,t}.respROIdff0_highcorr_M = logical(highcorr)
            clear('highcorr')
        end
    end
end 

% based on no stim df/f0
for t=1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_R_N')
        if length(allData{1,t}.resp_ROIs) > 0
            for r = 1:size(allData{1,t}.respROIdff0_maxR_sq_N, 1)
                for c = 1:size(allData{1,t}.respROIdff0_maxR_sq_N, 2)
                    if allData{1,t}.respROIdff0_maxR_sq_N(r,c) > 0.5
                        highcorr(r, c) = 1;
                    else
                        highcorr(r, c) = 0;
                    end
                end
            end
            allData{1,t}.respROIdff0_highcorr_N = logical(highcorr)
            clear('highcorr')
        end
    end
end 

%% Now assign ROIs to groups based on high correlations
% % % % First add empty rows to any tads where the last ROI has no data
% % % 
% % %     if isfield(allData{1,t}, 'respROIdff0_R_all')
% % %         if size(allData{1,t}.OKtrialpairs_M, 1) ~= length(allData{1,t}.resp_ROIs)
% % %             allData{1,t}.OKtrialpairs_M{length(allData{1,t}.resp_ROIs), 1} = []
% % %         end
% % %         if size(allData{1,t}.OKtrialpairs_M, 2) ~= length(allData{1,t}.resp_ROIs)
% % %             allData{1,t}.OKtrialpairs_M{1, length(allData{1,t}.resp_ROIs)} = []
% % %         end
% % %         for r1 = 1:length(allData{1,t}.resp_ROIs)
% % %             for r2 = 1:length(allData{1,t}.resp_ROIs)
% % %                 if ~isempty(allData{1,t}.OKtrialpairs_M{r1,r2})
% % %                 [allData{1,t}.respROIdff0_R_M{r1, r2}, allData{1,t}.respROIdff0_lag_M{r1, r2}] = xcorr(allData{1,t}.OKtrialpairs_M{r1,r2}(:,1), allData{1,t}.OKtrialpairs_M{r1,r2}(:,2), set_lag, 'coeff');
% % %                 end
% % %             end
% % %         end
% % %     end
% % %     %%%%%%%%%%%%%%%%%%%%%%%

% get indexes of ROIs that are highly correlated, by ROI
for t = 1:length(allData)
   if isfield(allData{1,t}, 'respROIdff0_R_all')
        if length(allData{1,t}.resp_ROIs) > 0
            for r = 1:size(allData{1,t}.respROIdff0_highcorr_MS,1)
                allData{1,t}.correlated_ROIs_alldff0{1,r} = find(allData{1,t}.respROIdff0_highcorr(r,:));
                allData{1,t}.correlated_ROIs_dff0_MS{1,r} = find(allData{1,t}.respROIdff0_highcorr_MS(r,:));
                allData{1,t}.correlated_ROIs_dff0_V{1,r} = find(allData{1,t}.respROIdff0_highcorr_V(r,:));
                allData{1,t}.correlated_ROIs_dff0_M{1,r} = find(allData{1,t}.respROIdff0_highcorr_M(r,:));
                allData{1,t}.correlated_ROIs_dff0_N{1,r} = find(allData{1,t}.respROIdff0_highcorr_N(r,:));
            end
        end
   end
end

% Determine overlap in which other ROIs are correlated with that ROI over all other ROIs
for t = 1:length(allData)
   if isfield(allData{1,t}, 'respROIdff0_R_all')
        if length(allData{1,t}.resp_ROIs) > 0
            for r = 1:size(allData{1,t}.correlated_ROIs_dff0_MS,2)
                for c = 1:size(allData{1,t}.correlated_ROIs_dff0_MS,2)
                    allData{1,t}.correlated_ROIs_alldff0_int{r,c} = intersect(allData{1,t}.correlated_ROIs_alldff0{1,r}, allData{1,t}.correlated_ROIs_alldff0{1,c});
                    allData{1,t}.correlated_ROIs_dff0_MS_int{r,c} = intersect(allData{1,t}.correlated_ROIs_dff0_MS{1,r}, allData{1,t}.correlated_ROIs_dff0_MS{1,c});
                    allData{1,t}.correlated_ROIs_dff0_V_int{r,c} = intersect(allData{1,t}.correlated_ROIs_dff0_V{1,r}, allData{1,t}.correlated_ROIs_dff0_V{1,c});
                    allData{1,t}.correlated_ROIs_dff0_M_int{r,c} = intersect(allData{1,t}.correlated_ROIs_dff0_M{1,r}, allData{1,t}.correlated_ROIs_dff0_M{1,c});
                    allData{1,t}.correlated_ROIs_dff0_MN_int{r,c} = intersect(allData{1,t}.correlated_ROIs_dff0_N{1,r}, allData{1,t}.correlated_ROIs_dff0_N{1,c});
                end
            end
        end
   end
end

%% Determine the ROIs that are the same across all cells with significant
% overlap with a given ROI
% significan overlap = at least 1/6*total num ROIS in correlated_ROIs_alldff0_int

%% Multisensory
clear('lens', 'int', 'roi_count')
for t = 1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_R_MS')
        if length(allData{1,t}.resp_ROIs) > 0
            roi_count = (1/6)*length(allData{1,t}.resp_ROIs);
            for row = 1:size(allData{1,t}.correlated_ROIs_dff0_MS_int,1)
                for ct = 1:size(allData{1,t}.correlated_ROIs_dff0_MS_int,2)
                    lens(ct) = length(allData{1,t}.correlated_ROIs_dff0_MS_int{ct,row});
                end
                first_roi = find((lens > roi_count), 1)
                if isempty(first_roi) 
                    continue
                else
                    int = allData{1,t}.correlated_ROIs_dff0_MS_int{row, first_roi}

                    for col = first_roi:size(allData{1,t}.correlated_ROIs_dff0_MS_int,2)
                        if lens(col) > roi_count
                            int = intersect(int, allData{1,t}.correlated_ROIs_dff0_MS_int{row, col});
                        else
                            continue
                        end
                    end
                allData{1,t}.correlated_ROIs_dff0_MS_common{row} = int;
                end
                clear('lens', 'int')
            end
            clear('roi_count')
        end
    end
end

% Index ROI numbers to the actual ROIs 
for t = 1:length(allData)
   % if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common')
            roi_list = allData{1,t}.resp_ROIs;
            for i = 1:length(allData{1,t}.correlated_ROIs_dff0_MS_common)
                if ~isempty(allData{1,t}.correlated_ROIs_dff0_MS_common(i))
                    allData{1,t}.correlated_ROIs_dff0_MS_common_AROI{i} = roi_list(allData{1,t}.correlated_ROIs_dff0_MS_common{i});
                end
            end
        end
        clear('roi_list')
    %end
end

%% Mechanosensory
% Multisensory
clear('lens', 'int', 'roi_count')
for t = 1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_R_M')
        if length(allData{1,t}.resp_ROIs) > 0
            roi_count = (1/6)*length(allData{1,t}.resp_ROIs);
            for row = 1:size(allData{1,t}.correlated_ROIs_dff0_M_int,1)
                for ct = 1:size(allData{1,t}.correlated_ROIs_dff0_M_int,2)
                    lens(ct) = length(allData{1,t}.correlated_ROIs_dff0_M_int{ct,row});
                end
                first_roi = find((lens > roi_count), 1)
                if isempty(first_roi) 
                    continue
                else
                    int = allData{1,t}.correlated_ROIs_dff0_M_int{row, first_roi}

                    for col = first_roi:size(allData{1,t}.correlated_ROIs_dff0_M_int,2)
                        if lens(col) > roi_count
                            int = intersect(int, allData{1,t}.correlated_ROIs_dff0_M_int{row, col});
                        else
                            continue
                        end
                    end
                allData{1,t}.correlated_ROIs_dff0_M_common{row} = int;
                end
                clear('lens', 'int')
            end
            clear('roi_count')
        end
    end
end

% Index ROI numbers to the actual ROIs 
for t = 1:length(allData)
   % if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_M_common')
            roi_list = allData{1,t}.resp_ROIs;
            for i = 1:length(allData{1,t}.correlated_ROIs_dff0_M_common)
                if ~isempty(allData{1,t}.correlated_ROIs_dff0_M_common(i))
                    allData{1,t}.correlated_ROIs_dff0_M_common_AROI{i} = roi_list(allData{1,t}.correlated_ROIs_dff0_M_common{i});
                end
            end
        end
        clear('roi_list')
    %end
end

%% Visual 
clear('lens', 'int', 'roi_count')
for t = 1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_R_V')
        if length(allData{1,t}.resp_ROIs) > 0
            roi_count = (1/6)*length(allData{1,t}.resp_ROIs);
            for row = 1:size(allData{1,t}.correlated_ROIs_dff0_V_int,1)
                for ct = 1:size(allData{1,t}.correlated_ROIs_dff0_V_int,2)
                    lens(ct) = length(allData{1,t}.correlated_ROIs_dff0_V_int{ct,row});
                end
                first_roi = find((lens > roi_count), 1)
                if isempty(first_roi) 
                    continue
                else
                    int = allData{1,t}.correlated_ROIs_dff0_V_int{row, first_roi}

                    for col = first_roi:size(allData{1,t}.correlated_ROIs_dff0_V_int,2)
                        if lens(col) > roi_count
                            int = intersect(int, allData{1,t}.correlated_ROIs_dff0_V_int{row, col});
                        else
                            continue
                        end
                    end
                allData{1,t}.correlated_ROIs_dff0_V_common{row} = int;
                end
                clear('lens', 'int')
            end
            clear('roi_count', 'lens')
        end
    end
end

% Index ROI numbers to the actual ROIs 
for t = 1:length(allData)
   % if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_V_common')
            roi_list = allData{1,t}.resp_ROIs;
            for i = 1:length(allData{1,t}.correlated_ROIs_dff0_V_common)
                if ~isempty(allData{1,t}.correlated_ROIs_dff0_V_common(i))
                    allData{1,t}.correlated_ROIs_dff0_V_common_AROI{i} = roi_list(allData{1,t}.correlated_ROIs_dff0_V_common{i});
                end
            end
        end
        clear('roi_list')
    %end
end

%% no stim 
clear('lens', 'int', 'roi_count')
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_N_int')
        if length(allData{1,t}.resp_ROIs) > 0
            roi_count = (1/6)*length(allData{1,t}.resp_ROIs);
            for row = 1:size(allData{1,t}.correlated_ROIs_dff0_N_int,1)
                for ct = 1:size(allData{1,t}.correlated_ROIs_dff0_N_int,2)
                    lens(ct) = length(allData{1,t}.correlated_ROIs_dff0_N_int{ct,row});
                end
                first_roi = find((lens > roi_count), 1)
                if isempty(first_roi) 
                    continue
                else
                    int = allData{1,t}.correlated_ROIs_dff0_N_int{row, first_roi}

                    for col = first_roi:size(allData{1,t}.correlated_ROIs_dff0_N_int,2)
                        if lens(col) > roi_count
                            int = intersect(int, allData{1,t}.correlated_ROIs_dff0_N_int{row, col});
                        else
                            continue
                        end
                    end
                allData{1,t}.correlated_ROIs_dff0_N_common{row} = int;
                end
                clear('lens', 'int')
            end
            clear('roi_count', 'lens')
        end
    end
end

% Index ROI numbers to the actual ROIs 
for t = 1:length(allData)
   % if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_N_common')
            roi_list = allData{1,t}.resp_ROIs;
            for i = 1:length(allData{1,t}.correlated_ROIs_dff0_N_common)
                if ~isempty(allData{1,t}.correlated_ROIs_dff0_N_common(i))
                    allData{1,t}.correlated_ROIs_dff0_N_common_AROI{i} = roi_list(allData{1,t}.correlated_ROIs_dff0_N_common{i});
                end
            end
        end
        clear('roi_list')
    %end
end

%% all trials 
clear('lens', 'int', 'roi_count')
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_all_int')
        if length(allData{1,t}.resp_ROIs) > 0
            roi_count = (1/6)*length(allData{1,t}.resp_ROIs);
            for row = 1:size(allData{1,t}.correlated_ROIs_dff0_all_int,1)
                for ct = 1:size(allData{1,t}.correlated_ROIs_dff0_all_int,2)
                    lens(ct) = length(allData{1,t}.correlated_ROIs_dff0_all_int{ct,row});
                end
                first_roi = find((lens > roi_count), 1)
                if isempty(first_roi) 
                    continue
                else
                    int = allData{1,t}.correlated_ROIs_dff0_all_int{row, first_roi}

                    for col = first_roi:size(allData{1,t}.correlated_ROIs_dff0_all_int,2)
                        if lens(col) > roi_count
                            int = intersect(int, allData{1,t}.correlated_ROIs_dff0_all_int{row, col});
                        else
                            continue
                        end
                    end
                allData{1,t}.correlated_ROIs_dff0_all_common{row} = int;
                end
                clear('lens', 'int')
            end
            clear('roi_count', 'lens')
        end
    end
end

% Index ROI numbers to the actual ROIs 
for t = 1:length(allData)
   % if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_all_common')
            roi_list = allData{1,t}.resp_ROIs;
            for i = 1:length(allData{1,t}.correlated_ROIs_dff0_all_common)
                if ~isempty(allData{1,t}.correlated_ROIs_dff0_all_common(i))
                    allData{1,t}.correlated_ROIs_dff0_all_common_AROI{i} = roi_list(allData{1,t}.correlated_ROIs_dff0_all_common{i});
                end
            end
        end
        clear('roi_list')
    %end
end

%%%%%%% Use highcorr_vs_nothighcorr code to create
%%%%%%% allData{1,t}.uniquehighcorrROI

%%%%%%% Use sfn2017_figures to generate tectum shaped scatter plots by
%%%%%%% tadpole of highcorr ROIs







%% Make plots of each xcorr ([maxR, lag of maxR, 0 lag R] for [all MS, V, M, N])

% Sort by average maxR in MS case, then use that order for everything
for t = 1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_maxR_sq_MS')
    avg = mean(allData{1,t}.respROIdff0_maxR_sq_MS)
    [B,I] = sort(avg, 'descend')
    allData{1,t}.respROI_maxR_MS_I = I;
    [ allData{1,t}.maxR_all_sorted ] = sort_ROIs( I, allData{1,t}.respROIdff0_maxR_sq_all )
    [ allData{1,t}.maxR_MS_sorted ] = sort_ROIs( I, allData{1,t}.respROIdff0_maxR_sq_MS )
    [ allData{1,t}.maxR_M_sorted ] = sort_ROIs( I, allData{1,t}.respROIdff0_maxR_sq_MS )
    [ allData{1,t}.maxR_V_sorted ] = sort_ROIs( I, allData{1,t}.respROIdff0_maxR_sq_V )
    [ allData{1,t}.maxR_N_sorted ] = sort_ROIs( I, allData{1,t}.respROIdff0_maxR_sq_N )
    [ allData{1,t}.lagmaxR_all_sorted ] = sort_ROIs( I, allData{1,t}.respROIdff0_maxRlag_sq_all )
    [ allData{1,t}.lagmaxR_MS_sorted ] = sort_ROIs( I, allData{1,t}.respROIdff0_maxRlag_sq_MS )
    [ allData{1,t}.lagmaxR_M_sorted ] = sort_ROIs( I, allData{1,t}.respROIdff0_maxRlag_sq_M )
    [ allData{1,t}.lagmaxR_V_sorted ] = sort_ROIs( I, allData{1,t}.respROIdff0_maxRlag_sq_V )
    [ allData{1,t}.lagmaxR_N_sorted ] = sort_ROIs( I, allData{1,t}.respROIdff0_maxRlag_sq_N )
    [ allData{1,t}.lag0R_all_sorted ] = sort_ROIs( I, allData{1,t}.respROIdff0_0lagR_sq_all )
    [ allData{1,t}.lag0R_MS_sorted ] = sort_ROIs( I, allData{1,t}.respROIdff0_0lagR_sq_MS )
    [ allData{1,t}.lag0R_M_sorted ] = sort_ROIs( I, allData{1,t}.respROIdff0_0lagR_sq_MS )
    [ allData{1,t}.lag0R_V_sorted ] = sort_ROIs( I, allData{1,t}.respROIdff0_0lagR_sq_V )
    [ allData{1,t}.lag0R_N_sorted ] = sort_ROIs( I, allData{1,t}.respROIdff0_0lagR_sq_N )
    end
end

% Plot xcorr maxR for all, MS, M, V, N
for t = 1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_maxR_sq_MS')
        plot_xcorr(allData{1,t}.maxR_all_sorted, 'hot')
        title(sprintf('Tad%d(t=%d) sorted maxR all', allData{1,t}.expnum, t))
        fig_filename = sprintf('Tad%d(t=%d) sorted maxR all', allData{1,t}.expnum, t)
        saveas(gcf, fig_filename, 'epsc2')
        close;
        plot_xcorr(allData{1,t}.maxR_MS_sorted, 'hot')
        title(sprintf('Tad%d(t=%d) sorted maxR MS', allData{1,t}.expnum, t))
        fig_filename = sprintf('Tad%d(t=%d) sorted maxR MS', allData{1,t}.expnum, t)
        saveas(gcf, fig_filename, 'epsc2')
        close;
        plot_xcorr(allData{1,t}.maxR_M_sorted, 'hot')
        title(sprintf('Tad%d(t=%d) sorted maxR M', allData{1,t}.expnum, t))
        fig_filename = sprintf('Tad%d(t=%d) sorted maxR M', allData{1,t}.expnum, t)
        saveas(gcf, fig_filename, 'epsc2')
        close;
        plot_xcorr(allData{1,t}.maxR_V_sorted, 'hot')
        title(sprintf('Tad%d(t=%d) sorted maxR V', allData{1,t}.expnum, t))
        fig_filename = sprintf('Tad%d(t=%d) sorted maxR V', allData{1,t}.expnum, t)
        saveas(gcf, fig_filename, 'epsc2')
        close;
        plot_xcorr(allData{1,t}.maxR_N_sorted, 'hot')
        title(sprintf('Tad%d(t=%d) sorted maxR N', allData{1,t}.expnum, t))
        fig_filename = sprintf('Tad%d(t=%d) sorted maxR N', allData{1,t}.expnum, t)
        saveas(gcf, fig_filename, 'epsc2')
        close;
    end    
end

% Plot xcorr lag of maxR for all, MS, M, V, N
for t = 1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_maxR_sq_MS')
        plot_xcorr(allData{1,t}.lagmaxR_all_sorted, 'jet')
        title(sprintf('Tad%d(t=%d) sorted lagmaxR all', allData{1,t}.expnum, t))
        fig_filename = sprintf('Tad%d(t=%d) sorted lagmaxR all', allData{1,t}.expnum, t)
        saveas(gcf, fig_filename, 'epsc2')
        close;
        plot_xcorr(allData{1,t}.lagmaxR_MS_sorted, 'jet')
        title(sprintf('Tad%d(t=%d) sorted lagmaxR MS', allData{1,t}.expnum, t))
        fig_filename = sprintf('Tad%d(t=%d) sorted lagmaxR MS', allData{1,t}.expnum, t)
        saveas(gcf, fig_filename, 'epsc2')
        close;
        plot_xcorr(allData{1,t}.lagmaxR_M_sorted, 'jet')
        title(sprintf('Tad%d(t=%d) sorted lagmaxR M', allData{1,t}.expnum, t))
        fig_filename = sprintf('Tad%d(t=%d) sorted lagmaxR M', allData{1,t}.expnum, t)
        saveas(gcf, fig_filename, 'epsc2')
        close;
        plot_xcorr(allData{1,t}.lagmaxR_V_sorted, 'jet')
        title(sprintf('Tad%d(t=%d) sorted lagmaxR V', allData{1,t}.expnum, t))
        fig_filename = sprintf('Tad%d(t=%d) sorted lagmaxR V', allData{1,t}.expnum, t)
        saveas(gcf, fig_filename, 'epsc2')
        close;
        plot_xcorr(allData{1,t}.lagmaxR_N_sorted, 'jet')
        title(sprintf('Tad%d(t=%d) sorted lagmaxR N', allData{1,t}.expnum, t))
        fig_filename = sprintf('Tad%d(t=%d) sorted lagmaxR N', allData{1,t}.expnum, t)
        saveas(gcf, fig_filename, 'epsc2')
        close;
    end    
end

% Plot xcorr R at 0lag for all, MS, M, V, N
for t = 1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_maxR_sq_MS')
        plot_xcorr(allData{1,t}.lag0R_all_sorted, 'hot')
        title(sprintf('Tad%d(t=%d) sorted lag0R all', allData{1,t}.expnum, t))
        fig_filename = sprintf('Tad%d(t=%d) sorted lag0R all', allData{1,t}.expnum, t);
        saveas(gcf, fig_filename, 'epsc2')
        close;
        plot_xcorr(allData{1,t}.lag0R_MS_sorted, 'hot')
        title(sprintf('Tad%d(t=%d) sorted lag0R MS', allData{1,t}.expnum, t))
        fig_filename = sprintf('Tad%d(t=%d) sorted lag0R MS', allData{1,t}.expnum, t)
        saveas(gcf, fig_filename, 'epsc2')
        close;
        plot_xcorr(allData{1,t}.lag0R_M_sorted, 'hot')
        title(sprintf('Tad%d(t=%d) sorted lag0R M', allData{1,t}.expnum, t))
        fig_filename = sprintf('Tad%d(t=%d) sorted lag0R M', allData{1,t}.expnum, t)
        saveas(gcf, fig_filename, 'epsc2')
        close;
        plot_xcorr(allData{1,t}.lag0R_V_sorted, 'hot')
        title(sprintf('Tad%d(t=%d) sorted lag0R V', allData{1,t}.expnum, t))
        fig_filename = sprintf('Tad%d(t=%d) sorted lag0R V', allData{1,t}.expnum, t)
        saveas(gcf, fig_filename, 'epsc2')
        close;
        plot_xcorr(allData{1,t}.lag0R_N_sorted, 'hot')
        title(sprintf('Tad%d(t=%d) sorted lag0R N', allData{1,t}.expnum, t))
        fig_filename = sprintf('Tad%d(t=%d) sorted lag0R N', allData{1,t}.expnum, t)
        saveas(gcf, fig_filename, 'epsc2')
        close;
    end    
end


%% Fix elimination criteria to be less strigent

% histogram of peaks by trial
for t = 1:length(allData)
    figure;
    hist(allData{1,t}.peak_bytrial)
    fig_filename = sprintf('hist peakbytrial tad%d(t=%d)', allData{1,t}.expnum, t)
    saveas(gcf, fig_filename, 'png')
end

% range is weird - 0-5 in some, up to 5000 in others.

% re-calculate peak_bytrial using df/f0 (this should be smoothed)

for t = 1:length(allData)
    start = floor((length(allData{1,t}.df_f0{1,1})/7)*2)
    [allData{1,t}.peak_bytrial2, allData{1,t}.peak_loc2] = calc_peak2(allData{1,t}.df_f0, start, 2);
end

% rerun histogram
for t = 1:length(allData)
    figure;
    hist(allData{1,t}.peak_bytrial2)
    fig_filename = sprintf('hist peakbytrial2 tad%d(t=%d)', allData{1,t}.expnum, t)
    saveas(gcf, fig_filename, 'png')
end
%these have smaller ranges - more like expected. Go back and rerun
%everything with this peak value. 



%% Characterize highcorr ROIs
















%% Characterize xcorr generally

%using old code from correlations_v2_xcorr

% distribution of R values for each stim type (1 figure per tadpole)
for t = 1:length(allData)
    %if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if length(allData{1,t}.resp_ROIs) > 0
        figure;
        subplot(2,2,1)
        hist(allData{1,t}.respROIdff0_maxR_sq_MS,40)
        title('Multi')
        subplot(2,2,2)
        hist(allData{1,t}.respROIdff0_maxR_sq_V,40)
        title('Vis')
        subplot(2,2,3)
        hist(allData{1,t}.respROIdff0_maxR_sq_M,40)
        title('Mech')
        subplot(2,2,4)
        hist(allData{1,t}.respROIdff0_maxR_sq_N,40)
        title('None')
        suptitle(sprintf('tad %d(t=%d) responding cells hist of xcorr R vals by stimtype', allData{1,t}.expnum, t))
        fig_filename = sprintf('tad %d(t=%d) responding cells hist of xcorr Rvals by stimtype', allData{1,t}.expnum, t)
        saveas(gcf,fig_filename,'png');
        close;
        %end
    end
end

% distribution of lag of maxR for each stimtype (1 figure per tadpole
for t = 1:length(allData)
    %if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if length(allData{1,t}.resp_ROIs) > 0
        lim_min = min(allData{1,t}.respROIdff0_lag_all{1,1});
        lim_max = max(allData{1,t}.respROIdff0_lag_all{1,1});
        figure;
        subplot(2,2,1)
        hist(allData{1,t}.respROIdff0_maxRlag_sq_MS,40)
        xlim([lim_min, lim_max])
        title('Multi')
        subplot(2,2,2)
        hist(allData{1,t}.respROIdff0_maxRlag_sq_V,40)
        xlim([lim_min, lim_max])
        title('Vis')
        subplot(2,2,3)
        hist(allData{1,t}.respROIdff0_maxRlag_sq_M,40)
        xlim([lim_min, lim_max])
        title('Mech')
        subplot(2,2,4)
        hist(allData{1,t}.respROIdff0_maxRlag_sq_N,40)
        xlim([lim_min, lim_max])
        title('None')
        suptitle(sprintf('tad %d(t=%d) responding cells hist of xcorr lag of maxR by stimtype', allData{1,t}.expnum, t))
        fig_filename = sprintf('tad %d(t=%d) responding cells hist of xcorr lag of maxR by stimtype', allData{1,t}.expnum, t)
        saveas(gcf,fig_filename,'png');
        close;
        end
    %end
end

% distribution of R values for 0 lag for each stim type (1 figure per tadpole)
for t = 1:length(allData)
    %if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if length(allData{1,t}.resp_ROIs) > 0
        figure;
        subplot(2,2,1)
        hist(allData{1,t}.respROIdff0_0lagR_sq_MS,40)
        title('Multi')
        subplot(2,2,2)
        hist(allData{1,t}.respROIdff0_0lagR_sq_V,40)
        title('Vis')
        subplot(2,2,3)
        hist(allData{1,t}.respROIdff0_0lagR_sq_M,40)
        title('Mech')
        subplot(2,2,4)
        hist(allData{1,t}.respROIdff0_0lagR_sq_N,40)
        title('None')
        suptitle(sprintf('tad %d(t=%d) responding cells hist of xcorr 0lag R vals by stimtype', allData{1,t}.expnum, t))
        fig_filename = sprintf('tad %d(t=%d) responding cells hist of xcorr 0lag Rvals by stimtype', allData{1,t}.expnum, t)
        saveas(gcf,fig_filename,'png');
        close;
        %end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%















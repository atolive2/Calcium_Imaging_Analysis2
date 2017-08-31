%% Since the data doesn't cluster neatly, let's look at empirical probability

%started with tadcluster_data_20170814_clean.mat which is all the basic ROI
%info

%% First collect all relevant parameters into 1 matrix (same as cluster_analysis_v2)

% Carlos suggests: avg peak, jitter, mean squared error, hindbrain bias, peaks by modality
% overall strategy: calculate for all cells, then extract just the cells I
% want to cluster at the end.

%% peaks by modality - mean and stddev of each stim type
% avg_peak is in tadcluster{1,t}.peak_avg

for t = 1:length(tadcluster)
    tadcluster{1,t}.stimmask = get_stimmask(tadcluster{1,t}.stimorder);
end
for t = 1:length(tadcluster)
    tadcluster{1,t}.peak_stddev_bymod = std_by_stimtype(tadcluster{1,t}.peak_bytrial, tadcluster{1,t}.stimmask);
end

%% jitter = variation in onset time (separate by modality)
% copy code from analyze_onset_time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%add stim onset time to tadcluster%%%%
% exps 1-9 have onset at 0.5s 
% exps 10+ have onset at 2s
for t = 1:length(tadcluster)
    if tadcluster{1,t}.expnum <=9
        tadcluster{1,t}.stim_onset = 0.5;
    elseif tadcluster{1,t}.expnum >9
        tadcluster{1,t}.stim_onset = 2;
    else
        fprintf('error exp %d, t')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get time to 50% of the peak (= response onset time)
% same as used in analyze_onset_time
for t = 1:length(tadcluster)
    for r = 1:size(tadcluster{1,t}.df_f0,1) %over each ROI
        for tt = 1:size(tadcluster{1,t}.df_f0,2) %over each trial
            half_peak = tadcluster{1,t}.peak_bytrial(r, tt) / 2;
            tmp_onset = find((tadcluster{1,t}.df_f0{r, tt} > half_peak) , 3);
            % get trial length (this takes care of exps with not 160 frames)
            trial_length = length(tadcluster{1,t}.df_f0{r, tt});
            %find the proper onset time (e.g. eliminate errors of 1 or 2)
            if length(tmp_onset) < 3 %this sets onset time to 0 for trials with no response
                onset_time(r, tt) = 0;
            else
                if tmp_onset(1) == (1 || 2)
                    if tmp_onset(2) == 2 
                        onset_time(r, tt) = tmp_onset(3) / (trial_length/7);
                    else
                        onset_time(r, tt) = tmp_onset(2) / (trial_length/7);
                    end
                else
                    onset_time(r, tt) = tmp_onset(1) / (trial_length/7);
                end
            end
        end
    end
    % put onset_time into tadcluster
    tadcluster{1,t}.onset_time = onset_time;    
end

% Get mean and sd of onset_time by ROI (this includes 0s - no response trials)
for t = 1:length(tadcluster)
    tadcluster{1,t}.onset_time_mean = mean_by_stimtype2(tadcluster{1,t}.onset_time, tadcluster{1,t}.stimmask);
    tadcluster{1,t}.onset_time_stdev = std_by_stimtype(tadcluster{1,t}.onset_time, tadcluster{1,t}.stimmask);
end

%% multisensory enhancement (recalculate by num responses)

% use tadcluster{1,t}.MSEnh_peak 

% number of responses
for  t = 1:length(tadcluster)
    for r = 1:size(tadcluster{1,t}.boolean_response, 1)
        multi = sum(tadcluster{1,t}.boolean_response(r, tadcluster{1,t}.stimmask(:,1))) / sum(tadcluster{1,t}.stimmask(:,1));
        vis = sum(tadcluster{1,t}.boolean_response(r, tadcluster{1,t}.stimmask(:,2))) / sum(tadcluster{1,t}.stimmask(:,2));
        mech = sum(tadcluster{1,t}.boolean_response(r, tadcluster{1,t}.stimmask(:,3))) / sum(tadcluster{1,t}.stimmask(:,3));
        if (vis + mech) == 0 
            if multi > 0
                tadcluster{1,t}.MSEnh_numresponses(r) = 1;
            else
                tadcluster{1,t}.MSEnh_numresponses(r) = 0; 
            end
        else
            tadcluster{1,t}.MSEnh_numresponses(r) = multi / (vis + mech);
        end
    end
end

%% unisensory bias = ratio of vis responses to mech responses

% get number of responses to vis and mech, divide vis by mech
for  t = 1:length(tadcluster)
    for r = 1:size(tadcluster{1,t}.boolean_response, 1)
        vis = sum(tadcluster{1,t}.boolean_response(r, tadcluster{1,t}.stimmask(:,2))) / sum(tadcluster{1,t}.stimmask(:,2));
        mech = sum(tadcluster{1,t}.boolean_response(r, tadcluster{1,t}.stimmask(:,3))) / sum(tadcluster{1,t}.stimmask(:,3));
        if (vis + mech) == 0 %doesn't respond to either
            tadcluster{1,t}.unibias_numresponses(r) = 0;
        else
            tadcluster{1,t}.unibias_numresponses(r) = vis/(vis+mech);
        end
    end
end

%% What is the empirical cumulative distribution function for each modality?

% calculate the response percentage for all ROIs for each modality
for  t = 1:length(tadcluster)
    for r = 1:size(tadcluster{1,t}.boolean_response, 1)
        tadcluster{1,t}.Prop_resp_MS(r) = sum(tadcluster{1,t}.boolean_response(r, tadcluster{1,t}.stimmask(:,1))) / sum(tadcluster{1,t}.stimmask(:,1));
        tadcluster{1,t}.Prop_resp_V(r) = sum(tadcluster{1,t}.boolean_response(r, tadcluster{1,t}.stimmask(:,2))) / sum(tadcluster{1,t}.stimmask(:,2));
        tadcluster{1,t}.Prop_resp_M(r) = sum(tadcluster{1,t}.boolean_response(r, tadcluster{1,t}.stimmask(:,3))) / sum(tadcluster{1,t}.stimmask(:,3));
        tadcluster{1,t}.Prop_resp_N(r) = sum(tadcluster{1,t}.boolean_response(r, tadcluster{1,t}.stimmask(:,4))) / sum(tadcluster{1,t}.stimmask(:,4));
    end
end

% calculate the empirical cumulative distribution function
% https://www.mathworks.com/help/stats/ecdf.html
for t = 1:length(tadcluster)
    if length(tadcluster{1,t}.resp_ROIs) > 17
        [tadcluster{1,t}.F_MS, tadcluster{1,t}.X_MS] = ecdf(tadcluster{1,t}.Prop_resp_MS(tadcluster{1,t}.resp_ROIs));
        [tadcluster{1,t}.F_V, tadcluster{1,t}.X_V] = ecdf(tadcluster{1,t}.Prop_resp_V(tadcluster{1,t}.resp_ROIs));
        [tadcluster{1,t}.F_M, tadcluster{1,t}.X_M] = ecdf(tadcluster{1,t}.Prop_resp_M(tadcluster{1,t}.resp_ROIs));
        [tadcluster{1,t}.F_N, tadcluster{1,t}.X_N] = ecdf(tadcluster{1,t}.Prop_resp_N(tadcluster{1,t}.resp_ROIs));
    end
end

% plot the ECDFs by tad
for t = 1:length(tadcluster)
    if length(tadcluster{1,t}.resp_ROIs) > 17
        figure;
        hold on
        plot(tadcluster{1,t}.F_MS, tadcluster{1,t}.X_MS, 'm')
        plot(tadcluster{1,t}.F_V, tadcluster{1,t}.X_V, 'r')
        plot(tadcluster{1,t}.F_M, tadcluster{1,t}.X_M, 'b')
        plot(tadcluster{1,t}.F_N, tadcluster{1,t}.X_N, 'k')
        hold off
        title(sprintf('tad %d ECDF of responses', t))
        xlabel('proportion responses')
        ylabel('ROI count')
        fig_filename = sprintf('tad %d ECDF of responding ROIs', t)
        saveas(gcf, fig_filename, 'png')
        close;
    end
end

% combine and plot ECDF of all respROIs 
% combine all Prop_resp of resp_ROIs into 1 vector 
resp_ROI_prop_respMS = [];
resp_ROI_prop_respV = [];
resp_ROI_prop_respM = [];
resp_ROI_prop_respN = [];
for t = 1:length(tadcluster)
    resp_ROI_prop_respMS = [resp_ROI_prop_respMS tadcluster{1,t}.Prop_resp_MS(tadcluster{1,t}.resp_ROIs)];
    resp_ROI_prop_respV = [resp_ROI_prop_respV tadcluster{1,t}.Prop_resp_V(tadcluster{1,t}.resp_ROIs)];
    resp_ROI_prop_respM = [resp_ROI_prop_respM tadcluster{1,t}.Prop_resp_M(tadcluster{1,t}.resp_ROIs)];
    resp_ROI_prop_respN = [resp_ROI_prop_respN tadcluster{1,t}.Prop_resp_N(tadcluster{1,t}.resp_ROIs)];
end

% calcualte ECDF
[F_resp_ROI_prop_respMS X_resp_ROI_prop_respMS, flo_MS, flu_MS] = ecdf(resp_ROI_prop_respMS);
[F_resp_ROI_prop_respV X_resp_ROI_prop_respV, flo_V, flu_V] = ecdf(resp_ROI_prop_respV);
[F_resp_ROI_prop_respM X_resp_ROI_prop_respM, flo_M, flu_M] = ecdf(resp_ROI_prop_respM);
[F_resp_ROI_prop_respN X_resp_ROI_prop_respN, flo_N, flu_N] = ecdf(resp_ROI_prop_respN);

%create figure
figure;
hold on
plot(F_resp_ROI_prop_respMS, X_resp_ROI_prop_respMS, 'm', 'LineWidth', 1)
plot(F_resp_ROI_prop_respV, X_resp_ROI_prop_respV, 'r', 'LineWidth', 1)
plot(F_resp_ROI_prop_respM, X_resp_ROI_prop_respM, 'b', 'LineWidth', 1)
plot(F_resp_ROI_prop_respN, X_resp_ROI_prop_respN, 'k', 'LineWidth', 1)
%hold off
xlabel('proportion responses')
ylabel('ROI count')
title('ECDF of all responding ROIs')
annotation('textbox', 'Position', [0.2 0.75 .1 .1], 'String', ['Multi'], 'Color', 'm', 'LineStyle', 'none' );
annotation('textbox', 'Position', [0.2 0.7 .1 .1], 'String', ['Vis'], 'Color', 'b', 'LineStyle', 'none' );
annotation('textbox', 'Position', [0.2 0.65 .1 .1], 'String', ['Mech'], 'Color', 'r', 'LineStyle', 'none' );
annotation('textbox', 'Position', [0.2 0.6 .1 .1], 'String', ['No stim'], 'Color', 'k', 'LineStyle', 'none' );
%add confidence intervals
%ciplot(flo_MS, flu_MS, F_resp_ROI_prop_respMS, 'm')

% stat test to see if they're different
all_mods = [resp_ROI_prop_respMS', resp_ROI_prop_respV', resp_ROI_prop_respM', resp_ROI_prop_respN'];
size(all_mods)
[p,tbl,stats] = kruskalwallis(all_mods);
c = multcompare(stats)
% yes they are! 


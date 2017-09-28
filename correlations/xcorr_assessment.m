%% Correlation Assessment of all stage 49 tads (exps 1-40)

% data is tadcluster from both exps 1-30 and 31, 34, 40
load('tadcluster_xcorr_31,34,40')
%rename tadcluster to tadcluster1
load('tadcluster_analysis_xcorr_20170818.mat') 

% add 3 exps to the rest
tadcluster = [tadcluster, tadcluster1];

%% Does general/average correlation change by stimtype?

% Collect maxR for each cell with each other cell into 1 vector
% multi
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxR_sq_MS')
        tmp_data = [];
        for i = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_MS,1)
            for j = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_MS,2)
                tmp_data = [tmp_data, tadcluster{1,t}.respROIdff0_maxR_sq_MS(i, (j+1):end)];
            end
        end                    
        tadcluster{1,t}.respROIdff0_maxR_MS_v = tmp_data;
    end
end

%vis
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxR_sq_V')
        tmp_data = [];
        for i = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_V,1)
            for j = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_V,2)
                tmp_data = [tmp_data, tadcluster{1,t}.respROIdff0_maxR_sq_V(i, (j+1):end)];
            end
        end                    
        tadcluster{1,t}.respROIdff0_maxR_V_v = tmp_data;
    end
end

%mech
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxR_sq_M')
        tmp_data = [];
        for i = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_M,1)
            for j = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_M,2)
                tmp_data = [tmp_data, tadcluster{1,t}.respROIdff0_maxR_sq_M(i, (j+1):end)];
            end
        end                    
        tadcluster{1,t}.respROIdff0_maxR_M_v = tmp_data;
    end
end

%no stim
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxR_sq_N')
        tmp_data = [];
        for i = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_N,1)
            for j = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_N,2)
                tmp_data = [tmp_data, tadcluster{1,t}.respROIdff0_maxR_sq_N(i, (j+1):end)];
            end
        end                    
        tadcluster{1,t}.respROIdff0_maxR_N_v = tmp_data;
    end
end

% Get cumulative distribution for each tad
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxR_MS_v')
    [tadcluster{1,t}.ecdf_maxR_MS_F, tadcluster{1,t}.ecdf_maxR_MS_X] = ecdf(tadcluster{1,t}.respROIdff0_maxR_MS_v);
    [tadcluster{1,t}.ecdf_maxR_V_F, tadcluster{1,t}.ecdf_maxR_V_X] = ecdf(tadcluster{1,t}.respROIdff0_maxR_V_v);
    [tadcluster{1,t}.ecdf_maxR_M_F, tadcluster{1,t}.ecdf_maxR_M_X] = ecdf(tadcluster{1,t}.respROIdff0_maxR_M_v);
    [tadcluster{1,t}.ecdf_maxR_N_F, tadcluster{1,t}.ecdf_maxR_N_X] = ecdf(tadcluster{1,t}.respROIdff0_maxR_N_v);
    end
 end

% plot ECDF for each tad
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxR_MS_v')
        figure;
        hold on
        plot(tadcluster{1,t}.ecdf_maxR_MS_X, tadcluster{1,t}.ecdf_maxR_MS_F, 'm')
        plot(tadcluster{1,t}.ecdf_maxR_V_X, tadcluster{1,t}.ecdf_maxR_V_F, 'r')
        plot(tadcluster{1,t}.ecdf_maxR_M_X, tadcluster{1,t}.ecdf_maxR_M_F, 'b')
        plot(tadcluster{1,t}.ecdf_maxR_N_X, tadcluster{1,t}.ecdf_maxR_N_F, 'k')
        hold off
        title(sprintf('tad %d ECDF of maxR', t))
        xlabel('maxR')
        ylabel('ROI proportion')
        annotation('textbox', 'Position', [0.2 0.75 .1 .1], 'String', ['Multi'], 'Color', 'm', 'LineStyle', 'none' );
        annotation('textbox', 'Position', [0.2 0.7 .1 .1], 'String', ['Vis'], 'Color', 'b', 'LineStyle', 'none' );
        annotation('textbox', 'Position', [0.2 0.65 .1 .1], 'String', ['Mech'], 'Color', 'r', 'LineStyle', 'none' );
        annotation('textbox', 'Position', [0.2 0.6 .1 .1], 'String', ['No stim'], 'Color', 'k', 'LineStyle', 'none' );

        fig_filename = sprintf('tad %d ECDF of maxR of responding ROIs', t)
        saveas(gcf, fig_filename, 'png')
        close;
    end
end

% Statistical difference in any tad?
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxR_MS_v')
        all_mods = [tadcluster{1,t}.respROIdff0_maxR_MS_v', tadcluster{1,t}.respROIdff0_maxR_V_v', tadcluster{1,t}.respROIdff0_maxR_M_v', tadcluster{1,t}.respROIdff0_maxR_N_v'];
        size(all_mods)
        [tadcluster{1,t}.maxR_p, tadcluster{1,t}.maxR_tbl, tadcluster{1,t}.maxR_stats] = kruskalwallis(all_mods);
        tadcluster{1,t}.maxR_c = multcompare(tadcluster{1,t}.maxR_stats)
    end
end
% collect P vals into 1 vector for easy reading
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'maxR_p')
        maxR_allP(t) = tadcluster{1,t}.maxR_p
    elseif ~isfield(tadcluster{1,t}, 'maxR_p')
        maxR_allP(t) = NaN 
    end
end

% if P is significant, collect the P vals of the multiple comparison into 1
% matrix for easy reading.
for t = 1:length(maxR_allP)
    if maxR_allP(t) < 0.05
        which_mod_diff(:,t) = tadcluster{1,t}.maxR_c(:,6)
    else %if maxR_allP(t) == NaN
        which_mod_diff(:,t) = NaN
    end
end
% how many experiments are significant (e.g. have values not NaN)
sig_tads = find(~isnan(which_mod_diff(1,:)))
num_sig_tads = length(sig_tads)
% 14 of 18 exps are significantly different. Exps 5, 7, 9, 13 are not. Exp
% 5 does not have stimtype = 1 (high multi). 13 has no respROIs. 

% Of the experiemnts with significant P val, how many tads are significant in each stim comparison?
for i = 1:size(which_mod_diff,1)
    num_tads_diff(i,1) = length( find((which_mod_diff(i,:) ~= 0) & (which_mod_diff(i,:) < 0.05)));
end

% bar plot the significant combos
labels = {'MS vs V', 'MS vs M', 'MS vs N', 'V vs M', 'V vs N', 'M vs N'}
bar(num_tads_diff) %, labels)
set(gca, 'XTickLabel',labels, 'XTick',1:numel(labels))
ylabel('Num Sig Pairs')
ylim([0 14])
title('Num of tads with P < 0.05 in mult comparison if P_{KW} < 0.05')

%% Combine across tadpoles to compare maxR vals across stim

%multi
tmp_data = [];
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxR_sq_MS')
        for i = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_MS,1)
            for j = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_MS,2)
                tmp_data = [tmp_data, tadcluster{1,t}.respROIdff0_maxR_sq_MS(i, (j+1):end)];
            end
        end                     
    end
end
maxR_allrespROI_MS = tmp_data;

% vis
tmp_data = [];
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxR_sq_V')
        for i = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_V,1)
            for j = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_V,2)
                tmp_data = [tmp_data, tadcluster{1,t}.respROIdff0_maxR_sq_V(i, (j+1):end)];
            end
        end                     
    end
end
maxR_allrespROI_V = tmp_data;

% mech
tmp_data = [];
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxR_sq_M')
        for i = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_M,1)
            for j = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_M,2)
                tmp_data = [tmp_data, tadcluster{1,t}.respROIdff0_maxR_sq_M(i, (j+1):end)];
            end
        end                     
    end
end
maxR_allrespROI_M = tmp_data;

%no stim
tmp_data = [];
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxR_sq_N')
        for i = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_N,1)
            for j = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_N,2)
                tmp_data = [tmp_data, tadcluster{1,t}.respROIdff0_maxR_sq_N(i, (j+1):end)];
            end
        end                     
    end
end
maxR_allrespROI_N = tmp_data;
clear('tmp_data')

% how many cells is this total?
roiCount = 0
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxR_sq_MS')
        roiCount = roiCount + length(tadcluster{1,t}.resp_ROIs);
    end
end

% make ECDF of all combined data
figure;
hold on
ecdf(maxR_allrespROI_MS)
ecdf(maxR_allrespROI_V)
ecdf(maxR_allrespROI_M)
ecdf(maxR_allrespROI_N)
title('ECDF of maxR all tads all respROIs')
xlabel('maxR')
ylabel('ROI proportion')
annotation('textbox', 'Position', [0.2 0.75 .1 .1], 'String', ['Multi'], 'Color', 'b', 'LineStyle', 'none' );
annotation('textbox', 'Position', [0.2 0.7 .1 .1], 'String', ['Vis'], 'Color', 'r', 'LineStyle', 'none' );
annotation('textbox', 'Position', [0.2 0.65 .1 .1], 'String', ['Mech'], 'Color', 'y', 'LineStyle', 'none' );
annotation('textbox', 'Position', [0.2 0.6 .1 .1], 'String', ['No stim'], 'Color', [0.5 0 0.5], 'LineStyle', 'none' );
annotation('textbox', 'Position', [0.2 0.55 .1 .1], 'String', 'N = 16, n = 474', 'Color', 'k', 'LineStyle', 'none');
saveas(gcf, 'ECDF of maxR all tads all respROIs', 'png')

% stat test all combined data
% Statistical difference in any tad?
all_mods = [maxR_allrespROI_MS', maxR_allrespROI_V', maxR_allrespROI_M', maxR_allrespROI_N'];
size(all_mods)
[maxR_p, maxR_tbl, maxR_stats] = kruskalwallis(all_mods)
maxR_c = multcompare(maxR_stats)

%% Repeat this analysis on lag times - specifically the lag time of maxR
% Collect lag time of maxR for each cell with each other cell into 1 vector
% multi
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxRlag_sq_MS')
        tmp_data = [];
        for i = 1:size(tadcluster{1,t}.respROIdff0_maxRlag_sq_MS,1)
            for j = 1:size(tadcluster{1,t}.respROIdff0_maxRlag_sq_MS,2)
                tmp_data = [tmp_data, tadcluster{1,t}.respROIdff0_maxRlag_sq_MS(i, (j+1):end)];
            end
        end                    
        tadcluster{1,t}.respROIdff0_maxRlag_MS_v = tmp_data;
    end
end

%vis
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxRlag_sq_V')
        tmp_data = [];
        for i = 1:size(tadcluster{1,t}.respROIdff0_maxRlag_sq_V,1)
            for j = 1:size(tadcluster{1,t}.respROIdff0_maxRlag_sq_V,2)
                tmp_data = [tmp_data, tadcluster{1,t}.respROIdff0_maxRlag_sq_V(i, (j+1):end)];
            end
        end                    
        tadcluster{1,t}.respROIdff0_maxRlag_V_v = tmp_data;
    end
end

%mech
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxRlag_sq_M')
        tmp_data = [];
        for i = 1:size(tadcluster{1,t}.respROIdff0_maxRlag_sq_M,1)
            for j = 1:size(tadcluster{1,t}.respROIdff0_maxRlag_sq_M,2)
                tmp_data = [tmp_data, tadcluster{1,t}.respROIdff0_maxRlag_sq_M(i, (j+1):end)];
            end
        end                    
        tadcluster{1,t}.respROIdff0_maxRlag_M_v = tmp_data;
    end
end

%no stim
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxRlag_sq_N')
        tmp_data = [];
        for i = 1:size(tadcluster{1,t}.respROIdff0_maxRlag_sq_N,1)
            for j = 1:size(tadcluster{1,t}.respROIdff0_maxRlag_sq_N,2)
                tmp_data = [tmp_data, tadcluster{1,t}.respROIdff0_maxRlag_sq_N(i, (j+1):end)];
            end
        end                    
        tadcluster{1,t}.respROIdff0_maxRlag_N_v = tmp_data;
    end
end

% Get cumulative distribution for each tad
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxR_MS_v')
    [tadcluster{1,t}.ecdf_maxRlag_MS_F, tadcluster{1,t}.ecdf_maxRlag_MS_X] = ecdf(tadcluster{1,t}.respROIdff0_maxRlag_MS_v);
    [tadcluster{1,t}.ecdf_maxRlag_V_F, tadcluster{1,t}.ecdf_maxRlag_V_X] = ecdf(tadcluster{1,t}.respROIdff0_maxRlag_V_v);
    [tadcluster{1,t}.ecdf_maxRlag_M_F, tadcluster{1,t}.ecdf_maxRlag_M_X] = ecdf(tadcluster{1,t}.respROIdff0_maxRlag_M_v);
    [tadcluster{1,t}.ecdf_maxRlag_N_F, tadcluster{1,t}.ecdf_maxRlag_N_X] = ecdf(tadcluster{1,t}.respROIdff0_maxRlag_N_v);
    end
 end

% plot ECDF for each tad
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxR_MS_v')
        figure;
        hold on
        plot(tadcluster{1,t}.ecdf_maxRlag_MS_X, tadcluster{1,t}.ecdf_maxRlag_MS_F, 'm')
        plot(tadcluster{1,t}.ecdf_maxRlag_V_X, tadcluster{1,t}.ecdf_maxRlag_V_F, 'r')
        plot(tadcluster{1,t}.ecdf_maxRlag_M_X, tadcluster{1,t}.ecdf_maxRlag_M_F, 'b')
        plot(tadcluster{1,t}.ecdf_maxRlag_N_X, tadcluster{1,t}.ecdf_maxRlag_N_F, 'k')
        hold off
        title(sprintf('tad %d ECDF of maxR lag time', t))
        xlabel('maxR')
        ylabel('ROI proportion')
        annotation('textbox', 'Position', [0.2 0.75 .1 .1], 'String', ['Multi'], 'Color', 'm', 'LineStyle', 'none' );
        annotation('textbox', 'Position', [0.2 0.7 .1 .1], 'String', ['Vis'], 'Color', 'b', 'LineStyle', 'none' );
        annotation('textbox', 'Position', [0.2 0.65 .1 .1], 'String', ['Mech'], 'Color', 'r', 'LineStyle', 'none' );
        annotation('textbox', 'Position', [0.2 0.6 .1 .1], 'String', ['No stim'], 'Color', 'k', 'LineStyle', 'none' );

        fig_filename = sprintf('tad %d ECDF of maxR lag time of responding ROIs', t)
        saveas(gcf, fig_filename, 'png')
        close;
    end
end

% Statistical difference in any tad?
all_mods = [];
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxRlag_MS_v')
        all_mods = [tadcluster{1,t}.respROIdff0_maxRlag_MS_v', tadcluster{1,t}.respROIdff0_maxRlag_V_v', tadcluster{1,t}.respROIdff0_maxRlag_M_v', tadcluster{1,t}.respROIdff0_maxRlag_N_v'];
        size(all_mods)
        [tadcluster{1,t}.maxRlag_p, tadcluster{1,t}.maxRlag_tbl, tadcluster{1,t}.maxRlag_stats] = kruskalwallis(all_mods);
        tadcluster{1,t}.maxRlag_c = multcompare(tadcluster{1,t}.maxRlag_stats)
    end
end
% collect P vals into 1 vector for easy reading
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'maxRlag_p')
        maxRlag_allP(t) = tadcluster{1,t}.maxRlag_p
    elseif ~isfield(tadcluster{1,t}, 'maxRlag_p')
        maxRlag_allP(t) = NaN 
    end
end

% if P is significant, collect the P vals of the multiple comparison into 1
% matrix for easy reading.
for t = 1:length(maxRlag_allP)
    if maxRlag_allP(t) < 0.05
        which_mod_diff_lag(:,t) = tadcluster{1,t}.maxRlag_c(:,6)
    else %if maxR_allP(t) == NaN
        which_mod_diff_lag(:,t) = NaN
    end
end
% how many experiments are significant (e.g. have values not NaN)
sig_tads_lag = find(~isnan(which_mod_diff_lag(1,:)))
num_sig_tads_lag = length(sig_tads_lag)
% 12 of 18 exps are significantly different. Exps 5, 7, 9, 13, 14, 15 are not. Exp
% 5 does not have stimtype = 1 (high multi). 13 has no respROIs. 

% Of the experiemnts with significant P val, how many tads are significant in each stim comparison?
for i = 1:size(which_mod_diff_lag,1)
    num_tads_diff_lag(i,1) = length( find((which_mod_diff_lag(i,:) ~= 0) & (which_mod_diff_lag(i,:) < 0.05)));
end

% bar plot the significant combos
labels = {'MS vs V', 'MS vs M', 'MS vs N', 'V vs M', 'V vs N', 'M vs N'}
bar(num_tads_diff_lag) %, labels)
set(gca, 'XTickLabel',labels, 'XTick',1:numel(labels))
ylabel('Num Sig Pairs')
ylim([0 12])
title('maxR lag time: Num of tads with P < 0.05 in mult comparison if P_{KW} < 0.05')

%% Combine across tadpoles to compare maxR vals across stim

%multi
tmp_data = [];
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxRlag_sq_MS')
        for i = 1:size(tadcluster{1,t}.respROIdff0_maxRlag_sq_MS,1)
            for j = 1:size(tadcluster{1,t}.respROIdff0_maxRlag_sq_MS,2)
                tmp_data = [tmp_data, tadcluster{1,t}.respROIdff0_maxRlag_sq_MS(i, (j+1):end)];
            end
        end                     
    end
end
maxRlag_allrespROI_MS = tmp_data;

% vis
tmp_data = [];
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxRlag_sq_V')
        for i = 1:size(tadcluster{1,t}.respROIdff0_maxRlag_sq_V,1)
            for j = 1:size(tadcluster{1,t}.respROIdff0_maxRlag_sq_V,2)
                tmp_data = [tmp_data, tadcluster{1,t}.respROIdff0_maxRlag_sq_V(i, (j+1):end)];
            end
        end                     
    end
end
maxRlag_allrespROI_V = tmp_data;

% mech
tmp_data = [];
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxRlag_sq_M')
        for i = 1:size(tadcluster{1,t}.respROIdff0_maxRlag_sq_M,1)
            for j = 1:size(tadcluster{1,t}.respROIdff0_maxRlag_sq_M,2)
                tmp_data = [tmp_data, tadcluster{1,t}.respROIdff0_maxRlag_sq_M(i, (j+1):end)];
            end
        end                     
    end
end
maxRlag_allrespROI_M = tmp_data;

%no stim
tmp_data = [];
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'respROIdff0_maxRlag_sq_N')
        for i = 1:size(tadcluster{1,t}.respROIdff0_maxRlag_sq_N,1)
            for j = 1:size(tadcluster{1,t}.respROIdff0_maxRlag_sq_N,2)
                tmp_data = [tmp_data, tadcluster{1,t}.respROIdff0_maxRlag_sq_N(i, (j+1):end)];
            end
        end                     
    end
end
maxRlag_allrespROI_N = tmp_data;
clear('tmp_data')

% make ECDF of all combined data
figure;
hold on
ecdf(maxRlag_allrespROI_MS)
ecdf(maxRlag_allrespROI_V)
ecdf(maxRlag_allrespROI_M)
ecdf(maxRlag_allrespROI_N)
title('ECDF of maxR lag time all tads all respROIs')
xlabel('maxR lag time')
ylabel('ROI proportion')
annotation('textbox', 'Position', [0.2 0.75 .1 .1], 'String', ['Multi'], 'Color', 'b', 'LineStyle', 'none' );
annotation('textbox', 'Position', [0.2 0.7 .1 .1], 'String', ['Vis'], 'Color', 'r', 'LineStyle', 'none' );
annotation('textbox', 'Position', [0.2 0.65 .1 .1], 'String', ['Mech'], 'Color', 'y', 'LineStyle', 'none' );
annotation('textbox', 'Position', [0.2 0.6 .1 .1], 'String', ['No stim'], 'Color', [0.5 0 0.5], 'LineStyle', 'none' );
annotation('textbox', 'Position', [0.2 0.55 .1 .1], 'String', 'N = 16, n = 474', 'Color', 'k', 'LineStyle', 'none');
saveas(gcf, 'ECDF of maxR lag time all tads all respROIs', 'png')

% stat test all combined data
% Statistical difference in any tad?
all_mods = [maxRlag_allrespROI_MS', maxRlag_allrespROI_V', maxRlag_allrespROI_M', maxRlag_allrespROI_N'];
size(all_mods)
[maxRlag_p, maxRlag_tbl, maxRlag_stats] = kruskalwallis(all_mods)
maxRlag_c = multcompare(maxRlag_stats)
% MS not diff from V, all other combos p < 0.05









%% Correlations by modality
% This code will take the basic correlation of each point in all trials
% with the same modality. 
% Importantly, this is "cleaned up" data. Previous iteration used all data
% from all trials. This version removes ROIs that don't respond to
% anything. 

% Start with tadcluster_analysis_20170809

%% Find and remove cells that don't respond at all
% response = peak is after stim onset (20 frames is a generous buffer), greater than 0.2 df/f0
%Criteria for elimination (e.g. these represent "bad" trials):
% -	Peak location is earlier than the stimulus onset
% -	Area is negative
% -	Peak is negative
% - peak is greater than 10

for t = 1:length(tadcluster)
[ tadcluster{1,t}.boolean_response, tadcluster{1,t}.sum_responses ] = get_respondingROIs3( tadcluster{1,t}.area_bytrial, tadcluster{1,t}.peak_bytrial, tadcluster{1,t}.peakloc_bytrial )
end

% if sum_responses = 0, then eliminate the ROI from further analysis. 
% This leaves some experiments with very few ROIs. Maybe eliminate them?

%% Reapply the same analysis as found in cluster_analysis_v1, lines 158-203 (over all trials)
% only use ROIs that have at least 1 response. 

% All ROIs have a vector of all trials in tadcluster{1,t}.alldff0

% Index the ROIs with responses
for t=1:length(tadcluster)
    tadcluster{1,t}.resp_ROIs = find(tadcluster{1,t}.sum_responses)
end

% calculate correlation coefficient for responding ROIs using all trials
for t = 1:length(tadcluster)
    [tadcluster{1,t}.respROIdff0_R, tadcluster{1,t}.respROIdff0_P] = corrcoef(tadcluster{1,t}.alldff0(tadcluster{1,t}.resp_ROIs,:)');
end

% any signficant?
% Is p < 0.05 for any ROI combinations?
for t = 1:length(tadcluster)
    tadcluster{1,t}.respROIsigP = find(tadcluster{1,t}.respROIdff0_P < 0.05);
    tadcluster{1,t}.respROIsigP_total = length(tadcluster{1,t}.respROIsigP)
end

for t = 1:length(tadcluster)
    respROIsigPcount(t) = tadcluster{1,t}.respROIsigP_total
end

% Plot correlation coefficients as a matrix
for t = 1:length(tadcluster)
    figure;
    colormap('hot')
    colormap(flipud(colormap))
    imagesc(tadcluster{1,t}.respROIdff0_R)
    colorbar
    title(sprintf('tad %d responding cells correlations (R) using all df/f0', t))
    fig_filename = sprintf('tad %d Responding ROIs correlation', t);
    saveas(gcf,fig_filename,'png');
    close;
end

% Plot P vals of correlation coefficients as a matrix
for t = 1:length(tadcluster)
    figure;
    colormap('hot')
    %colormap(flipud(colormap))
    imagesc(tadcluster{1,t}.respROIdff0_R)
    colorbar
    title(sprintf('tad %d responding cells P vals of correlations using all df/f0', t))
    fig_filename = sprintf('tad %d Responding ROIs correlation P val', t);
    saveas(gcf,fig_filename,'png');
    close;
end

%% Reapply the same analysis as found in cluster_analysis_v1, lines 235-377 (over each modality type separately)

% trial data is found in tadcluster{1,t}.dff0_multi ...vis, ...mech,
% ...none
% this only uses high/high trials
    
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
    [tadcluster{1,t}.respROIdff0_bystimtype_MS_R, tadcluster{1,t}.respROIdff0_bystimtype_MS_P] = corrcoef(tadcluster{1,t}.dff0_multi(tadcluster{1,t}.resp_ROIs,:)');
    [tadcluster{1,t}.respROIdff0_bystimtype_V_R, tadcluster{1,t}.respROIdff0_bystimtype_V_P] = corrcoef(tadcluster{1,t}.dff0_vis(tadcluster{1,t}.resp_ROIs,:)');
    [tadcluster{1,t}.respROIdff0_bystimtype_M_R, tadcluster{1,t}.respROIdff0_bystimtype_M_P] = corrcoef(tadcluster{1,t}.dff0_mech(tadcluster{1,t}.resp_ROIs,:)');
    [tadcluster{1,t}.respROIdff0_bystimtype_N_R, tadcluster{1,t}.respROIdff0_bystimtype_N_P] = corrcoef(tadcluster{1,t}.dff0_none(tadcluster{1,t}.resp_ROIs,:)');
    end
end

% Plot P vals of correlation coefficients as a matrix
% multisensory
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        figure;
        colormap('hot')
        %colormap(flipud(colormap))
        imagesc(tadcluster{1,t}.respROIdff0_bystimtype_MS_P)
        colorbar
        title(sprintf('tad %d responding cells P vals correlations using Multi df/f0', t))
        fig_filename = sprintf('tad %d responding cells Multi correlation P val', t)
        saveas(gcf,fig_filename,'png');
        close;
    end
end

%visual
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        figure;
        colormap('hot')
        %colormap(flipud(colormap))
        imagesc(tadcluster{1,t}.respROIdff0_bystimtype_V_P)
        colorbar
        title(sprintf('tad %d responding cells P vals correlations using Vis df/f0', t))
        fig_filename = sprintf('tad %d responding cells Vis correlation P val', t)
        saveas(gcf,fig_filename,'png');
        close;
    end
end

% mech
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        figure;
        colormap('hot')
        %colormap(flipud(colormap))
        imagesc(tadcluster{1,t}.respROIdff0_bystimtype_M_P)
        colorbar
        title(sprintf('tad %d responding cells P vals correlations using Mech df/f0', t))
        fig_filename = sprintf('tad %d responding cells Mech correlation P val', t)
        saveas(gcf,fig_filename,'png');
        close;
    end
end

% none
        colormap('hot')
        %colormap(flipud(colormap))
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        figure;
        colormap('hot')
        %colormap(flipud(colormap))
        imagesc(tadcluster{1,t}.respROIdff0_bystimtype_N_P)
        colorbar
        title(sprintf('tad %d responding cells P vals correlations using no stim df/f0', t))
        fig_filename = sprintf('tad %d responding cells no stim correlation P val', t)
        saveas(gcf,fig_filename,'png');
        close;
    end
end

% distribution of P values for each stim type (1 figure per tadpole)
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        figure;
        subplot(2,2,1)
        hist(tadcluster{1,t}.respROIdff0_bystimtype_MS_P,40)
        title('Multi')
        subplot(2,2,2)
        hist(tadcluster{1,t}.respROIdff0_bystimtype_V_P,40)
        title('Vis')
        subplot(2,2,3)
        hist(tadcluster{1,t}.respROIdff0_bystimtype_M_P,40)
        title('Mech')
        subplot(2,2,4)
        hist(tadcluster{1,t}.respROIdff0_bystimtype_N_P,40)
        title('None')
        suptitle(sprintf('tad %d responding cells hist of P vals by stimtype', t))
        fig_filename = sprintf('tad %d responding cells hist of Pvals of correlation by stimtype', t)
        saveas(gcf,fig_filename,'png');
        close;
    end
end


% Plot correlation coefficients (R) as a matrix
% multisensory
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        figure;
        colormap('hot')
        colormap(flipud(colormap))
        imagesc(tadcluster{1,t}.respROIdff0_bystimtype_MS_R)
        colorbar
        title(sprintf('tad %d responding cells correlations (R) using Multi df/f0', t))
        fig_filename = sprintf('tad %d responding cells Multi correlation', t)
        saveas(gcf,fig_filename,'png');
        close;
    end
end

%visual
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        figure;
        colormap('hot')
        colormap(flipud(colormap))
        imagesc(tadcluster{1,t}.respROIdff0_bystimtype_V_R)
        colorbar
        title(sprintf('tad %d responding cells correlations (R) using Vis df/f0', t))
        fig_filename = sprintf('tad %d responding cells Vis correlation', t)
        saveas(gcf,fig_filename,'png');
        close;
    end
end

% mech
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        figure;
        colormap('hot')
        colormap(flipud(colormap))
        imagesc(tadcluster{1,t}.respROIdff0_bystimtype_M_R)
        colorbar
        title(sprintf('tad %d responding cells correlations (R) using Mech df/f0', t))
        fig_filename = sprintf('tad %d responding cells Mech correlation', t)
        saveas(gcf,fig_filename,'png');
        close;
    end
end

% none
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        figure;
        colormap('hot')
        colormap(flipud(colormap))
        imagesc(tadcluster{1,t}.respROIdff0_bystimtype_N_R)
        colorbar
        title(sprintf('tad %d responding cells correlations (R) using no stim df/f0', t))
        fig_filename = sprintf('tad %d responding cells no stim correlation', t)
        saveas(gcf,fig_filename,'png');
        close;
    end
end

% distribution of R values for each stim type (1 figure per tadpole)
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        figure;
        subplot(2,2,1)
        hist(tadcluster{1,t}.respROIdff0_bystimtype_MS_R,40)
        title('Multi')
        subplot(2,2,2)
        hist(tadcluster{1,t}.respROIdff0_bystimtype_V_R,40)
        title('Vis')
        subplot(2,2,3)
        hist(tadcluster{1,t}.respROIdff0_bystimtype_M_R,40)
        title('Mech')
        subplot(2,2,4)
        hist(tadcluster{1,t}.respROIdff0_bystimtype_N_R,40)
        title('None')
        suptitle(sprintf('tad %d responding cells hist of R vals by stimtype', t))
        fig_filename = sprintf('tad %d responding cells hist of R vals of correlation by stimtype', t)
        saveas(gcf,fig_filename,'png');
        close;
    end
end

%% Plot correlation (R) vs. Euclidean distance for responding ROIs

% distance is in tadcluster{1,t}.ROIdist

%%%%%%This doesn't work. The fit is in the wrong place and all the P vals
%%%%%%seem to be very close together. Also dimension issues with t > 1

% plot distance vs P val for each ROI pair with linear trendline
for t = 1:length(tadcluster)
    figure;
for i = 1:(size(tadcluster{1,t}.respROIdff0_P,1)-1)
    hold on
    plot(tadcluster{1,t}.ROIdist(tadcluster{1,t}.resp_ROIs,:), tadcluster{1,t}.respROIdff0_P(i+1,:), 'ko')
    %all_vals(i,:) = tadcluster{1,t}.ROIdist(tadcluster{1,t}.resp_ROIs(i),:), tadcluster{1,t}.respROIdff0_P(i+1,:)
end
    % generate trendline
    %myfit = polyfit(all_vals(:,1), all_vals(:,2), 1)
    %plot(myfit, 'r', 'LineWidth', 5)
    hold off
    title(sprintf('tad %d responding cells ROI distance vs P value df/f0', t));
    fig_filename = sprintf('tad %d responding cells ROI distance vs P value', t);
    xlabel('ROI distance (pixels)')
    ylabel('P value')
    saveas(gcf,fig_filename,'png');
    close;
end

%% Divide correlations into groups (low/mid/high) and plot distance

for t = 1:length(tadcluster)
    % get Pvals for all modalities combined
    tadcluster{1,t}.respROI_Pvals_low = get_subset_pvals(0, 0.01, tadcluster{1,t}.respROIdff0_P);
    tadcluster{1,t}.respROI_Pvals_med = get_subset_pvals(0.4, 0.6, tadcluster{1,t}.respROIdff0_P);
    tadcluster{1,t}.respROI_Pvals_high = get_subset_pvals(0.8, 1, tadcluster{1,t}.respROIdff0_P);
    
    % get pairwise distances for each Pval subset
    tadcluster{1,t}.respROI_Pvals_low_pairdist = get_pairwise_distances(tadcluster{1,t}.ROIdist, tadcluster{1,t}.respROI_Pvals_low, 1, tadcluster{1,t}.resp_ROIs);
    tadcluster{1,t}.respROI_Pvals_med_pairdist = get_pairwise_distances(tadcluster{1,t}.ROIdist, tadcluster{1,t}.respROI_Pvals_med, 1, tadcluster{1,t}.resp_ROIs);
    tadcluster{1,t}.respROI_Pvals_high_pairdist = get_pairwise_distances(tadcluster{1,t}.ROIdist, tadcluster{1,t}.respROI_Pvals_high, 1, tadcluster{1,t}.resp_ROIs);
end 

% Plot the distance between the cells in each category
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.respROI_Pvals_low_pairdist)) > 1 && sum(size(tadcluster{1,t}.respROI_Pvals_med_pairdist)) > 1 && sum(size(tadcluster{1,t}.respROI_Pvals_high_pairdist)) > 1
    % create X vals
    xval_low = ones(size(tadcluster{1,t}.respROI_Pvals_low_pairdist,1), 1)
    xval_med = 2* ones(size(tadcluster{1,t}.respROI_Pvals_med_pairdist,1), 1)
    xval_high = 3* ones(size(tadcluster{1,t}.respROI_Pvals_high_pairdist,1), 1)
    figure;
    hold on
    scatter(xval_low, tadcluster{1,t}.respROI_Pvals_low_pairdist(:,3))
    scatter(xval_med, tadcluster{1,t}.respROI_Pvals_med_pairdist(:,3))
    scatter(xval_high, tadcluster{1,t}.respROI_Pvals_high_pairdist(:,3))
    %errorbar(avgof_avgpeakloc(1:4), stdof_avgpeakloc(1:4), 'ob', 'linewidth', 2, 'Markersize', 10, 'MarkerFaceColor', 'b') % mean with stdev

    hold off
    title(sprintf('ROI distances by P value exp %d', t))
    ax=gca;
    ax.XTick = [1, 2, 3]
    ax.XTickLabel = {'low', 'med', 'high'}
    xlabel('P val category')
    ylabel('distance (pixels)')
    xlim([0.5 3.5])
    fig_filename = sprintf('ROI distances by P value exp %d', t)
    saveas(gcf,fig_filename,'png');
    close;
    end
end

%% Repeat dividing correlations into groups by p values but using by-modality stimulus presented data

for t = 1:length(tadcluster)
    % get Pvals for all modalities combined
    if length(find(tadcluster{1,t}.stimorder == 1)) > 1 %this eliminates experiments without high/high MS trials
        %Multi
        tadcluster{1,t}.respROI_bystimtype_MS_Pvals_low = get_subset_pvals(0, 0.01, tadcluster{1,t}.respROIdff0_bystimtype_MS_P);
        tadcluster{1,t}.respROI_bystimtype_MS_Pvals_med = get_subset_pvals(0.4, 0.6, tadcluster{1,t}.respROIdff0_bystimtype_MS_P);
        tadcluster{1,t}.respROI_bystimtype_MS_Pvals_high = get_subset_pvals(0.8, 1, tadcluster{1,t}.respROIdff0_bystimtype_MS_P);
        %Vis
        tadcluster{1,t}.respROI_bystimtype_V_Pvals_low = get_subset_pvals(0, 0.01, tadcluster{1,t}.respROIdff0_bystimtype_V_P);
        tadcluster{1,t}.respROI_bystimtype_V_Pvals_med = get_subset_pvals(0.4, 0.6, tadcluster{1,t}.respROIdff0_bystimtype_V_P);
        tadcluster{1,t}.respROI_bystimtype_V_Pvals_high = get_subset_pvals(0.8, 1, tadcluster{1,t}.respROIdff0_bystimtype_V_P);
        %Mech
        tadcluster{1,t}.respROI_bystimtype_M_Pvals_low = get_subset_pvals(0, 0.01, tadcluster{1,t}.respROIdff0_bystimtype_M_P);
        tadcluster{1,t}.respROI_bystimtype_M_Pvals_med = get_subset_pvals(0.4, 0.6, tadcluster{1,t}.respROIdff0_bystimtype_M_P);
        tadcluster{1,t}.respROI_bystimtype_M_Pvals_high = get_subset_pvals(0.8, 1, tadcluster{1,t}.respROIdff0_bystimtype_M_P);
        %none
        tadcluster{1,t}.respROI_bystimtype_N_Pvals_low = get_subset_pvals(0, 0.01, tadcluster{1,t}.respROIdff0_bystimtype_N_P);
        tadcluster{1,t}.respROI_bystimtype_N_Pvals_med = get_subset_pvals(0.4, 0.6, tadcluster{1,t}.respROIdff0_bystimtype_N_P);
        tadcluster{1,t}.respROI_bystimtype_N_Pvals_high = get_subset_pvals(0.8, 1, tadcluster{1,t}.respROIdff0_bystimtype_N_P);

        % get pairwise distances for each Pval subset
        %multi
        tadcluster{1,t}.respROI_bystimtype_MS_Pvals_low_pairdist = get_pairwise_distances(tadcluster{1,t}.ROIdist, tadcluster{1,t}.respROI_bystimtype_MS_Pvals_low, 1, tadcluster{1,t}.resp_ROIs);
        tadcluster{1,t}.respROI_bystimtype_MS_Pvals_med_pairdist = get_pairwise_distances(tadcluster{1,t}.ROIdist, tadcluster{1,t}.respROI_bystimtype_MS_Pvals_med, 1, tadcluster{1,t}.resp_ROIs);
        tadcluster{1,t}.respROI_bystimtype_MS_Pvals_high_pairdist = get_pairwise_distances(tadcluster{1,t}.ROIdist, tadcluster{1,t}.respROI_bystimtype_MS_Pvals_high, 1, tadcluster{1,t}.resp_ROIs);
        %vis
        tadcluster{1,t}.respROI_bystimtype_V_Pvals_low_pairdist = get_pairwise_distances(tadcluster{1,t}.ROIdist, tadcluster{1,t}.respROI_bystimtype_V_Pvals_low, 1, tadcluster{1,t}.resp_ROIs);
        tadcluster{1,t}.respROI_bystimtype_V_Pvals_med_pairdist = get_pairwise_distances(tadcluster{1,t}.ROIdist, tadcluster{1,t}.respROI_bystimtype_V_Pvals_med, 1, tadcluster{1,t}.resp_ROIs);
        tadcluster{1,t}.respROI_bystimtype_V_Pvals_high_pairdist = get_pairwise_distances(tadcluster{1,t}.ROIdist, tadcluster{1,t}.respROI_bystimtype_V_Pvals_high, 1, tadcluster{1,t}.resp_ROIs);
        %mech
        tadcluster{1,t}.respROI_bystimtype_M_Pvals_low_pairdist = get_pairwise_distances(tadcluster{1,t}.ROIdist, tadcluster{1,t}.respROI_bystimtype_M_Pvals_low, 1, tadcluster{1,t}.resp_ROIs);
        tadcluster{1,t}.respROI_bystimtype_M_Pvals_med_pairdist = get_pairwise_distances(tadcluster{1,t}.ROIdist, tadcluster{1,t}.respROI_bystimtype_M_Pvals_med, 1, tadcluster{1,t}.resp_ROIs);
        tadcluster{1,t}.respROI_bystimtype_M_Pvals_high_pairdist = get_pairwise_distances(tadcluster{1,t}.ROIdist, tadcluster{1,t}.respROI_bystimtype_M_Pvals_high, 1, tadcluster{1,t}.resp_ROIs);
        %none
        tadcluster{1,t}.respROI_bystimtype_N_Pvals_low_pairdist = get_pairwise_distances(tadcluster{1,t}.ROIdist, tadcluster{1,t}.respROI_bystimtype_N_Pvals_low, 1, tadcluster{1,t}.resp_ROIs);
        tadcluster{1,t}.respROI_bystimtype_N_Pvals_med_pairdist = get_pairwise_distances(tadcluster{1,t}.ROIdist, tadcluster{1,t}.respROI_bystimtype_N_Pvals_med, 1, tadcluster{1,t}.resp_ROIs);
        tadcluster{1,t}.respROI_bystimtype_N_Pvals_high_pairdist = get_pairwise_distances(tadcluster{1,t}.ROIdist, tadcluster{1,t}.respROI_bystimtype_N_Pvals_high, 1, tadcluster{1,t}.resp_ROIs);
    end
end

% Plot the distance between the cells in each category
for t = 1:length(tadcluster)
    if length(find(tadcluster{1,t}.stimorder == 1)) > 1 %this eliminates experiments without high/high MS trials
        if ~isempty(tadcluster{1,t}.respROI_bystimtype_N_Pvals_low_pairdist) && ~isempty(tadcluster{1,t}.respROI_bystimtype_N_Pvals_med_pairdist) && ~isempty(tadcluster{1,t}.respROI_bystimtype_N_Pvals_high_pairdist)...
                && ~isempty(tadcluster{1,t}.respROI_bystimtype_MS_Pvals_low_pairdist) && ~isempty(tadcluster{1,t}.respROI_bystimtype_MS_Pvals_med_pairdist) && ~isempty(tadcluster{1,t}.respROI_bystimtype_MS_Pvals_high_pairdist) % eliminate exps with no P vals in a category.
            figure;
            %multi
            subplot(2,2,1)
            hold on
            scatter(ones(size(tadcluster{1,t}.respROI_bystimtype_MS_Pvals_low_pairdist,1), 1), tadcluster{1,t}.respROI_bystimtype_MS_Pvals_low_pairdist(:,3))
            scatter((2* ones(size(tadcluster{1,t}.respROI_bystimtype_MS_Pvals_med_pairdist,1), 1)), tadcluster{1,t}.respROI_bystimtype_MS_Pvals_med_pairdist(:,3))
            scatter((3* ones(size(tadcluster{1,t}.respROI_bystimtype_MS_Pvals_high_pairdist,1), 1)), tadcluster{1,t}.respROI_bystimtype_MS_Pvals_high_pairdist(:,3))
            hold off
            title('multi')
            ax=gca;
            ax.XTick = [1, 2, 3]
            ax.XTickLabel = {'low', 'med', 'high'}
            xlabel('P val category')
            ylabel('distance (pixels)')
            xlim([0.5 3.5])

            %vis
            subplot(2,2,2)
            hold on
            scatter(ones(size(tadcluster{1,t}.respROI_bystimtype_V_Pvals_low_pairdist,1), 1), tadcluster{1,t}.respROI_bystimtype_V_Pvals_low_pairdist(:,3))
            scatter((2* ones(size(tadcluster{1,t}.respROI_bystimtype_V_Pvals_med_pairdist,1), 1)), tadcluster{1,t}.respROI_bystimtype_V_Pvals_med_pairdist(:,3))
            scatter((3* ones(size(tadcluster{1,t}.respROI_bystimtype_V_Pvals_high_pairdist,1), 1)), tadcluster{1,t}.respROI_bystimtype_V_Pvals_high_pairdist(:,3))
            hold off
            title('visual')
            ax=gca;
            ax.XTick = [1, 2, 3]
            ax.XTickLabel = {'low', 'med', 'high'}
            xlabel('P val category')
            ylabel('distance (pixels)')
            xlim([0.5 3.5])

            %mech
            subplot(2,2,3)
            hold on
            scatter(ones(size(tadcluster{1,t}.respROI_bystimtype_M_Pvals_low_pairdist,1), 1), tadcluster{1,t}.respROI_bystimtype_M_Pvals_low_pairdist(:,3))
            scatter((2* ones(size(tadcluster{1,t}.respROI_bystimtype_M_Pvals_med_pairdist,1), 1)), tadcluster{1,t}.respROI_bystimtype_M_Pvals_med_pairdist(:,3))
            scatter((3* ones(size(tadcluster{1,t}.respROI_bystimtype_M_Pvals_high_pairdist,1), 1)), tadcluster{1,t}.respROI_bystimtype_M_Pvals_high_pairdist(:,3))
            hold off
            title('mech')
            ax=gca;
            ax.XTick = [1, 2, 3]
            ax.XTickLabel = {'low', 'med', 'high'}
            xlabel('P val category')
            ylabel('distance (pixels)')
            xlim([0.5 3.5])

            %none
            subplot(2,2,4)
            hold on
            scatter(ones(size(tadcluster{1,t}.respROI_bystimtype_N_Pvals_low_pairdist,1), 1), tadcluster{1,t}.respROI_bystimtype_N_Pvals_low_pairdist(:,3))
            scatter((2* ones(size(tadcluster{1,t}.respROI_bystimtype_N_Pvals_med_pairdist,1), 1)), tadcluster{1,t}.respROI_bystimtype_N_Pvals_med_pairdist(:,3))
            scatter((3* ones(size(tadcluster{1,t}.respROI_bystimtype_N_Pvals_high_pairdist,1), 1)), tadcluster{1,t}.respROI_bystimtype_N_Pvals_high_pairdist(:,3))
            hold off
            title('none')
            ax=gca;
            ax.XTick = [1, 2, 3]
            ax.XTickLabel = {'low', 'med', 'high'}
            xlabel('P val category')
            ylabel('distance (pixels)')
            xlim([0.5 3.5])

            suptitle(sprintf('ROI distances by P value split by modality presented exp %d', t))
            fig_filename = sprintf('ROI distances by P value split by modality presented exp %d', t)
            saveas(gcf,fig_filename,'png');
            close;
        end
    end
    
end








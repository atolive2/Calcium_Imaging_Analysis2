%% SFN 2017 Poster Figures
% begin with load('allData_workspace_20171107.mat')

%% Tectum-shaped scatter plot of highcorr ROIs
% Make a movie of the 3 modalities

% multi tectum shaped scatterplot (show resp_ROIs and highcorr ROIs only)
for t = 1:length(allData)
    if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
            figure;
            hold on
            scatter(allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs, 1), allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs, 2), 65, 'o', 'MarkerEdgeColor', [0.5 0.5 0.5])
            scatter(allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI, 1), allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI, 2), 100, 'o', 'filled', 'MarkerEdgeColor', [0.5 0 0.5], 'MarkerFaceColor', [0.5 0 0.5])
            hold off
            xlim([0 500])
            ylim([0 500])
            fig_filename = sprintf('tectum shaped plot of highcorr ROIs tad %d stage %d MS', allData{1,t}.expnum, allData{1,t}.stage)
            saveas(gcf, fig_filename, 'tiff')
            close;
        end
    end
end

%% Visual tectum shaped scatterplot (show resp_ROIs and highcorr ROIs only)

tmp =  [];
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_V_common_AROI')
    for r = 1:length(allData{1,t}.correlated_ROIs_dff0_V_common_AROI)
        tmp = [tmp; cell2mat(allData{1,t}.correlated_ROIs_dff0_V_common_AROI(r))];
    end
    allData{1,t}.uniqueHighCorrROI_V =  unique(tmp)
    tmp = [];
    end
end

for t = 1:length(allData)
    if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_V_common_AROI')
            figure;
            hold on
            scatter(allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs, 1), allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs, 2), 65, 'o', 'MarkerEdgeColor', [0.5 0.5 0.5])
            scatter(allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_V, 1), allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_V, 2), 100, 'o', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
            hold off
            xlim([0 500])
            ylim([0 500])
            fig_filename = sprintf('tectum shaped plot of highcorr ROIs tad %d stage %d V', allData{1,t}.expnum, allData{1,t}.stage)
            saveas(gcf, fig_filename, 'tiff')
            close;
        end
    end
end

%% Mechanosensory tectum shaped scatterplot (show resp_ROIs and highcorr ROIs only)

tmp =  [];
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_M_common_AROI')
    for r = 1:length(allData{1,t}.correlated_ROIs_dff0_M_common_AROI)
        tmp = [tmp; cell2mat(allData{1,t}.correlated_ROIs_dff0_M_common_AROI(r))];
    end
    allData{1,t}.uniqueHighCorrROI_M =  unique(tmp)
    tmp = [];
    end
end

for t = 1:length(allData)
    if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_M_common_AROI')
            figure;
            hold on
            scatter(allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs, 1), allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs, 2), 65, 'o', 'MarkerEdgeColor', [0.5 0.5 0.5])
            scatter(allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_M, 1), allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_M, 2), 100, 'o', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b')
            hold off
            xlim([0 500])
            ylim([0 500])
            fig_filename = sprintf('tectum shaped plot of highcorr ROIs tad %d stage %d M', allData{1,t}.expnum, allData{1,t}.stage)
            saveas(gcf, fig_filename, 'tiff')
            close;
        end
    end
end

%% Overlpping cells plotted together
% Plot all 3 modalities, then note ROIs that are in multiple groups

% Find cells that are high corr in multiple stimulus conditions
for t = 1:length(allData)
    if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
            if isfield(allData{1,t}, 'correlated_ROIs_dff0_M_common_AROI')
                if isfield(allData{1,t}, 'correlated_ROIs_dff0_V_common_AROI')
                    V_rois = [];
                    MS_rois = [];
                    M_rois = [];
                    for i = 1:length(allData{1,t}.correlated_ROIs_dff0_V_common_AROI)
                        V_rois = [V_rois; cell2mat(allData{1,t}.correlated_ROIs_dff0_V_common_AROI(i))];
                    end
                    for i = 1:length(allData{1,t}.correlated_ROIs_dff0_MS_common_AROI)
                        MS_rois = [MS_rois; cell2mat(allData{1,t}.correlated_ROIs_dff0_MS_common_AROI(i))];
                    end
                    for i = 1:length(allData{1,t}.correlated_ROIs_dff0_M_common_AROI)
                        M_rois = [M_rois; cell2mat(allData{1,t}.correlated_ROIs_dff0_M_common_AROI(i))];
                    end
                    V_roisU = unique(V_rois)
                    MS_roisU = unique(MS_rois)
                    M_roisU = unique(M_rois)
                    allData{1,t}.highcorr_M_V = intersect(V_roisU, M_roisU);
                    allData{1,t}.highcorr_M_MS = intersect(MS_roisU, M_roisU);
                    allData{1,t}.highcorr_MS_V = intersect(V_roisU, MS_roisU);
                    allData{1,t}.highcorr_M_V_MS = intersect(allData{1,t}.highcorr_MS_V, allData{1,t}.highcorr_M_MS);
                end
            end
        end
    end
end



% Scatter plot data for tads with highcorr in all 3 modalities
for t = 1:length(allData)
    if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
            if isfield(allData{1,t}, 'correlated_ROIs_dff0_M_common_AROI')
                if isfield(allData{1,t}, 'correlated_ROIs_dff0_V_common_AROI')
                    figure;
                    hold on
                    %all resp_ROIs as gray empty circles
                    scatter(allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs, 1), allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs, 2), 65, 'o', 'MarkerEdgeColor', [0.5 0.5 0.5])
                    %all highcorr_M ROIs as blue filled circles
                    scatter(allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_M, 1), allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_M, 2), 100, 'o', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b')
                    %all highcorr_V ROIs as red filled circles
                    scatter(allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_V, 1), allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_V, 2), 100, 'o', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
                    %all highcorr_MS ROIs as yellow filled circles
                    scatter(allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI, 1), allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI, 2), 100, 'o', 'filled', 'MarkerEdgeColor', 'y', 'MarkerFaceColor', 'y')
                    % all M and V highcorr ROIs as purple filled circles
                    scatter(allData{1,t}.ROIcenters(allData{1,t}.highcorr_M_V, 1), allData{1,t}.ROIcenters(allData{1,t}.highcorr_M_V, 2), 100, 'o', 'filled', 'MarkerEdgeColor', [0.5 0 0.5], 'MarkerFaceColor', [0.5 0 0.5])
                    % all M and MS highcorr ROIs as green filled circles
                    scatter(allData{1,t}.ROIcenters(allData{1,t}.highcorr_M_MS, 1), allData{1,t}.ROIcenters(allData{1,t}.highcorr_M_MS, 2), 100, 'o', 'filled', 'MarkerEdgeColor', [0.4 0.5 0.2], 'MarkerFaceColor', [0.4 0.5 0.2])
                    % all V and MS highcorr ROIs as orange filled circles
                    scatter(allData{1,t}.ROIcenters(allData{1,t}.highcorr_M_V, 1), allData{1,t}.ROIcenters(allData{1,t}.highcorr_M_V, 2), 100, 'o', 'filled', 'MarkerEdgeColor', [1 0.6 0.2], 'MarkerFaceColor', [1 0.6 0.2])
                    % all M, V and MS highcorr ROIs as black filled circles
                    scatter(allData{1,t}.ROIcenters(allData{1,t}.highcorr_M_V_MS, 1), allData{1,t}.ROIcenters(allData{1,t}.highcorr_M_V_MS, 2), 100, 'o', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
                    hold off
                    xlim([0 500])
                    ylim([0 500])
                    fig_filename = sprintf('tectum shaped plot of highcorr ROIs tad %d (t%d) stage %d overlap', allData{1,t}.expnum, t, allData{1,t}.stage)
                    saveas(gcf, fig_filename, 'tiff')
                    close;
                end
            end
        end
    end
end

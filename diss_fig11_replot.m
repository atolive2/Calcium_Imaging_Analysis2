%% Diss figure 11 replot

%MS
for t = 1:length(allData)
    if isfield(allData{1,t}, 'uniqueHighCorrROI_MS')
        figure;
        hold on
        scatter(allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs,1), allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs,2), 150, [0.7 0.7 0.7], 'filled') %all respROIs
        scatter(allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_MS, 1), allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_MS, 2), 150, [226/256 195/256 15/256], 'filled', 'MarkerEdgeColor', 'k')
        hold off
        set(gca, 'Visible', 'off')
        fig_filename = sprintf('Tad %d (t=%d) tectum shaped plot HC MS', allData{1,t}.expnum, t)
        saveas(gcf, fig_filename, 'png')
        saveas(gcf, fig_filename, 'epsc2')
        close;
    end
end

%V
for t = 1:length(allData)
    if isfield(allData{1,t}, 'uniqueHighCorrROI_V')
        figure;
        hold on
        scatter(allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs,1), allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs,2), 150, [0.7 0.7 0.7], 'filled') %all respROIs
        scatter(allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_V, 1), allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_V, 2), 150, 'r', 'filled', 'MarkerEdgeColor', 'k')
        hold off
        set(gca, 'Visible', 'off')
        fig_filename = sprintf('Tad %d (t=%d) tectum shaped plot HC V', allData{1,t}.expnum, t)
        saveas(gcf, fig_filename, 'png')
        saveas(gcf, fig_filename, 'epsc2')
        close;
    end
end

%M
for t = 1:length(allData)
    if isfield(allData{1,t}, 'uniqueHighCorrROI_M')
        figure;
        hold on
        scatter(allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs,1), allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs,2), 150, [0.7 0.7 0.7], 'filled') %all respROIs
        scatter(allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_M, 1), allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_M, 2), 150, [0.3010, 0.7450, 0.9330], 'filled', 'MarkerEdgeColor', 'k')
        hold off
        set(gca, 'Visible', 'off')
        fig_filename = sprintf('Tad %d (t=%d) tectum shaped plot HC M', allData{1,t}.expnum, t)
        saveas(gcf, fig_filename, 'png')
        saveas(gcf, fig_filename, 'epsc2')
        close;
    end
end

% Overlap
for t = 1:length(allData)
    if length(allData{1,t}.resp_ROIs) > 3
        % all resp ROIs in gray
        figure;
        hold on
        scatter(allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs,1), allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs,2), 150, [0.7 0.7 0.7], 'filled') %all respROIs
    %plot all HC - each mod
    if isfield(allData{1,t}, 'uniqueHighCorrROI_M')
        scatter(allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_M, 1), allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_M, 2), 150, [0.3010, 0.7450, 0.9330], 'filled', 'MarkerEdgeColor', 'k')
    end
    if isfield(allData{1,t}, 'uniqueHighCorrROI_V')
        scatter(allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_V, 1), allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_V, 2), 150, 'r', 'filled', 'MarkerEdgeColor', 'k')
    end
    if isfield(allData{1,t}, 'uniqueHighCorrROI_MS')
        scatter(allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_MS, 1), allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_MS, 2), 150, [226/256 195/256 15/256], 'filled', 'MarkerEdgeColor', 'k')
    end
    %plot overlaps
    if isfield(allData{1,t}, 'highcorr_M_V')
        scatter(allData{1,t}.ROIcenters(allData{1,t}.highcorr_M_V, 1), allData{1,t}.ROIcenters(allData{1,t}.highcorr_M_V, 2), 150, [0.5 0 0.5], 'filled', 'MarkerEdgeColor', 'k')
    end
    if isfield(allData{1,t}, 'highcorr_M_MS')
        scatter(allData{1,t}.ROIcenters(allData{1,t}.highcorr_M_V, 1), allData{1,t}.ROIcenters(allData{1,t}.highcorr_M_V, 2), 150, 'g', 'filled', 'MarkerEdgeColor', 'k')
    end
    if isfield(allData{1,t}, 'highcorr_MS_V')
        scatter(allData{1,t}.ROIcenters(allData{1,t}.highcorr_MS_V, 1), allData{1,t}.ROIcenters(allData{1,t}.highcorr_MS_V, 2), 150, [0.8500, 0.3250, 0.0980], 'filled', 'MarkerEdgeColor', 'k')
    end    
    if isfield(allData{1,t}, 'highcorr_M_V_MS')
        scatter(allData{1,t}.ROIcenters(allData{1,t}.highcorr_M_V_MS, 1), allData{1,t}.ROIcenters(allData{1,t}.highcorr_M_V_MS, 2), 150, 'k', 'filled', 'MarkerEdgeColor', 'k')
    end
    
    % finish the plot
    hold off
    set(gca, 'Visible', 'off')
    fig_filename = sprintf('Tad %d (t=%d) tectum shaped plot HC overlap', allData{1,t}.expnum, t)
    saveas(gcf, fig_filename, 'png')
    saveas(gcf, fig_filename, 'epsc2')
    close;
    end
end



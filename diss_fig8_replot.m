% replot diss figure 9 (HC topography)

%% Topographical distance ECDF

% get topographical distance between each pair of respROIs
for t = 1:length(allData)
    for r = 1:length(allData{1,t}.resp_ROIs)
        for s = 1:length(allData{1,t}.resp_ROIs)
            rr = allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs(r), :);
            ss = allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs(s), :);
            allData{1,t}.topo_dist(r,s) = sqrt( (rr(1) - ss(1)) ^2 + (rr(2) - ss(2)) ^2);
        end
    end
end

% sort by high corr status into 2 vectors
HC_dist = [];
SC_dist = [];
for t = 1:length(allData)
    if length(allData{1,t}.resp_ROIs) > 2
        for r = 1:length(allData{1,t}.resp_ROIs)
            for s = 1:length(allData{1,t}.resp_ROIs)
                if ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.allhighcorrROI)
                    if ismember(allData{1,t}.resp_ROIs(s), allData{1,t}.allhighcorrROI)
                        HC_dist = [HC_dist, allData{1,t}.topo_dist(r,s)];
                    else
                        SC_dist = [SC_dist, allData{1,t}.topo_dist(r,s)];
                    end
                else
                    SC_dist = [SC_dist, allData{1,t}.topo_dist(r,s)];
                end
            end
        end
    end
end

% Plot ECDF of distances
[f1, x1] = ecdf(HC_dist);
[f2, x2] = ecdf(SC_dist);
figure;
hold on
plot(x1, f1, 'Color', [0/256 73/256 73/256], 'LineWidth', 3)
plot(x2, f2, 'Color', [146/256 73/256 0/256], 'LineWidth', 3)
hold off
ylabel('ROI proportion')
xlabel('Distance (pixels)')
xlim([-inf inf])
set(gca,'FontSize',30)
fig_filename = 'ecdf hc_nhc topo distance'
saveas(gcf, fig_filename, 'png')
saveas(gcf, fig_filename, 'epsc2')

% histogram
figure;
histogram(HC_dist)
figure;
histogram(SC_dist)

% stats
[h, p] = kstest2(HC_dist, SC_dist)
% significant, p = 0

%sooo...this is backwards from the original. WTF???

%% plot a topographical map by HC status 

% run them all and pick a pretty one
% size of dot is proportional to the number of ROIs it is HC with
% gray = resp_ROI that is not HC. don't plot nonresponders

% get a value for how many HC with (use multisensory only)
for t = 1:length(allData)
    %if isfield(allData{1,t}, '
    for r = 1:length(allData{1,t}.resp_ROIs)
        allData{1,t}.HC_count(r) = length(find(allData{1,t}.respROIdff0_maxR_sq_MS(r,:) > 0.5));
    end
end

for t = 1:length(allData)
    if isfield(allData{1,t}, 'HC_count')
        if isfield(allData{1,t}, 'uniqueHighCorrROI_MS')
        lg_ct = max(allData{1,t}.HC_count)
        dot_size = (allData{1,t}.HC_count/lg_ct)+1;
        for r = 1:length(allData{1,t}.resp_ROIs)
            if ~ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.uniqueHighCorrROI_MS)
                dot_size(r) = 0;
            end
        end
        [row col dot_sizeHC] = find(dot_size)
        figure;
        hold on
        scatter(allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs,1), allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs,2), 50, [0.7 0.7 0.7], 'filled') %all respROIs
        scatter(allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_MS, 1), allData{1,t}.ROIcenters(allData{1,t}.uniqueHighCorrROI_MS, 2), dot_sizeHC * 100, [0/256 73/256 73/256], 'filled')
        hold off
        %annotation('textbox', 'Position', [0.14 0.16 .1 .1], 'String', ['All responding ROIs'], 'Color', [0.7 0.7 0.7], 'LineStyle', 'none', 'FontSize', 30 );
        %annotation('textbox', 'Position', [0.14 0.11 .1 .1], 'String', ['Highly correlated ROIs'], 'Color', 'b', 'LineStyle', 'none', 'FontSize', 30 );
        title(sprintf('tad %d (t=%d) high corr cells by size', allData{1,t}.expnum, t))
        %set(gca,'FontSize',30)
        set(gca, 'Visible', 'off')
        fig_filename = sprintf('Tad %d (t=%d) tectum shaped plot HC size by num_HC', allData{1,t}.expnum, t)
        saveas(gcf, fig_filename, 'png')
        saveas(gcf, fig_filename, 'epsc2')
        clear('dot_size', 'dot_sizeHC', 'lg_ct')
        close;
        end
    end
end



t = 12
for t = 1:length(allData)
    if isfield(allData{1,t}, 'highcorr_numROIs_MS')
        figure;
        hold on
        scatter(allData{1,t}.ROIcenters(:,1), allData{1,t}.ROIcenters(:,2), '+', 'MarkerEdgeColor', [0.5 0.5 0.5])
        scatter(allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs,1), allData{1,t}.ROIcenters(allData{1,t}.resp_ROIs,2), 50, [0.7 0.7 0.7], 'filled')
        scatter(allData{1,t}.ROIcenters(allData{1,t}.high_corr_ROIs,1), allData{1,t}.ROIcenters(allData{1,t}.high_corr_ROIs, 2), 50, 'm', 'filled')
        scatter(allData{1,t}.ROIcenters(allData{1,t}.low_corr_ROIs,1), allData{1,t}.ROIcenters(allData{1,t}.low_corr_ROIs, 2), 50, 'g', 'filled')
        hold off
        annotation('textbox', 'Position', [0.14 0.21 .1 .1], 'String', ['All recorded ROIs (+)'], 'Color', [0.5 0.5 0.5], 'LineStyle', 'none' );
        annotation('textbox', 'Position', [0.14 0.16 .1 .1], 'String', ['All responding ROIs (o)'], 'Color', [0.7 0.7 0.7], 'LineStyle', 'none' );
        annotation('textbox', 'Position', [0.14 0.11 .1 .1], 'String', ['High correlation ROIs'], 'Color', 'm', 'LineStyle', 'none' );
        annotation('textbox', 'Position', [0.14 0.07 .1 .1], 'String', ['Low correlation ROIs'], 'Color', 'g', 'LineStyle', 'none' );
        title(sprintf('tad %d high and low correlation ROIs', t))
        fig_filename = sprintf('tad %d high and low correlation ROIs tectum shaped', t)
        saveas(gcf, fig_filename, 'png')
        close;
    end
end

























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


%% Combining all cells, are there differences by stage in basic params?
st46 = find(allRespROIs(:,28) == 46);
st49 = find(allRespROIs(:,28) == 49);
labels = {'area MS'; 'area V'; 'area M'; 'area NS'; 'peak MS'; 'peak V'; 'peak M'; 'peak NS'; ... 
    'peak loc MS'; 'peak loc V'; 'peak loc M'; 'peak loc NS'; 'MSEnh peak'; 'unimax peak'; 'unimax stimtype'; ... 
    'onset time MS'; 'onset time V'; 'onset time M'; 'onset time NS';...
    'onset time SD MS'; 'onset time SD V'; 'onset time SD M'; 'onset time SD NS'; ...
     'MSEnh num resp'; 'Uni bias num resp'};
 
% make ECDF plots for each var 
for i = 3:27
    s46 = allRespROIs(st46, i);
    s49 = allRespROIs(st49, i);
    figure;
    hold on
    ecdf(s46)

    %set(h, 'Color', 'g')
    ecdf(s49)
    h = get(gca, 'children')
    set(h, 'LineWidth', 3)
    %j = get(gca, 'children')
    %set(j, 'LineWidth', 3)
    %set(j, 'Color', [0.5 0 0.5])
    hold off
    title(sprintf('%s ECDF by stage', labels{i-2}))
    xlabel(labels{i-2})
    ylabel('ROI count')
    fig_filename = sprintf('st 46 vs 49 ECDF of %s', labels{i-2})
    saveas(gcf, fig_filename, 'png')
    close;
end


%% Does PCA differ by stage?

% this is a small subset of the parameters - only the ones that we think
% make sense to include. 

% copied from PCA_loadings
% This is based on Arseny's Elife paper plots
% looking for PCA loadings

% run PCA
% doc: https://www.mathworks.com/help/stats/pca.html

% get subset of allRespROIs data (eliminate area and peak loc)
allRespROIs_PCAsub = allRespROIs(:, [7:10, 15:27]);
labels_PCAsub = labels([5:8, 13:end])
% run PCA on all responding ROIs together
[coeff_all,score_all,latent_all,tsquared_all,explained_all,mu_all] = pca( [allRespROIs_PCAsub(:, :)]);
% run PCA on each stage separately
[coeff_46,score_46,latent_46,tsquared_46,explained_46,mu_46] = pca( [allRespROIs_PCAsub(st46, :)]);
[coeff_49,score_49,latent_49,tsquared_49,explained_49,mu_49] = pca( [allRespROIs_PCAsub(st49, :)]);

% the features (variables) put into allRespROIs(:, 3:27) are stored in labels

% plot coefficients and data onto PC1 x PC2 space
% https://www.mathworks.com/help/stats/biplot.html for documentation
figure;
biplot(coeff_all(:,1:2), 'Scores', score_all(:,1:2), 'VarLabels', labels_PCAsub)
title('PCA loadings all data')
figure;
biplot(coeff_46(:,1:2), 'Scores', score_46(:,1:2), 'VarLabels', labels_PCAsub)
title('PCA Loadings st 46 only')
figure;
biplot(coeff_49(:,1:2), 'Scores', score_49(:,1:2), 'VarLabels', labels_PCAsub)
title('PCA loadings st 49 only')

% plot amount of variance explained by each PC
figure;
bar(explained_all)
title('Variance explained by PC component all data')
xlabel('components')
ylabel('percent')
figure;
bar(explained_46)
title('Variance explained by PC component st46 only')
xlabel('components')
ylabel('percent')
figure;
bar(explained_49)
title('Variance explained by PC component st49 only')
xlabel('components')
ylabel('percent')

% plot MSEnh peak against PC1 value for all cells
figure;
scatter(score_all(:,1), allRespROIs_PCAsub(:,5))
title('MSEnh peak vs PC1 All Data')
xlabel('PC1')
ylabel('MSEnh peak')
saveas(gcf, 'MSEnh peak vs PC1 allRespROIs', 'png')

figure;
scatter(score_46(:,1), allRespROIs_PCAsub(st46,5))
title('MSEnh peak vs PC1 St46 only')
xlabel('PC1')
ylabel('MSEnh peak')
saveas(gcf, 'MSEnh peak vs PC1 st46 only', 'png')

figure;
scatter(score_49(:,1), allRespROIs_PCAsub(st49, 5))
title('MSEnh peak vs PC1 St49 only')
xlabel('PC1')
ylabel('MSEnh peak')
saveas(gcf, 'MSEnh peak vs PC1 st49 only', 'png')

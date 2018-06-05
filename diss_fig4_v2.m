%% diss_fig4: PCA

%% Revised PCA (based on 3/36 comments from CA)

% create new data pieces
% unisensory peak, onset, onset SD. And add overall response reliability

% Collapse unisensory information by averaging M and V 
for t = 1:length(allData) %tads
    for r = 1:size(allData{1,t}.avg_peak2, 2) %rois
        allData{1,t}.unimean_peak(r) = nanmean([allData{1,t}.avg_peak2(2,r), allData{1,t}.avg_peak2(3,r)]);
        allData{1,t}.unimean_onsettime(r) = nanmean([allData{1,t}.avg_onsettime2(2,r), allData{1,t}.avg_onsettime2(3,r)]);
        allData{1,t}.unistd_onsettime(r) = nanmean([allData{1,t}.std_onsettime2(2,r), allData{1,t}.std_onsettime2(3,r)]);
    end
end

% response reliability = boolean resp/sum resp
for t = 1:length(allData) %tads
    for r = 1:size(allData{1,t}.avg_peak2, 2) %rois
        allData{1,t}.resp_reliability(r) = allData{1,t}.sum_responses(r) / sum(sum(allData{1,t}.include(:,r,:)));
    end
end




%get list of all highcorr ROIs to ID for PCA
for t = 1:length(allData)
    if isfield(allData{1,t}, 'uniqueHighCorrROI_MS')
        if isfield(allData{1,t}, 'uniqueHighCorrROI_M')
           if isfield(allData{1,t}, 'uniqueHighCorrROI_V') 
                allData{1,t}.allhighcorrROIs = unique([allData{1,t}.uniqueHighCorrROI_MS; allData{1,t}.uniqueHighCorrROI_M; allData{1,t}.uniqueHighCorrROI_V]);
           else
               allData{1,t}.allhighcorrROIs = unique([allData{1,t}.uniqueHighCorrROI_MS; allData{1,t}.uniqueHighCorrROI_M]);
           end
        else 
            if isfield(allData{1,t}, 'uniqueHighCorrROI_V')
                allData{1,t}.allhighcorrROIs = unique([allData{1,t}.uniqueHighCorrROI_MS; allData{1,t}.uniqueHighCorrROI_V]);
            else
                allData{1,t}.allhighcorrROIs = allData{1,t}.uniqueHighCorrROI_MS;
            end
        end
    else
        if isfield(allData{1,t}, 'uniqueHighCorrROI_M')
           if isfield(allData{1,t}, 'uniqueHighCorrROI_V') 
                allData{1,t}.allhighcorrROIs = unique([allData{1,t}.uniqueHighCorrROI_M; allData{1,t}.uniqueHighCorrROI_V]);
           else
               allData{1,t}.allhighcorrROIs = unique([allData{1,t}.uniqueHighCorrROI_M]);
           end
        else 
            if isfield(allData{1,t}, 'uniqueHighCorrROI_V')
                allData{1,t}.allhighcorrROIs = unique([ allData{1,t}.uniqueHighCorrROI_V]);
            else
                allData{1,t}.allhighcorrROIs = [];
            end
        end
    end
end

% combine data points
idx = 1
for t = 1:length(allData)
    for roi = 1:length(allData{1,t}.resp_ROIs)
        r = allData{1,t}.resp_ROIs(roi);
        allRespROIs(idx, 1) = t;
        allRespROIs(idx, 2) = r;
        allRespROIs(idx, 3) = allData{1,t}.stage;
        allRespROIs(idx, 4) = allData{1,t}.unimean_peak(r);
        allRespROIs(idx, 5) = allData{1,t}.unimean_onsettime(r);
        allRespROIs(idx, 6) = allData{1,t}.unistd_onsettime(r); 
        allRespROIs(idx, 7) = allData{1,t}.resp_reliability(r); 
        allRespROIs(idx, 8) = allData{1,t}.avg_peak2(1,r); %multi peak
        allRespROIs(idx, 9) = allData{1,t}.avg_onsettime2(1,r); %multi avg onsettime
        allRespROIs(idx, 10) = allData{1,t}.std_onsettime2(1,r); %multi std onset time
        allRespROIs(idx, 11) = allData{1,t}.MSInd_peak2(r);
        allRespROIs(idx, 12) = allData{1,t}.MSInd_onsettime2(r);
        allRespROIs(idx, 13) = ismember(r,  allData{1,t}.allhighcorrROIs);
        idx = idx + 1
    end
end

% range normalize all data to [0,1]
for v = 4:12 % over all vars, skip cell info (1:3)
    range = nanmax(allRespROIs(:,v)) - nanmin(allRespROIs(:,v));
    allRespROIs_PCAsubRN(:,(v-3)) = (allRespROIs(:,v) - nanmin(allRespROIs(:,v))) / range;
end
% check to make sure that worked
largest = max(allRespROIs_PCAsubRN)
smallest = min(allRespROIs_PCAsubRN)



%% Run PCA on RN data
labels_PCAsub = {'Uni mean peak', 'Uni mean onset time', 'Uni std onset time', 'Response reliability', 'Multi mean peak', 'Multi mean onset time', 'Multi std onset time', 'MSInd peak', 'MSInd onset time'}
% run PCA on all responding ROIs together - range normalized
[coeff_all,score_all,latent_all,tsquared_all,explained_all,mu_all] = pca( [allRespROIs_PCAsubRN(:, :)]);

% plot coefficients and data onto PC1 x PC2 space
% https://www.mathworks.com/help/stats/biplot.html for documentation
% loading plot
figure;
biplot(coeff_all(:,1:2), 'VarLabels', labels_PCAsub)
%title('PCA loadings all data')
xlim([-0.1 1])
ylim([-0.4 0.8])
set(gca,'FontSize',30)
%saveas(gcf, 'PCA loadings all data fig3', 'png')
saveas(gcf, 'PCA loadings all data new 20180328', 'epsc2')

% C1 and C2 scores by stage
% stage is stored in allREspROIs_PCAsub, but PCA was done on
% allREspROIs_PCAsubRN
s46_forPCA = find(allRespROIs(:,3) == 46)
s49_forPCA = find(allRespROIs(:,3) == 49)
figure;
hold on 
plot(score_all(s46_forPCA, 1), score_all(s46_forPCA, 2), 'go')
plot(score_all(s49_forPCA, 1), score_all(s49_forPCA, 2), 'mo')
hold off
legend({'46', '49'})
xlabel('Component 1')
ylabel('Component 2')
set(gca, 'FontSize', 30)
saveas(gcf, 'PCA scores C1 C2 by stg new 20180328', 'epsc2')

% C1 and C2 scores colored by HC/SC
HC = find(allRespROIs(:,13) == 1)
SC = find(allRespROIs(:,13) == 0)
figure;
hold on 
plot(score_all(HC, 1), score_all(HC, 2), 'ro')
plot(score_all(SC, 1), score_all(SC, 2), 'ko')
hold off
legend({'HC', 'SC'})
xlabel('Component 1')
ylabel('Component 2')
set(gca, 'FontSize', 30)
saveas(gcf, 'PCA scores C1 C2 by corr new 20180328', 'epsc2')

% C1 and C2 scores colored by stage AND correlation 
HC_46 = intersect(HC, s46_forPCA)
SC_46 = intersect(SC, s46_forPCA)
HC_49 = intersect(HC, s49_forPCA)
SC_49 = intersect(SC, s49_forPCA)

figure;
hold on 
plot(score_all(HC_46, 1), score_all(HC_46, 2), 'bo')
plot(score_all(SC_46, 1), score_all(SC_46, 2), 'go')
plot(score_all(HC_49, 1), score_all(HC_49, 2), 'r*')
plot(score_all(SC_49, 1), score_all(SC_49, 2), '*', 'Color', [0.9290, 0.6940, 0.1250])
legend({'46 HC', '46 SC', '49 HC', '49 SC'})
hold off
xlabel('Component 1')
ylabel('Component 2')
set(gca, 'FontSize', 30)
saveas(gcf, 'PCA scores C1 C2 by stg and corr new 20180328', 'epsc2')

% stat test HC vs SC for C1 and C2
test = [1 2];
for i = 1:length(test)
    tmp = score_all(hc, test(i));
    dataHC = tmp(~isnan(tmp));
    tmp1 = score_all(nhc, test(i));
    dataNHC = tmp1(~isnan(tmp1));
    [HC_h(i), HC_p(i)] = kstest2(dataHC, dataNHC)
end
% both significant: 
%1: p = 0.048, 2: p = 0.024

%%%%%%%%% stat test 46 vs 49 fot C1 and C2
test = [1 2];
for i = 1:length(test)
    tmp = score_all(st46r, test(i));
    dataHC = tmp(~isnan(tmp));
    tmp1 = score_all(st49r, test(i));
    dataNHC = tmp1(~isnan(tmp1));
    [HC_h, HC_p] = kstest2(dataHC, dataNHC)
end
% 1 is not sig (p = 0.23). 2 is sig (p < 0.00001)

%% Plot distributions of C1 and C2 values split by HC/SC

filenames = {'PCA_HCvsSC C1', 'PCA_HCvsSC C2'}
Labels = {'C1 value', 'C2 value'}

% calculate median value
vars = [1, 2];
for i = 1:length(vars)
    hc_meds(i,1) = nanmedian(score_all(hc,vars(i)));
    hc_meds(i,2) = nanmedian(score_all(nhc,vars(i)));
end


list = [1, 2];
limits = [-0.4, 0.88; -0.65, 0.7];
%limits = [0 3];
for i = 1:length(list)
        figure;
        hold on
        histogram(score_all(hc, list(i)), 60, 'FaceColor', [0/256 73/256 73/256], 'EdgeColor', [0/256 73/256 73/256]) % hc in blue
        plot([hc_meds(i,1) hc_meds(i,1)], ylim, 'k', 'LineWidth', 3)
        hold off
        ylabel('ROI Count')
        xlabel(Labels{i})
        xlim(limits(i, :))
        set(gca,'FontSize',30)
        fig_filename = sprintf('hist hc_nhc HC %s', filenames{i})
        saveas(gcf, fig_filename, 'png')
        saveas(gcf, fig_filename, 'epsc2')
        %close;
        figure;
        hold on
        histogram(score_all(nhc, list(i)), 60, 'FaceColor', [146/256 73/256 0/256], 'EdgeColor', [146/256 73/256 0/256]) %nhc in purple
        plot([hc_meds(i,2) hc_meds(i,2)], ylim, 'k', 'LineWidth', 3)
        hold off
        ylabel('ROI Count')
        xlabel(Labels{i})
        xlim(limits(i, :))
        set(gca,'FontSize',30)
        fig_filename = sprintf('hist hc_nhc NHC %s', filenames{i})
        saveas(gcf, fig_filename, 'png')
        saveas(gcf, fig_filename, 'epsc2')
        %close;
        
end

for i = 1:length(list)
        [f1, x1] = ecdf(score_all(hc, list(i))); 
        [f2, x2] = ecdf(score_all(nhc, list(i))); 
        figure;
        hold on
        plot(x1, f1, 'Color', [0/256 73/256 73/256], 'LineWidth', 3) % hc in blue
        plot(x2, f2, 'Color', [146/256 73/256 0/256], 'LineWidth', 3) %nhc in purple
        hold off
        ylabel('ROI proportion')
        xlabel(Labels{i})
        xlim(limits(i, :))
        set(gca,'FontSize',30)
        fig_filename = sprintf('ecdf hc_nhc NHC %s', filenames{i})
        saveas(gcf, fig_filename, 'png')
        saveas(gcf, fig_filename, 'epsc2')
end

% figure out x limits
       % hcm = max(score_all(hc, list(i)))
        %nhcm = max(score_all(nhc, list(i)))


%% Plot distributions of C1 and C2 values split by stage 

filenames = {'PCA_46vs49 C1', 'PCA_46vs49 C2'}
Labels = {'C1 value', 'C2 value'}

% calculate median value
vars = [1, 2];
for i = 1:length(vars)
    hc_meds(i,1) = nanmedian(score_all(st46r,vars(i)));
    hc_meds(i,2) = nanmedian(score_all(st49r,vars(i)));
end


list = [1, 2];
limits = [-0.4, 0.88; -0.65, 0.7];
%limits = [0 3];
for i = 1:length(list)
        figure;
        hold on
        histogram(score_all(st46r, list(i)), 60, 'FaceColor', [255/256 0/256 255/256], 'EdgeColor', [255/256 0/256 255/256]) 
        plot([hc_meds(i,1) hc_meds(i,1)], ylim, 'k', 'LineWidth', 3)
        hold off
        ylabel('ROI Count')
        xlabel(Labels{i})
        xlim(limits(i, :))
        set(gca,'FontSize',30)
        fig_filename = sprintf('hist st46 %s', filenames{i})
        saveas(gcf, fig_filename, 'png')
        saveas(gcf, fig_filename, 'epsc2')
        %close;
        figure;
        hold on
        histogram(score_all(st49r, list(i)), 60, 'FaceColor', [0/256 128/256 0/256], 'EdgeColor', [0/256 128/256 0/256]) 
        plot([hc_meds(i,2) hc_meds(i,2)], ylim, 'k', 'LineWidth', 3)
        hold off
        ylabel('ROI Count')
        xlabel(Labels{i})
        xlim(limits(i, :))
        set(gca,'FontSize',30)
        fig_filename = sprintf('hist st49 %s', filenames{i})
        saveas(gcf, fig_filename, 'png')
        saveas(gcf, fig_filename, 'epsc2')
        %close;
        
end

for i = 1:length(list)
        [f1, x1] = ecdf(score_all(st46r, list(i))); 
        [f2, x2] = ecdf(score_all(st49r, list(i))); 
        figure;
        hold on
        plot(x1, f1, 'Color', [255/256 0/256 255/256], 'LineWidth', 3)
        plot(x2, f2, 'Color', [0/256 128/256 0/256], 'LineWidth', 3) 
        hold off
        ylabel('ROI proportion')
        xlabel(Labels{i})
        xlim(limits(i, :))
        set(gca,'FontSize',30)
        fig_filename = sprintf('ecdf %s', filenames{i})
        saveas(gcf, fig_filename, 'png')
        saveas(gcf, fig_filename, 'epsc2')
end

% figure out x limits
       % hcm = max(score_all(hc, list(i)))
        %nhcm = max(score_all(nhc, list(i)))
        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Old stuff 
%% 3C PCA
%%% Subset to use from allRespROIs: [7:9, 11, 13, 18:20, 22:24]
%%%% Removed 2 outlier cells with VERY LARGE MS indexes: Cells 41 and 66 
%%%% and VERY SMALL MS indexes: cells 235 and 1002 in
%%%% allRespROIs_PCAsub. Just deleted the rows so the matrix is smaller.

% Does PCA differ by stage?
% copied from PCA_loadings
% This is based on Arseny's Elife paper plots
% looking for PCA loadings

% run PCA
% doc: https://www.mathworks.com/help/stats/pca.html

% get subset of allRespROIs data (eliminate area and peak loc)
allRespROIs_PCAsub = allRespROIs([1:40, 42:65, 67:234, 236:1001, 1003:end], [7:9, 11, 13, 18:20, 22:24, 28]);
st46_noout = find(allRespROIs_PCAsub(:,12) == 46);
st49_noout = find(allRespROIs_PCAsub(:,12) == 49);
labels_PCAsub = labels([5:7, 9, 11, 16:18, 20:22])
% range normalize all data to [0,1]
for v = 1:(size(allRespROIs_PCAsub, 2)-1) % over all vars, skip stage
    range = max(allRespROIs_PCAsub(:,v)) - min(allRespROIs_PCAsub(:,v));
    allRespROIs_PCAsubRN(:,v) = (allRespROIs_PCAsub(:,v) - min(allRespROIs_PCAsub(:,v))) / range;
end
% check to make sure that worked
largest = max(allRespROIs_PCAsubRN)
smallest = min(allRespROIs_PCAsubRN)
% YES!

% run PCA on all responding ROIs together - range normalized
[coeff_all,score_all,latent_all,tsquared_all,explained_all,mu_all] = pca( [allRespROIs_PCAsubRN(:, :)]);
% run PCA on each stage separately
[coeff_46,score_46,latent_46,tsquared_46,explained_46,mu_46] = pca( allRespROIs_PCAsubRN(st46_noout, :));
[coeff_49,score_49,latent_49,tsquared_49,explained_49,mu_49] = pca( allRespROIs_PCAsubRN(st49_noout, :));

% the features (variables) put into allRespROIs(:, 3:27) are stored in labels

% plot coefficients and data onto PC1 x PC2 space
% https://www.mathworks.com/help/stats/biplot.html for documentation
% loading plot
figure;
biplot(coeff_all(:,1:2), 'VarLabels', labels_PCAsub)
%title('PCA loadings all data')
xlim([-0.05 1])
ylim([-0.5 0.8])
set(gca,'FontSize',30)
saveas(gcf, 'PCA loadings all data fig3', 'png')
saveas(gcf, 'PCA loadings all data fig3', 'epsc2')

% C1 and C2 scores by stage
% stage is stored in allREspROIs_PCAsub, but PCA was done on
% allREspROIs_PCAsubRN
s46_forPCA = find(allRespROIs_PCAsub(:,12) == 46)
s49_forPCA = find(allRespROIs_PCAsub(:,12) == 49)
figure;
hold on 
plot(score_all(s46_forPCA, 1), score_all(s46_forPCA, 2), 'go')
plot(score_all(s49_forPCA, 1), score_all(s49_forPCA, 2), 'mo')
hold off
xlabel('Component 1')
ylabel('Component 2')
set(gca, 'FontSize', 30)
saveas(gcf, 'PCA scores C1 C2 by stg', 'epsc2')


%% test for outliers

for v = 1:(size(allRespROIs_PCAsub, 2)-1) % for each variable included in PCA
    [PCA_outliers(:,v), PCA_outliers_stat(1,v), PCA_outliers_stat(2,v), PCA_outliers_stat(3,v)] = isoutlier(allRespROIs_PCAsub(:,v));
end



%% PCA on original data (outliers in)

allRespROIs_PCAsub_all = allRespROIs(:, [7:9, 11, 13, 18:20, 22:24, 28]);
for v = 1:(size(allRespROIs_PCAsub_all, 2)-1) % over all vars, skip stage
    range = max(allRespROIs_PCAsub_all(:,v)) - min(allRespROIs_PCAsub_all(:,v));
    allRespROIs_PCAsubRN_all(:,v) = (allRespROIs_PCAsub_all(:,v) - min(allRespROIs_PCAsub_all(:,v))) / range;
end

[coeff_all_all,score_all_all,latent_all_all,tsquared_all_all,explained_all_all,mu_all_all] = pca( [allRespROIs_PCAsubRN_all(:, :)]);

% loading plot
figure;
biplot(coeff_all_all(:,1:2), 'VarLabels', labels_PCAsub)
%title('PCA loadings all data')
xlim([-0.05 1])
ylim([-0.5 0.8])
set(gca,'FontSize',30)
saveas(gcf, 'PCA loadings all data fig3', 'png')
saveas(gcf, 'PCA loadings all data fig3', 'epsc2')

% C1 and C2 scores by stage
% stage is stored in allREspROIs_PCAsub, but PCA was done on
% allREspROIs_PCAsubRN
s46_forPCA_all = find(allRespROIs_PCAsub_all(:,12) == 46)
s49_forPCA_all = find(allRespROIs_PCAsub_all(:,12) == 49)
figure;
hold on 
plot(score_all_all(s46_forPCA_all, 1), score_all_all(s46_forPCA_all, 2), 'go')
plot(score_all_all(s49_forPCA_all, 1), score_all_all(s49_forPCA_all, 2), 'mo')
hold off
xlabel('Component 1')
ylabel('Component 2')
set(gca, 'FontSize', 30)
saveas(gcf, 'PCA scores C1 C2 by stg', 'epsc2')



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
biplot(coeff_46(:,1:2), 'Scores', score_46(:,1:2), 'VarLabels', labels_PCAsub)
title('PCA Loadings st 46 only')
set(gca,'FontSize',20)
saveas(gcf, 'PCA loadings 46 only fig3', 'png')
saveas(gcf, 'PCA loadings 46 only fig3', 'epsc2')
figure;
biplot(coeff_49(:,1:2), 'Scores', score_49(:,1:2), 'VarLabels', labels_PCAsub)
title('PCA loadings st 49 only')
set(gca,'FontSize',20)
saveas(gcf, 'PCA loadings 49 only fig3', 'png')
saveas(gcf, 'PCA loadings 49 only fig3', 'epsc2')

% plot amount of variance explained by each PC
figure;
bar(explained_all)
title('Variance explained by PC component all data')
xlabel('components')
ylabel('percent')
set(gca,'FontSize',20)
saveas(gcf, 'Variance explained by PC component all data fig3', 'png')
saveas(gcf, 'Variance explained by PC component all data fig3', 'epsc2')
figure;
bar(explained_46)
title('Variance explained by PC component st46 only')
xlabel('components')
ylabel('percent')
set(gca,'FontSize',20)
saveas(gcf, 'Variance explained by PC component 46 only fig3', 'png')
saveas(gcf, 'Variance explained by PC component 46 only fig3', 'epsc2')
figure;
bar(explained_49)
title('Variance explained by PC component st49 only')
xlabel('components')
ylabel('percent')
set(gca,'FontSize',20)
saveas(gcf, 'Variance explained by PC component 49 only fig3', 'png')
saveas(gcf, 'Variance explained by PC component 49 only fig3', 'epsc2')

% plot onset time MS against PC1 value for all cells [9, 10, 11]
figure;
scatter(score_all(:,1), allRespROIs_PCAsubRN(:,9))
title('SD of MS onset time vs PC1 All Data')
xlabel('PC1')
ylabel('SD of MS onset time')
set(gca,'FontSize',20)
saveas(gcf, 'SD of MS onset time vs PC1 allRespROIs', 'png')
saveas(gcf, 'SD of MS onset time vs PC1 allRespROIs', 'epsc2')
figure;
scatter(score_46(:,1), allRespROIs_PCAsub(st46_noout,9))
title('SD of MS onset time vs PC1 St46 only')
xlabel('PC1')
ylabel('SD of MS onset time')
set(gca,'FontSize',20)
saveas(gcf, 'SD of MS onset time vs PC1 st46 only', 'png')
saveas(gcf, 'SD of MS onset time vs PC1 st46 only', 'epsc2')
figure;
scatter(score_49(:,1), allRespROIs_PCAsub(st49_noout, 9))
title('SD of MS onset time vs PC1 St49 only')
xlabel('PC1')
ylabel('SD of MS onset time')
set(gca,'FontSize',20)
saveas(gcf, 'SD of MS onset time vs PC1 st49 only', 'png')
saveas(gcf, 'SD of MS onset time vs PC1 st49 only', 'epsc2')

% plot onset time M against PC1 value for all cells [9, 10, 11]
figure;
scatter(score_all(:,1), allRespROIs_PCAsubRN(:,10))
title('SD of M onset time vs PC1 All Data')
xlabel('PC1')
ylabel('SD of M onset time')
set(gca,'FontSize',20)
saveas(gcf, 'SD of M onset time vs PC1 allRespROIs', 'png')
saveas(gcf, 'SD of M onset time vs PC1 allRespROIs', 'epsc2')
figure;
scatter(score_46(:,1), allRespROIs_PCAsub(st46_noout,10))
title('SD of M onset time vs PC1 St46 only')
xlabel('PC1')
ylabel('SD of M onset time')
set(gca,'FontSize',20)
saveas(gcf, 'SD of M onset time vs PC1 st46 only', 'png')
saveas(gcf, 'SD of M onset time vs PC1 st46 only', 'epsc2')
figure;
scatter(score_49(:,1), allRespROIs_PCAsub(st49_noout, 10))
title('SD of M onset time vs PC1 St49 only')
xlabel('PC1')
ylabel('SD of M onset time')
set(gca,'FontSize',20)
saveas(gcf, 'SD of M onset time vs PC1 st49 only', 'png')
saveas(gcf, 'SD of M onset time vs PC1 st49 only', 'epsc2')

% plot onset time V against PC1 value for all cells [9, 10, 11]
figure;
scatter(score_all(:,1), allRespROIs_PCAsubRN(:,11))
title('SD of V onset time vs PC1 All Data')
xlabel('PC1')
ylabel('SD of V onset time')
set(gca,'FontSize',20)
saveas(gcf, 'SD of V onset time vs PC1 allRespROIs', 'png')
saveas(gcf, 'SD of V onset time vs PC1 allRespROIs', 'epsc2')
figure;
scatter(score_46(:,1), allRespROIs_PCAsub(st46_noout,11))
title('SD of V onset time vs PC1 St46 only')
xlabel('PC1')
ylabel('SD of V onset time')
set(gca,'FontSize',20)
saveas(gcf, 'SD of V onset time vs PC1 st46 only', 'png')
saveas(gcf, 'SD of V onset time vs PC1 st46 only', 'epsc2')
figure;
scatter(score_49(:,1), allRespROIs_PCAsub(st49_noout, 11))
title('SD of V onset time vs PC1 St49 only')
xlabel('PC1')
ylabel('SD of V onset time')
set(gca,'FontSize',20)
saveas(gcf, 'SD of V onset time vs PC1 st49 only', 'png')
saveas(gcf, 'SD of V onset time vs PC1 st49 only', 'epsc2')

% Plot onset time MS vs PC1 with 46 and 49 in different colors
figure;
hold on
scatter(score_46(:,1), allRespROIs_PCAsub(st46_noout,9), 'g')
scatter(score_49(:,1), allRespROIs_PCAsub(st49_noout, 9), 'm')
hold off
xlabel('PC1')
ylabel('SD of onset time')
set(gca,'FontSize',20)
saveas(gcf, 'SD of MS onset time vs PC1 alldata color spl', 'png')
saveas(gcf, 'SD of MS onset time vs PC1 alldata color spl', 'epsc2')
% onset time V vs PC1
figure;
hold on
scatter(score_46(:,1), allRespROIs_PCAsub(st46_noout,10), 'g')
scatter(score_49(:,1), allRespROIs_PCAsub(st49_noout, 10), 'm')
hold off
xlabel('PC1')
ylabel('SD of onset time')
set(gca,'FontSize',20)
saveas(gcf, 'SD of V onset time vs PC1 alldata color spl', 'png')
saveas(gcf, 'SD of V onset time vs PC1 alldata color spl', 'epsc2')
% onset time V vs PC1
figure;
hold on
scatter(score_46(:,1), allRespROIs_PCAsub(st46_noout,11), 'g')
scatter(score_49(:,1), allRespROIs_PCAsub(st49_noout, 11), 'm')
hold off
xlabel('PC1')
ylabel('SD of onset time')
set(gca,'FontSize',20)
saveas(gcf, 'SD of M onset time vs PC1 alldata color spl', 'png')
saveas(gcf, 'SD of M onset time vs PC1 alldata color spl', 'epsc2')


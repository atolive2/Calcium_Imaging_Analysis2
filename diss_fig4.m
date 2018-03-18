%% diss_fig4: PCA

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
figure;
biplot(coeff_all(:,1:2), 'Scores', score_all(:,1:2), 'VarLabels', labels_PCAsub)
title('PCA loadings all data')
set(gca,'FontSize',20)
saveas(gcf, 'PCA loadings all data fig3', 'png')
saveas(gcf, 'PCA loadings all data fig3', 'epsc2')
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


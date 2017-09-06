%% Principle Component Analysis new version
% This is based on Arseny's Elife paper plots
% looking for PCA loadings

% run PCA
% doc: https://www.mathworks.com/help/stats/pca.html
[coeff,score,latent,tsquared,explained,mu] = pca( [cluster_data_respROI(:, 3:21)])

% these are the features I put in cluster_data_respROI cols 3:end
labels = {'peak MS'; 'peak V'; 'peak M'; 'peak NS'; 'peak SD MS'; 'peak SD V'; 'peak SD M'; 'peak SD NS'; 'onset time MS'; 'onset time V'; 'onset time M'; 'onset time NS';...
    'onset time SD MS'; 'onset time SD V'; 'onset time SD M'; 'onset time SD NS'; ...
    'MSEnh peak'; 'MSEnh num resp'; 'Uni bias num resp'}

% plot coefficients and data onto PC1 x PC2 space
% https://www.mathworks.com/help/stats/biplot.html for documentation
figure;
biplot(coeff(:,1:2), 'Scores', score(:,1:2), 'VarLabels', labels)

% plot amount of variance explained by each PC
figure;
bar(explained)
title('Variance explained by PC component')
xlabel('components')
ylabel('percent')

% plot MSEnh peak against PC1 value for all cells
figure;
scatter(score(:,1), cluster_data_respROI(:,19))
title('MSEnh peak vs PC1')
xlabel('PC1')
ylabel('MSEnh peak')
saveas(gcf, 'MSEnh peak vs PC1', 'png')

%plot(latent)
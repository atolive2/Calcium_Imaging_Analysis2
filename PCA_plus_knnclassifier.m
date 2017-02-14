%% PCA for data with trials of length 160 only.
% this is to used the filtered data, to use the whole thing to find most
% useful PCs. 

figure_filepath = 'F:/Calcium_Imaging_Analysis/classifiers/figures/'
for t = 1:length(tadpole)
    % get only trials with 160 data points
    data_allcell = [];
    cols_incl = [];
    for i = 1:size(tadpole{1,t}.filtered,2) %number of trials
        if size(tadpole{1,t}.filtered{1,i})  == [160 1]
            for k = 1:size(tadpole{1,t}.filtered,1) %number of ROIs
                cell_data(k,:) = tadpole{1,t}.filtered{k,i};
            end
            data_allcell = cat(3, data_allcell, cell_data);
            cols_incl=[cols_incl i];
        end
        clear('cell_data')
    end
    trial_type = tadpole{1,t}.stimorder(1,cols_incl);

% DIMENSIONALITY REDUCTION (PCA)
    for c = 1:size(data_allcell,1)
        data_onecell(:,:) = data_allcell(c,:,:)
        % use principal component analysis
        [coeff, score] = pca((data_onecell'));
        % Isolate the projections from the first 3 PCs
        tadpole{1,t}.data_onecell3D(c,:,:) = score(:,1:3);
        data_onecell3D = score(:,1:3);
        % generate 1 fig per cell per tad
        figure
        f1 = find(trial_type == 1 | trial_type == 5 | trial_type == 8 | trial_type == 9 | trial_type ==10 | trial_type ==11)
        h = plot3(data_onecell3D(f1,1),data_onecell3D(f1,2),data_onecell3D(f1,3),'ro'); set(h,'MarkerFaceColor','r')
        xlabel('PC1')
        ylabel('PC2')
        zlabel('PC3')
        hold on
        f2 = find(trial_type == 2 | trial_type == 6);
        h = plot3(data_onecell3D(f2,1),data_onecell3D(f2,2),data_onecell3D(f2,3),'bo'); set(h,'MarkerFaceColor','b')
        f3 = find(trial_type == 3 | trial_type == 7 | trial_type == 12);
        h = plot3(data_onecell3D(f3,1),data_onecell3D(f3,2),data_onecell3D(f3,3),'ko'); set(h,'MarkerFaceColor','k')
        f4 = find(trial_type == 4 | trial_type == 12);
        h = plot3(data_onecell3D(f4,1),data_onecell3D(f4,2),data_onecell3D(f4,3),'ko'); set(h,'MarkerFaceColor','g')
        legend([{'multi'} {'vis'} {'mech'} {'no stim'}])
        fig_filename=sprintf([figure_filepath 'filteredexp%dcell%d.png'], t, c)
        saveas(gcf,fig_filename,'png');
        close;
        clear('fig_filename')
        clear('data_onecell', 'data_onecell3D', 'coeff', 'score', 'f1', 'f2', 'f3', 'f4')
    end
end

%% Run k-NN classifer on each cell %%%%%%%%%%%DOESN'T WORK%%%%%%%%%%%%%%
% Using filtered data, run classifier
% data is in tadpole{1,t}.data_onecell3D(c,:,:) = score(:,1:3);

for t = 1:length(tadpole)
    sprintf('start tadpole %d', t)
    % get only trials with 160 data points
    cols_incl = [];
    for i = 1:size(tadpole{1,t}.filtered,2) %number of trials
        if size(tadpole{1,t}.filtered{1,i})  == [160 1]
           cols_incl=[cols_incl i];
        end
    end
    tadpole{1,t}.trial_type = tadpole{1,t}.stimorder(1,cols_incl);
    trial_type = tadpole{1,t}.stimorder(1,cols_incl);
    if length(trial_type) > 1
        % run classifier on 1 cell at a time
        for c = 1:size(tadpole{1,t}.data_onecell3D ,1)
            data(:,:)=tadpole{1,t}.data_onecell3D(c,:,:);
            cellID = sprintf('cell_%d',c)
            tadpole{1,t}.KNNClass.(cellID) =fitcknn(data,trial_type);
            tadpole{1,t}.CVKNN.(cellID) = crossval(tadpole{1,t}.KNNClass.(cellID),'kfold',10);

            % and report the mean and sdev for classification accuracy
            tadpole{1,t}.klossKNN.(cellID) = kfoldLoss(tadpole{1,t}.CVKNN.(cellID),'mode','individual');
            tadpole{1,t}.accuracyM.(cellID) = mean(1-tadpole{1,t}.klossKNN.(cellID)) *100;
            tadpole{1,t}.accuracySTD.(cellID) = std(1-tadpole{1,t}.klossKNN.(cellID)) *100;

            % What should we expect by chance? 
            % Lets' shuffle the trial labels using 'randperm'
            % using 1000 trials - slow but correct.
            for r = 1:100
                tadpole{1,t}.CHKNNClass.(cellID) =fitcknn(data,trial_type(randperm(numel(trial_type))));
                tadpole{1,t}.CHRCVKNN.(cellID) = crossval(tadpole{1,t}.CHKNNClass.(cellID),'kfold',10);
                tadpole{1,t}.CHklossKNN.(cellID) = kfoldLoss(tadpole{1,t}.CHRCVKNN.(cellID),'mode','individual');
                tadpole{1,t}.CHaccuracyM.(cellID)(r,1) = mean(1-tadpole{1,t}.CHklossKNN.(cellID)) *100;
                tadpole{1,t}.CHaccuracySTD.(cellID)(r,1) = std(1-tadpole{1,t}.CHklossKNN.(cellID)) *100;
                %sprintf(['randperm %d complete'], r)
            end
            clear('data', 'cellID')
        end
    end
    clear('trial_type')
    sprintf('end tadpole %d', t)
end

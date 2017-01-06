%% Plots for PIE paper

% get data
myFolder = 'F:/Calcium_Imaging_Analysis/tadpoles_byexp/'; % May need to correct this.
mkdir([myFolder 'figures']);
if ~isdir(myFolder)
	errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
	uiwait(warndlg(errorMessage));
	return;
end
filePattern = fullfile(myFolder, 'exp*.mat');
matFiles = dir(filePattern)

for k = 1:length(matFiles)
	matFilename = fullfile(myFolder, matFiles(k).name)
	matData = load(matFilename); % Retrieves a structure.
%     for i = 1:length(matData.tadpole.somaticROICenters)
%         tmp = matData.tadpole.somaticROICenters{1,i}.Centroid;
%         ROIlocs(i,:) = [tmp(1,1), tmp(1,2)+(2*(250-tmp(1,2)))];
%     end
    
    % generate scatterplot of ROIs
    % color = primary modality
    % size = size of largest response

    % locations are stored in tadpole.somaticROICenters{1,:}.Centroid
    for i = 1:length(matData.tadpole.somaticROICenters)
        tmp = matData.tadpole.somaticROICenters{1,i}.Centroid;
        ROIlocs(i,:) = [tmp(1,1), tmp(1,2)+(2*(250-tmp(1,2)))];
    end

    %size
    sz = max(matData.tadpole.peak_avg(1:4,:))*500;

    %color
    % vis = blue, mech = red, MS = purple, no stim = white (e.g. no show)
    for i = 1:size(matData.tadpole.peak_avg,2)
        [m(i), primary_modality(i)] = max(matData.tadpole.peak_avg(1:4,i));
    end

    % change numbers into color
    for i = 1:length(primary_modality)
        if primary_modality(i) == 1
            color(i, :) = [0.5 0 0.5];
        elseif primary_modality(i) == 2
            color(i,:) = [0 0 1];
        elseif primary_modality(i) == 3
            color(i,:) = [1 0 0];
        elseif primary_modality(i) == 4
            color(i,:) = [1 1 1];
        else 
            color(i,:) = [0 1 1];
        end
    end

    % make plot
    figure;
    scatter(ROIlocs(:,1), ROIlocs(:,2), sz, color, 'filled')
    clear('matData', 'ROIlocs', 'primary_modality', 'color', 'm', 'sz')
end
    
% example experiment is #5. 
matFilename = fullfile(myFolder, matFiles(7).name)
matData = load(matFilename);

% pie chart of primary modality. 
    for i = 1:size(matData.tadpole.peak_avg,2)
        [m(i), primary_modality(i)] = max(matData.tadpole.peak_avg(1:4,i));
    end
    stims = unique(primary_modality)
    for i = 1:length(stims)
        piechart_data(i) = length(find(primary_modality == stims(i)))
    end
    
    percentages = round(piechart_data/sum(piechart_data)*100)
    %stim_labels = {'multisensory', 'visual', 'mechanosensory', 'no response'}
    figure;
    txt = {'multisensory: 51%'; 'visual: 13%'; 'mechanosensory: 28%'; 'no response: 9%'}; % strings
    h = pie(piechart_data, txt);
    hp = findobj(h, 'Type', 'patch');
    set(hp(1), 'FaceColor', [0.5 0 0.5]);
    set(hp(2), 'FaceColor', [0 0 1]);
    set(hp(3), 'FaceColor', [1 0 0]);
    set(hp(4), 'FaceColor', [1 1 1]);
    
hText = findobj(h,'Type','text'); % text object handles
%percentValues = get(hText,'String'); % percent values

%combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(hText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell); % numeric array

%     hText(1).String = combinedtxt(1);
% hText(2).String = combinedtxt(2);
% hText(3).String = combinedtxt(3);
% hText(3).String = combinedtxt(4);

% %     hText(1).String = txt(1);
% % hText(2).String = txt(2);
% % hText(3).String = txt(3);
% % hText(4).String = txt(4);

    newExtents_cell = get(hText,'Extent'); % cell array
newExtents = cell2mat(newExtents_cell); % numeric array
width_change = newExtents(:,3)-oldExtents(:,3);
signValues = sign(oldExtents(:,1));
offset = signValues.*(width_change/2);
textPositions_cell = get(hText,{'Position'}); % cell array
textPositions = cell2mat(textPositions_cell); % numeric array
textPositions(:,1) = textPositions(:,1) + offset; % add offset

hText(1).Position = textPositions(1,:);
hText(2).Position = textPositions(2,:);
hText(3).Position = textPositions(3,:);
hText(4).Position = textPositions(4,:);


%% Plot all responses from 1 trialblock (first one) by response type
figure;
plot(cell2mat(matData.tadpole.df_f0(1:20,1))', 'color', [0.5 0 0.5])
figure;
plot(cell2mat(matData.tadpole.df_f0(1:20,22))', 'color', [0 0 1])
figure; plot(cell2mat(matData.tadpole.df_f0(1:20,3))', 'color', [1 0 0])
figure; plot(cell2mat(matData.tadpole.df_f0(1:20,4))', 'color', [0 0 0])

for k = 1:length(matData.tadpole.stimorder)
    figure;
    if matData.tadpole.stimorder(k) == 1
        plot(cell2mat(matData.tadpole.df_f0(:,k))', 'color', [0.5 0 0.5])
    elseif matData.tadpole.stimorder(k) == 2
        plot(cell2mat(matData.tadpole.df_f0(:,k))', 'color', [0 0 1])
    elseif matData.tadpole.stimorder(k) == 3
        plot(cell2mat(matData.tadpole.df_f0(:,k))', 'color', [1 0 0])
    elseif matData.tadpole.stimorder(k) == 4
        plot(cell2mat(matData.tadpole.df_f0(:,k))', 'color', [0 0 0])
    end
    fig_filename=sprintf('F:/Calcium_Imaging_Analysis/Ca_exp_%d/figures/exp%dtrial%d_all_%d',matData.tadpole.expnum,matData.tadpole.expnum, k, matData.tadpole.stimorder(k));
    saveas(gcf,fig_filename,'png');
    close
end

    
fig_filename=sprintf('F:/Calcium_Imaging_Analysis/Ca_exp_%d/figures/exp%dMSI_PIE_area2_5_unisplit.png',expnum,expnum);
    saveas(gcf,fig_filename,'png');
    
    
%% If you have already opened alltads_exps1-11_PIEpaper.mat

% this is the new graph for paper. Uses experiment 5, which is 7 in
% matFiles and tadpole. 
    % generate scatterplot of ROIs
    % color = primary modality 
    % size = size of largest response (make non responding ROIs size
    % consistent. 
for t = 1:length(tadpole)
    % locations are stored in tadpole.somaticROICenters{1,:}.Centroid
    for i = 1:length(tadpole{1,t}.somaticROICenters)
        tmp = tadpole{1,t}.somaticROICenters{1,i}.Centroid;
        ROIlocs(i,:) = [tmp(1,1), tmp(1,2)+(2*(250-tmp(1,2)))];
    end

    %size
    %sz = max(tadpole{1,5}.peak_avg(1:4,:))*500;

    %color
    % vis = blue, mech = red, MS = purple, no stim = black dot (e.g. minimize)
    for i = 1:size(tadpole{1,t}.peak_avg,2)
        [m(i), primary_modality(i)] = max(tadpole{1,8}.peak_avg(1:4,i));
    end

    % change numbers into color
    for i = 1:length(primary_modality)
        if primary_modality(i) == 1
            color(i, :) = [0.5 0 0.5];
            sz(i,:) = m(i)*500;
        elseif primary_modality(i) == 2
            color(i,:) = [0 0 1];
            sz(i,:) = m(i)*500;
        elseif primary_modality(i) == 3
            color(i,:) = [1 0 0];
            sz(i,:) = m(i)*500;
        elseif primary_modality(i) == 4
            color(i,:) = [0 0 0];
            sz(i,:) = 0.005*500;
        else 
            color(i,:) = [0 1 1];
            sz(i,:) = 0.006*500;
        end
    end

    % make plot
    figure;
    scatter(ROIlocs(:,1), ROIlocs(:,2), sz, color, 'filled')
    title(matFiles(t).name)
    clear('matData', 'ROIlocs', 'primary_modality', 'color', 'm', 'sz')
end

%% make pie charts of each tadpole
for t = 1:length(tadpole)
    for i = 1:size(tadpole{1,t}.peak_avg,2)
        [m(i), primary_modality(i)] = max(tadpole{1,8}.peak_avg(1:4,i));
    end
    props(1,t) = sum(primary_modality==1)/length(primary_modality) %ms
    props(2,t) = sum(primary_modality==2)/length(primary_modality) %vis
    props(3,t) = sum(primary_modality==3)/length(primary_modality) %mech
    props(4,t) = sum(primary_modality==4)/length(primary_modality) %no stim
    props(5,t) = sum(props(1:4,t))
    clear('m', 'primary_modality')
end

percentages = round((props(1:4,:)*100))

for t = 1:size(props,2)
    txt_ms = sprintf('multisensory: %d%%', percentages(1,t))
    txt_v = sprintf('visual: %d%%', percentages(2,t))
    txt_mech = sprintf('mechanosensory: %d%%', percentages(3,t))
    txt_ns = sprintf('no stimulus: %d%%', percentages(4,t))
    all_txt = {txt_ms, txt_v, txt_mech, txt_ns}
    figure;
    h = pie(props(1:4,t), all_txt);
    hp = findobj(h, 'Type', 'patch');
    set(hp(1), 'FaceColor', [0.5 0 0.5]);
    set(hp(2), 'FaceColor', 'r');
    set(hp(3), 'FaceColor', 'b');
    set(hp(4), 'FaceColor', 'w');
    title(matFiles(t).name)
    clear('fig_filename')
    fig_filename=sprintf('F:/Calcium_Imaging_Analysis/tadpoles_byexp/figures/sfn2016poster/exp%dpiechart.eps', tadpole{1,t}.expnum)
    saveas(gcf,fig_filename,'psc2');
end

avg_prop = mean(props, 2)
for i = 1:(size(props,1)-1)
    stdev_prop(i) = std(props(i,:))
end

figure;
hold on
bar(avg_prop(1:4,1), 0.6)
errorbar(avg_prop(1:4,1), stdev_prop, '.')
hold off


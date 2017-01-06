%% Proximity analysis version 3 - based on SFN feedback 11/2016
load('F:\Calcium_Imaging_Analysis\tadpoles_byexp\alltads_exps1-11_proximity.mat')

%%%%%%% This assumes you have opened workspace 'alltads_exps1-11_proximity.mat'

%% local clustering: 
% set a threshold distance from an ROI and calculate proportion of ROIs with 
% same and different primary modality within the radius

% Find and set the threshold

% find maximum distance from ROI to ROI
for t = 1:length(roiEucDist_byTadpole)
    %max_distance(1,t) = max(roiEucDist_byTadpole{1,t})
    %min_distance(1,t) = min(roiEucDist_byTadpole{1,t})
    median_distance(1,t) = median(roiEucDist_byTadpole{1,t})
end
avg_max_dist = mean(max_distance)
avg_min_dist = mean(min_distance)
avg_med_dist = mean(median_distance)

% threshold = some percentage of avg_med_dist
% lets start by testing 20% and 50% (and 100% as a control)
threshold01 = 0.01 * avg_med_dist
threshold05 = 0.05 * avg_med_dist
threshold20 = 0.2 * avg_med_dist
threshold50 = 0.5 * avg_med_dist
threshold100 = 1 * avg_med_dist

% establish which ROIs are within threshold distance of a given ROI. 

% use pairMap2d to locate ROI-ROI values
% this creates a matrix that is numrois by (numrois-1), with the rows
% indicating each pair location for the roi corresponding to that row. The
% cols, however, are right until the col num = row num (e.g. the roi to
% itself), then they are off by one because there is no calculation for roi
% to itself included in the pairMap. 

for t = 1:length(tadpole)
    numrois = max(max(pairMap2d{1,t}))
    for i = 1:numrois %over each roi
        [rowi, coli] = find(pairMap2d{1,t} == i); % coli gives value needed
        tadpole{1,t}.distlocs_byroi(i,:) = coli;
    end
end

% Add NaN to all diagonal entries (thereby putting all ROIs in the correct
% location (row and col) to identify them. 
for t = 1:length(tadpole)
    for i = 1:size(tadpole{1,t}.distlocs_byroi, 1) %over each roi
        if i == 1
        tadpole{1,t}.distlocs_byroi_corr(i,:)=[NaN, tadpole{1,t}.distlocs_byroi(i,1:end)]; 
        elseif i > 1
        idx=i-1;
        tadpole{1,t}.distlocs_byroi_corr(i,:) = [tadpole{1,t}.distlocs_byroi(i,1:idx), NaN, tadpole{1,t}.distlocs_byroi(i,(idx+1):end)];
        end
    end
end

% if you fuck up above
% for t = 1:length(tadpole)
%     tadpole{1,t}.distlocs_byroi = [];
% end

% use distlocs_byroi to get distances from roiEucDist_byTadpole sorted by
% roi
for t = 1:length(tadpole)
    for i = 1:size(tadpole{1,t}.distlocs_byroi,1)
        for j = 1:size(tadpole{1,t}.distlocs_byroi,2)
            if isnan(tadpole{1,t}.distlocs_byroi(i,j))
               tadpole{1,t}.dist_byroi(i,j) = NaN
            elseif ~isnan(tadpole{1,t}.distlocs_byroi(i,j))
            tadpole{1,t}.dist_byroi(i,j) = roiEucDist_byTadpole{1,t}(1, tadpole{1,t}.distlocs_byroi(i,j));
            end
        end
    end
end

% histogram of distances 
for t = 1:length(tadpole)
    figure;
    hist(tadpole{1,t}.dist_byroi(1,:),50)
    title(sprintf('tadpole %d', t))
end

% threshold distances and generate a logical type matrix to identify
% "nearby" ROIs for a given ROI.

for t = 1:length(tadpole)
    tadpole{1,t}.log_byroi50 = tadpole{1,t}.dist_byroi < threshold50
    tadpole{1,t}.log_byroi20 = tadpole{1,t}.dist_byroi < threshold20
    tadpole{1,t}.sumlog_byroi20=sum(log_byroi20,2)
    tadpole{1,t}.sumlog_byroi50=sum(log_byroi50,2)
end

% Assess proportion with each primary modality of "nearby" ROIs for each
% ROI. 
% primary modality is stored in tadpole{1,t}.peak_avg_primary_modality
% 20%
for t = 1:length(tadpole)
    for i = 1:size(tadpole{1,t}.log_byroi20,1) % over each roi
        [row, close_rois] = find(tadpole{1,t}.log_byroi20(i,:));
        prim_mod = tadpole{1,t}.peak_avg_primary_modality(1,close_rois);
        tadpole{1,t}.count_primarymodality_byroi20(i,1) = length(find(prim_mod == 1));
        tadpole{1,t}.count_primarymodality_byroi20(i,2) = length(find(prim_mod == 2));
        tadpole{1,t}.count_primarymodality_byroi20(i,3) = length(find(prim_mod == 3));
        tadpole{1,t}.count_primarymodality_byroi20(i,4) = length(find(prim_mod == 4));
        clear('row', 'close_rois', 'prim_mod')
    end
end

% 50%
for t = 1:length(tadpole)
    for i = 1:size(tadpole{1,t}.log_byroi50,1) % over each roi
        [row, close_rois] = find(tadpole{1,t}.log_byroi50(i,:));
        prim_mod = tadpole{1,t}.peak_avg_primary_modality(1,close_rois);
        tadpole{1,t}.count_primarymodality_byroi50(i,1) = length(find(prim_mod == 1));
        tadpole{1,t}.count_primarymodality_byroi50(i,2) = length(find(prim_mod == 2));
        tadpole{1,t}.count_primarymodality_byroi50(i,3) = length(find(prim_mod == 3));
        tadpole{1,t}.count_primarymodality_byroi50(i,4) = length(find(prim_mod == 4));
        clear('row', 'close_rois', 'prim_mod')
    end
end

% What proportion of nearby ROIs are same primary modality as given ROI?








% Do statistics. 



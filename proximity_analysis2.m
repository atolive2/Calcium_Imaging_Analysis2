%% Proximity analysis

% this code assumes you have opened alltads_exps1-11_proximity
% distance formula: sqrt( (x2 - x1) ^2 + (y2 - y1) ^2 )

% ROI center info is in tadpole{1,i}.somaticROICenters{1,1}.Centroid and is
% type double

%% calculate the distance from each ROI to every other ROI in a tadpole
% code assumes you have tadpole{1,:} struct in workspace.
% get euclidian distance (in pixels) for Torrey's Tadpoles. (Written by
% Chris Deister, 11/7/2016). 

pixelToDistanceScale=2.26; % this is the number of microns per pixel. I assume 1 for this code. N microns per pixel based on a calibration.

for n=1:numel(tadpole)
    pairNum=1;
    tadpole_totalPairs(:,n)=nchoosek(numel(tadpole{n}.somaticROICenters),2);
    for k=1:numel(tadpole{n}.somaticROICenters)
        for j=k+1:numel(tadpole{n}.somaticROICenters)

            roiDXs_byTadpole{n}(:,pairNum)=fix(abs(tadpole{n}.somaticROICenters{j}(1).Centroid(2)-tadpole{n}.somaticROICenters{k}(1).Centroid(2)));
            % faster but not perfect
            roiDYs_byTadpole{n}(:,pairNum)=fix(abs(tadpole{n}.somaticROICenters{j}(1).Centroid(1)-tadpole{n}.somaticROICenters{k}(1).Centroid(1)));
            % we round because only integers should be allowed because we
            % have not transformed from pixels yet; the interpolation done
            % by the roi centroid alogorithm can lead to decimals.
            pairMap{n}(:,:,pairNum)=[k,j];
            pairNum=pairNum+1;
        end
    end
    % assign distance as dy/dx, but handle cases where they are either on
    % the same x or y pixel.
    roiEucDist_byTadpole{n}=zeros(size(roiDXs_byTadpole{n}));
    dNzInds=find(roiDYs_byTadpole{n}~=0 & roiDXs_byTadpole{n}~=0);
    xZInds=find(roiDXs_byTadpole{n}==0);
    yZInds=find(roiDYs_byTadpole{n}==0);
    roiEucDist_byTadpole{n}(dNzInds)=roiDYs_byTadpole{n}(dNzInds)./roiDXs_byTadpole{n}(dNzInds);
    roiEucDist_byTadpole{n}(yZInds)=roiDXs_byTadpole{n}(yZInds);
    roiEucDist_byTadpole{n}(xZInds)=roiDYs_byTadpole{n}(xZInds);
    
end

% makes a histogram of all distances for a given tadpole.
figure,nhist(roiEucDist_byTadpole{7},'box','minbins',100,'proportion')

%% group cells by primary modality in each tadpole, ttest


% first identify primary modality (max of peak_avg for stim 1-4 (high MS,
% M, V and no stim)
for t = 1:length(tadpole)
    for i = 1:size(tadpole{1,t}.peak_avg,2) %over each ROI
        [tadpole{1,t}.peak_avg_max(i), tadpole{1,t}.peak_avg_primary_modality(i)] = max(tadpole{1,t}.peak_avg(1:4,i));
    end
end

% Divide rois of same primary modality and of different primary modailty
% for each roi.
t = 1
for t = 1:length(tadpole)
    for r = 1:length(tadpole{1,t}.peak_avg_primary_modality)
        tadpole{1,t}.same_prim_mod{r,:} = find(tadpole{1,t}.peak_avg_primary_modality(1,:) == tadpole{1,t}.peak_avg_primary_modality(1,r))
        tadpole{1,t}.diff_prim_mod{r,:} = find(tadpole{1,t}.peak_avg_primary_modality(1,:) ~= tadpole{1,t}.peak_avg_primary_modality(1,r))
    end
end

% Gather distances for the rois in .same_prim_mod and .diff_prim_mod

% change pairMap from 3D to 2D by collapsing along dim1 (size =1 anyway)
for t = 1:length(pairMap)
    pairMap2d{1,t}(:,:) = pairMap{1,t}(1,:,:)
end

% get logical array of locations in pairMap2d where the rois of interest
% are
for t = 1:length(tadpole)
    for r = 1:length(tadpole{1,t}.same_prim_mod)
        roimap_same{1,t}{r,1} = ismember(pairMap2d{1,t}(:,:), tadpole{1,t}.same_prim_mod{r,1}(1,:));
        roimap_diff{1,t}{r,1} = ismember(pairMap2d{1,t}(:,:), tadpole{1,t}.diff_prim_mod{r,1}(1,:));
    end
end

% Create roimap2 to consolidate roimap into 1 vector
for t = 1:length(tadpole)
    for r = 1:length(tadpole{1,t}.same_prim_mod)
        roimap2_same{1,t}{r,1} = roimap_same{1,t}{r,1}(1,:) | roimap_same{1,t}{r,1}(2,:);
        roimap2_diff{1,t}{r,1} = roimap_diff{1,t}{r,1}(1,:) | roimap_diff{1,t}{r,1}(2,:);
    end
end

% Get distances for each relevant roi
% same
for t = 1:length(tadpole)
    for r = 1:length(tadpole{1,t}.same_prim_mod)
        same_prim_mod_vals=[];
        for k = 1:length(roimap2_same{1,t}{r,1})
            if roimap2_same{1,t}{r,1}(1,k) 
                same_prim_mod_vals=[same_prim_mod_vals; roiEucDist_byTadpole{1,t}(1,k)];
            end
        end
        tadpole{1,t}.same_prim_mod_vals{r,1} = same_prim_mod_vals;
    end
end
% diff
for t = 1:length(tadpole)
    for r = 1:length(tadpole{1,t}.diff_prim_mod)
        diff_prim_mod_vals=[];
        for k = 1:length(roimap2_diff{1,t}{r,1})
            if roimap2_diff{1,t}{r,1}(1,k) 
                diff_prim_mod_vals=[diff_prim_mod_vals; roiEucDist_byTadpole{1,t}(1,k)];
            end
        end
        tadpole{1,t}.diff_prim_mod_vals{r,1} = diff_prim_mod_vals;
    end
end

% if you fuck up above
% for t = 1:length(tadpole)
%     tadpole{1,t}.diff_prim_mod_vals = {};
%     tadpole{1,t}.same_prim_mod_vals = {};
% end

% stats on each roi - is same_prim_mod_vals different from
% diff_prim_mod_vals?
for t = 1:length(tadpole)
    clear('ttest2_h', 'ttest2_p', 'ttest2_ci', 'ttest2_stats');
    for r = 1:length(tadpole{1,t}.same_prim_mod_vals)
        %[h(r,1), p(r,1), ci(r,1), stats(r,1)] = ttest2(tadpole{1,t}.same_prim_mod_vals{r,1}, tadpole{1,t}.diff_prim_mod_vals{r,1});
        [h, p, ci, stats] = ttest2(tadpole{1,t}.same_prim_mod_vals{r,1}, tadpole{1,t}.diff_prim_mod_vals{r,1});
        ttest2_h(r,1)=h;
        ttest2_p(r,1)=p;
        ttest2_ci{r,1}=ci;
        ttest2_stats(r,1)=stats;
    end
    tadpole{1,t}.h = ttest2_h;
    tadpole{1,t}.p = ttest2_p;
    tadpole{1,t}.ci = ttest2_ci;
    tadpole{1,t}.stats = ttest2_stats;
end

% how many cells and what proportion of ROIs in each tadpole are different?
% use h (type double, with 1 = reject null, 0 = don't reject null). 
% Null = null hypothesis that the data in vectors x and y comes from independent 
% random samples from normal distributions with equal means and equal but unknown variances
for t = 1:length(tadpole)
    sum_h(1,t) = sum(tadpole{1,t}.h)
    prop_h(1,t) = sum_h(1,t)/length(tadpole{1,t}.h)
end

% what's the primary modality of the different ROIs?
for t = 1:length(tadpole)
    compare_t2_1{1,t}(:,1) = tadpole{1,t}.h
    compare_t2_1{1,t}(:,2) = tadpole{1,t}.peak_avg_primary_modality
end


for t = 1:length(tadpole)
end

% stats on all rois in a tadpole

%% collapse across tadpoles (graph 1 point by modality)

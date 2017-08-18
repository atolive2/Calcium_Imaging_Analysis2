function pairwise_distances = get_pairwise_distances(ROIdist, roi_pairs, YrespROI, respROI_list)
% this function takes ROI pairs and gets the distances between them,
% assuming you have a matrix of pairwise distances and a 2 col matrix with
% the ROI pairs. 

% If YrespROI = 1, then this section executes to map ROIs back onto real
% ROIs. If YrespROI = 0, then this section is skipped and the ROI labels in
% roi_pairs is used directly. 
if YrespROI
    ROI_list = respROI_list(roi_pairs);
end
if size(ROI_list,2) == 1
    ROI_list = ROI_list';
end
if ~isempty(ROI_list)
    for i = 1:size(ROI_list, 1)
        dists_tmp(i) = ROIdist(ROI_list(i,1), ROI_list(i,2));
    end
    pairwise_distances =  [ROI_list dists_tmp';];
else
    pairwise_distances = [];
end
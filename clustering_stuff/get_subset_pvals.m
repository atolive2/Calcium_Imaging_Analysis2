function [rois] = get_subset_pvals(lower_bound, upper_bound, all_pvals)
% This function takes P values from a correlation matrix in all_pvals and finds ROIs
% within the range lower_bound:upper_bound. 
% for respondingROI based data, this is the ROI assignments made once
% responding ROIs were extracted and therefore a second step is necessary
% to get back to the actual ROIs. 
if sum(size(all_pvals)) > 0 
    [tmp_sm(:,1), tmp_sm(:,2)] = find(all_pvals >= lower_bound);
    [tmp_lg(:,1), tmp_lg(:,2)] = find(all_pvals < upper_bound);
    result = ismember(tmp_sm, tmp_lg, 'rows');
    tmp_rois = tmp_sm(result, :);
    % elimnate duplicates (each cell is accounted for twice in correlation
    % matrix, but only want it once when doing stats and plotting).
    tmp_dups = tmp_rois(:,1) < tmp_rois(:,2);
    rois = tmp_rois(tmp_dups, :);
else
    rois = [];
end
end


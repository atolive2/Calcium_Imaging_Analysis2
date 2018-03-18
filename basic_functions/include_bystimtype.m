function [ include ] = include_bystimtype( stimorder, trialmask )
%respct_bystimtype defines good trials by ROI and stimtype
%   stimorder = vector with ID of trial type
%   trialmask = same size as data, 1=include

stims = unique(stimorder)
stimmask = get_stimmask(stimorder)
for s = 1:length(stims)
    for r = 1:size(trialmask, 1) %roi
        include(s,r,:) = stimmask(:,s)' & trialmask(r,:);
        
    end
end

    


end


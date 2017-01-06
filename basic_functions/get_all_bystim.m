function [ stim_vals ] = get_all_bystim( stimmask, val_of_interest, data )
%get_all_bystim takes stimmask and returns all values associated with that
%stimmask in a vector
%   inputs: stimmask = a logical array of size num_trials by num stim types
%           data = an array of size num_trial_types by num_rois, 1 data per
%           stim per roi.
%   output: a cell array of matrices of values (num_stim_type by num_rois) that meet the stimulus condition 

%stimvals = cell{

if length(val_of_interest) > 1
    for i = 1:length(val_of_interest)
        stimmask_single = stimmask(:,i);
        idx=1;
        for j = 1:length(stimmask_single)
            if stimmask_single(j)==1
            stimvals{i,idx} = data(:, j);
            idx = idx+1;
            end
        end
    end
    
else length(val_of_interest == 1)
    stimmask_single = stimmask(:,val_of_interest(1,1));
    idx=1;
    for j = 1:length(stimmask_single)
        if stimmask_single(j)==1
        stimvals{1,idx} = data(:, j);
        idx = idx+1;
        end
    end
end

for i = 1:size(stimvals,1)
    cat_data=[];
    for j = 1:length(stimvals(i,:))
        to_add = stimvals{i,j};
        cat_data = cat(1,cat_data,to_add);
    end
    stim_vals{i}=cat_data;
end
end

%hist(cell2mat(stim_vals{1,1}), 100)
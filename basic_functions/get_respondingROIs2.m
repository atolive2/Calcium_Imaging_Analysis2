function [ boolean_response, sum_responses ] = get_respondingROIs( area, peak, peakloc )
%get_respondingROIs takes the area, peak and peakloc information for a set of trials and ROIS
% and reports 1 if response and 0 if no response in a boolean array with same
% dims as original. 
%   inputs: area, peak and peakloc are all cell arrays with the same dims
%   of num_rois by num_trials. Must input data for the same cells.
%   output: boolean_response = summary stat, if ROI in a trial meets all
%   criteria, a 1 is present.

%Criteria for elimination:
% -	Peak location is earlier than the stimulus onset
% -	Area is negative
% -	Peak is negative
% - peak is greater than 10

boolean_response = zeros(size(area));
sum_responses = zeros(size(area,1), 1);
for i = 1:size(area,1)
    for j = 1:size(area,2)
        if area{i,j} > 0 && peak(i,j) > 0.2 && peakloc(i,j) > 20 && peak(i,j) < 10
            boolean_response(i,j) = 1;
        end
    end
    sum_responses(i)=sum(boolean_response(i,:));
end

end
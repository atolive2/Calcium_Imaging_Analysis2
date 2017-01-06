function [ MSenhancement ] = calc_MSenhancement( data )
%calc_MSenhancement calculated the proportion enhancement in multisensory
%vs largest unisensory.
%   inputs: stimmask = a logical array of dims num_trials by 9 (trial
%   types)
%       data = an array of dims numtrialtypes(9) by num_rois, with 1 value per roi
%   per trial (e.g. area, peak). values near 1 = no enhancement, above 1 =
%   enhancement, below 1 = supression
%   output: a 4 x num_rois matrix. 
%           row 1 = high M/high V, 
%           row 2 = high M/low V,
%           row 3 = low M/high V,
%           row 4 = low M/low V

    %1 = multisensory high M / high V
    %2 = visual high/crash
    %3 = mechanosensory high
    %4 = no stimulus
    
    %5 = multisensory high M / low V
    %6 = visual low/scrambled crash
    %7 = mechanosensory low
    
    %8 = multisensory low M / high V
    %9 = multisensory low M / low V

    %10 = multisensory med M / high V
    %11 = multisensory med M / low V
    %12 = mechanosensoy med
    
    % this matrix aligns the 3 conditions - row = MS, V, M, 
    % col = MS type, V type that matches MS, M type that matches
if size(data,1) == 4
    MSstims = [1; 2; 3];
elseif size(data,1) == 6
    MSstims = [1, 5;
               2, 6;
               3, 3];
elseif size(data,1) == 9
    MSstims = [1, 5, 8, 9;
               2, 6, 2, 6;
               3, 3, 7, 7];
elseif size(data,1) == 12
        MSstims = [1, 5, 8, 9, 10;
                   2, 6, 2, 6, 6;
                   3, 3, 7, 7, 12];
else 
    sprintf('error: dim of data not 4, 8 or 12')
    return
end
       
for i = 1:size(MSstims, 2) % over each stim type
    for j = 1:size(data,2) % over each ROI
        MSenhancement(i,j) = data(MSstims(1,i), j) / max([data(MSstims(2,i),j), data(MSstims(3,i),j)]);
    end
end


end


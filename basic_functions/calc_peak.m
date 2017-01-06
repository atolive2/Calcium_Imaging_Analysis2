function [ meanpeak_bytrial, peakloc_bytrial, meanpeak_bytrial_errors ] = calc_peak( signal )
%calc_peak takes df/f0 and returns the location of the peak and its value 
%   input = data = df/f0, a cell array with dims num_ROIs by num_trials
%   output: peak_bytrial is a cell array with dims num_ROIs by num_trials,
%   col 1 is peak, col 2 is location
meanpeak_bytrial_errors = [];
for i = 1:size(signal,1)
    for j = 1:size(signal,2)
        [peak{i,j}, peakloc_bytrial{i,j}] = max(signal{i,j});
        lrange = peakloc_bytrial{i,j} - 1;
        urange = peakloc_bytrial{i,j} + 1;
        if lrange > 1 && urange > 0 && urange < length(signal{i,j}(1,:))
            meanpeak_bytrial{i,j} = mean(signal{i,j}(1, lrange:urange) );
        else
            meanpeak_bytrial{i,j} = peak{i,j};
            sprintf('ROI %d trial %d has bad peak', i, j)
            meanpeak_bytrial_errors = [meanpeak_bytrial_errors; [i, j]];
        end
    end
end

end

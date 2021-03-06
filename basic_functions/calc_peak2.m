function [ meanpeak_bytrial, peakloc_bytrial] = calc_peak2( signal, startframe, endframe )
%[ meanpeak_bytrial, peakloc_bytrial, meanpeak_bytrial_errors ] = calc_peak( signal )

%calc_peak takes df/f0 and returns the location of the peak and its value 
%   input = data = df/f0, a cell array with dims num_ROIs by num_trials
%   output: peak_bytrial is a cell array with dims num_ROIs by num_trials,
%   col 1 is peak, col 2 is location

for i = 1:size(signal,1)
    for j = 1:size(signal,2)
        [peak_bytrial(i,j), peakloc(i,j)] = max(signal{i,j}(startframe:(end-endframe)));
        lrange = peakloc(i,j) - 1 + startframe;
        urange = peakloc(i,j) + 1 + startframe;
        if size(signal{1,j}, 1) ~= 1 % for this test data, this == 1, so this if statement doesn't run. WHY?
            meanpeak_bytrial(i,j) = mean(signal{i,j}(lrange:urange),1 );
        elseif size(signal{1,j}, 1) == 1
            meanpeak_bytrial(i,j) = mean(signal{i,j}(1,lrange:urange));
        end
    end
end
peakloc_bytrial_raw = peakloc + startframe
for i = 1:size(signal,1)
    for j = 1:size(signal,2)
        if length(signal{i,j}) == 160
            peakloc_bytrial(i,j) = peakloc_bytrial_raw(i,j);
        else
            len = length(signal{i,j})
            peakloc_bytrial(i,j) = peakloc_bytrial_raw(i, j) * (160/len);
        end
    end
end

end


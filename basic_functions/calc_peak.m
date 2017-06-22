function [ meanpeak_bytrial, peakloc_bytrial, peak_bytrial ] = calc_peak( signal, startframe)
%[ meanpeak_bytrial, peakloc_bytrial, meanpeak_bytrial_errors ] = calc_peak( signal )

%calc_peak takes df/f0 and returns the location of the peak and its value 
%   input = data = df/f0, a cell array with dims num_ROIs by num_trials
%   output: peak_bytrial is a cell array with dims num_ROIs by num_trials,
%   col 1 is peak, col 2 is location

for i = 1:size(signal,1)
    i
    for j = 1:size(signal,2)
        j
        [peak_bytrial(i,j), peakloc(i,j)] = max(signal{i,j}(1, startframe:(end-2)));
        lrange = peakloc(i,j) - 1 + startframe
        urange = peakloc(i,j) + 1 + startframe
        meanpeak_bytrial(i,j) = mean(signal{i,j}(1,lrange:urange) );
    end
end
peakloc_bytrial = peakloc + startframe


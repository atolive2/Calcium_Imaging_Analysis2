%% testing to figure out peakloc = 1 problem

t = 10
figure;
hold on
for k = 1:size(tadpole{1,t}.peakloc_bytrial_sm, 1)
    for j = 1:size(tadpole{1,t}.peakloc_bytrial_sm, 2)
        if tadpole{1,t}.peakloc_bytrial_sm{k,j} == 1
            plot(tadpole{1,t}.smoothed{k, j}, 'r')
        elseif tadpole{1,t}.peakloc_bytrial_sm{k,j} > 1
            plot(tadpole{1,t}.smoothed{k, j}, 'k')
        end
    end
end
hold off

%% Smooth df/f0

for t = 1:length(tadpole)
    for i = 1:size(tadpole{1,t}.df_f0, 1)
        for j = 1:size(tadpole{1,t}.df_f0, 2)
            tadpole{1,t}.smoothed_sgolay{i,j} = smooth(tadpole{1,t}.df_f0{i,j}(1,2:(end-1)), 5, 'sgolay');
        end
    end
end

t = 12
figure;
hold on
for k = 1:size(tadpole{1,t}.peakloc_bytrial_sm, 1)
    for j = 1:size(tadpole{1,t}.peakloc_bytrial_sm, 2)
        if tadpole{1,t}.boolean_response(k, j) == 1
            plot(tadpole{1,t}.smoothed_sgolay{k, j}, 'k')
        end
    end
end
hold off

figure;
hold on
for k = 1:size(tadpole{1,t}.peakloc_bytrial_sm, 1)
    for j = 1:size(tadpole{1,t}.peakloc_bytrial_sm, 2)
        if tadpole{1,t}.boolean_response(k, j) == 1
            plot(tadpole{1,t}.df_f0{k, j}, 'k')
        end
    end
end
hold off


%% Find a good ROI and work with it to make sure that good data is coming from the raw mean fluorescence calcs
% Exp 5 ROI 51
% in tadpole{1,11}
% exp 8 ROI 43 in tadpole{1,12}
t = 11
roi = 51
% raw somatic signal
figure;
hold on
for i = 1:size(tadpole{1,t}.signal, 2) 
    plot(tadpole{1,t}.signal{roi,i})
end
hold off
title(sprintf('tad %d roi %d signal', t, roi))

% raw df/f0
figure;
hold on
for i = 1:size(tadpole{1,t}.df_f0, 2) 
    plot(tadpole{1,t}.df_f0{roi,i})
end
hold off
title(sprintf('tad %d roi %d df/f0', t, roi))

% smoothed using moving average
figure;
hold on
for i = 1:size(tadpole{1,t}.smoothed, 2) 
    plot(tadpole{1,t}.smoothed{roi,i})
end
hold off
title(sprintf('tad %d roi %d df/f0 smoothed moving avg', t, roi)) 

% smoothed using sgolay
figure;
hold on
for i = 1:size(tadpole{1,t}.smoothed_sgolay, 2) 
    plot(tadpole{1,t}.smoothed_sgolay{roi,i})
end
hold off
title(sprintf('tad %d roi %d df/f0 smoothed sgolay', t, roi)) 

% smoothed using loess
figure;
hold on
for i = 1:size(tadpole{1,t}.smoothed_loess, 2) 
    plot(tadpole{1,t}.smoothed_loess{roi,i})
end
hold off
title(sprintf('tad %d roi %d df/f0 smoothed loess', t, roi)) 
% this looks worse than the raw data.

% smooth with a smaller span moving average (4 vs 8)
figure;
hold on
for i = 1:size(tadpole{1,t}.smoothed_smallmoving, 2) 
    plot(tadpole{1,t}.smoothed_smallmoving{roi,i})
end
hold off
title(sprintf('tad %d roi %d df/f0 smoothed moving 4 span', t, roi)) 

% smooth with a smaller span moving average (6 vs 8)
figure;
hold on
for i = 1:size(tadpole{1,t}.smoothed_smallmoving6, 2) 
    plot(tadpole{1,t}.smoothed_smallmoving6{roi,i})
end
hold off
title(sprintf('tad %d roi %d df/f0 smoothed moving 6 span', t, roi)) 

%% Smoothing algoritms
for t = 1:length(tadpole)
    for i = 1:size(tadpole{1,t}.df_f0, 1)
        for j = 1:size(tadpole{1,t}.df_f0, 2)
            tadpole{1,t}.smoothed_smallmoving6{i,j} = smooth(tadpole{1,t}.df_f0{i,j}(:,:), 6, 'moving');
        end
    end
end

%% None of the algorithms are better than moving with span=8
% therefore eliminate first 5 and last 5 data points from peak analysis
signal = tadpole{1,t}.smoothed
startframe=5
endframe=5
for i = 1:size(signal,1)
    for j = 1:size(signal,2)
        [peak(i,j), peakloc(i,j)] = max(signal{i,j}(startframe:(end-endframe),1));
        lrange = peakloc(i,j) - 1 + startframe;
        urange = peakloc(i,j) + 1 + startframe;
        meanpeak_bytrial(i,j) = mean(signal{i,j}(lrange:urange, 1) );
    end
end
peakloc_bytrial = peakloc + startframe
figure;
hist(reshape(peakloc_bytrial, [], 1),50)

resp_peakloc = [];
for i = 1:size(peak,1)
    for j = 1:size(peak,2)
        if peak(i,j) > 0.15
            resp_peakloc = [ resp_peakloc; peakloc_bytrial(i,j) ];
        end
    end
end
figure;
hist(reshape(resp_peakloc, [], 1),50)

% summary: moving average smoothing is fine, but must avoid first 5 points
% and last 5 points when calculating values. If all peak locs are plotted
% as a histogram, large number are ~1. But when peaklocs are only included
% when peak > 0.1, it looks more like expected (with highest count ~60,
% very few before ~40). Therefore, all data analysis of peaks must take
% into account this peak size requirement. 
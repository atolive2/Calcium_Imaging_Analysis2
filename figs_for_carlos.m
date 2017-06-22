%% Cumulative distribution of SD of time to peak

% use goodTadpole ROIs only
% std_peakloc contains the data needed for cum dist

%x = [0:(78/319):(78-78/319)]
A = std_peakloc(:,1)
AV = std_peakloc(:,2)
AM = std_peakloc(:,3)
AN = std_peakloc(:,4)
%y = cdf('normal', x, A')

% mean
BMS = avg_peakloc(:,1)
BV = avg_peakloc(:,2)
BM = avg_peakloc(:,3)
BN = avg_peakloc(:,4)


%% Is bigger peak response correlated with faster time to peak?
for i = 1:length(goodtadpole)
    meanpeak(i,:) = goodtadpole(i).peak_avg(1:4,:)
    meanpeakloc(i,:) = goodtadpole(i).peakloc_avg(1:4,:)
end

for i = 1:4
    figure;
    plot(meanpeakloc(:,i), meanpeak(:,i), 'o')
    title(sprintf('MeanPeakLoc vs MeanPeak %d', i))
end

% NO


%% Is bigger peak response correlated with more consistent time to peak?

for i = 1:4
    figure;
    plot(std_peakloc(:,i), meanpeak(:,i), 'o')
    title(sprintf('SD of PeakLoc vs MeanPeak %d', i))
end

% NO


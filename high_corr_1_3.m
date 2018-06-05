%% Re run high corr analysis with 1/3 cells threshold (to please Carlos)

% Basically, take info in allData for each tad and start over from after
% calculating xocrr. 
% Current allData is now allData1
% Run new analysis on allData, which only contains the fields up to xcorr. 

for t = 1:length(allData1)
    allData{1,t}.expnum = allData1{1,t}.expnum;
    allData{1,t}.stage = allData1{1,t}.stage;
    allData{1,t}.

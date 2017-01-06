%% plotting of calcium imaging data. Should use publish command. 


% plot area multisensory vs area of largest unisensory
Y = [-4:1:8]
X = [-4:1:8]
figure;
hold on
plot(tadpole.area_avg(1,:),tadpole.area_avg(2,:), 'b*')
plot(tadpole.area_avg(5,:),tadpole.area_avg(6,:), 'g*')
plot(tadpole.area_avg(8,:),tadpole.area_avg(9,:), 'm*')
plot(Y, X)
legend('high/high', 'high/low', 'low/low', 'Y=X')
xlabel('largest unisensory')
ylabel('multisensory')
title('PIE based on avg area')
hold off





figure;
hold on
for i = 1:size(tadpole.df_f0,1)
    if ~any(i == [tadpole.badtrials(:,:)])
        for j = 1:size(tadpole.df_f0,2)
             plot(tadpole.df_f0{i,j})
        end
    end
end
hold off

figure;
hold on
for i = 1:size(tadpole.df_f0,1)
    if tadpole.boolean_response(i,j)
        for j = 1:size(tadpole.df_f0,2)
             plot(tadpole.df_f0{i,j})
        end
    end
end
hold off

list = [1:24, 29:72]
hold on
for i = 1:length(list)
    k = list(i)
    for j = 1:size(tadpole.df_f0,2)
         plot(tadpole.df_f0{j,k})
    end
end
hold off

%% plot all no response trials with confidence intervals. WTF do these #'s
% represent???
[MUHAT,SIGMAHAT,MUCI,SIGMACI] = normfit(cell2mat(tadpole.stim_vals_area{1,4}))
hold on
hax=axes;
hist(cell2mat(tadpole.stim_vals_area{1,4}),100)
line([SIGMACI(1,1) SIGMACI(1,1)], get(hax,'YLim'),'Color',[1 0 0])
line([SIGMACI(2,1) SIGMACI(2,1)], get(hax,'YLim'),'Color',[0 1 0])
title('All area measurements for no stimulus trials')
xlabel('area')
ylabel('counts')
hold off

%% histogram of peak 
figure;
hist(cell2mat(tadpole.peak_bytrial(:,:)), 20)

% histogram of meanpeak (in theory more consistent)
hist(cell2mat(tadpole.meanpeak_bytrial(:,:)), 20)

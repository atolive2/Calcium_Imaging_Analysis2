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
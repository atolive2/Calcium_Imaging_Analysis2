%% PIE for multiple stim intensities

% Right now - exp 20 only. 
C = unique(tadpole.stimorder);
% 10 - 7; 11 - 8; 12 - 9;
figure;

% subplot 1: area - lgst uni vs MSenh all
% get data to plot, by stimulus type
area_highMhighV = [max(tadpole.area_avg(2:3,:),[],1); tadpole.area_avg(1,:)]; 
area_highMlowV = [max(tadpole.area_avg([3, 6],:),[],1); tadpole.area_avg(5,:)];
area_medMhighV = [max(tadpole.area_avg([2, 9],:),[],1); tadpole.area_avg(7,:)];
area_medMlowV = [max(tadpole.area_avg([6, 9],:),[],1); tadpole.area_avg(8,:)];
Y = -5:1:8
% get values for plotting Y=X line

subplot(2,2,1)
hold on
plot(area_highMhighV(1,:), area_highMhighV(2,:), 'r*')
plot(area_highMlowV(1,:), area_highMlowV(2,:), 'mo')
plot(area_medMhighV(1,:), area_medMhighV(2,:), 'b*')
plot(area_medMlowV(1,:), area_medMlowV(2,:), 'co')

plot(Y, Y, 'k-')
title('Area: largest uni vs MS enh')
xlabel('largest unisensory mean')
ylabel('MS enhancement mean')
legend('highMhighV', 'highMlowV', 'medMhighV', 'medMlowV', 'Y=X', 'location', 'best')
hold off

% subplot 2: area - linear sum vs MS response
% get data to plot 
area_highMhighVls = [(tadpole.area_avg(2,:)+tadpole.area_avg(3,:)); tadpole.area_avg(1,:)]; 
area_highMlowVls = [(tadpole.area_avg(3,:)+tadpole.area_avg(6,:)); tadpole.area_avg(5,:)];
area_medMhighVls = [(tadpole.area_avg(2,:)+tadpole.area_avg(9,:)); tadpole.area_avg(7,:)];
area_medMlowVls = [(tadpole.area_avg(6,:)+tadpole.area_avg(9,:)); tadpole.area_avg(8,:)];


subplot(2,2,2)
hold on
plot(area_highMhighVls(1,:), area_highMhighVls(2,:), 'r*')
plot(area_highMlowVls(1,:), area_highMlowVls(2,:), 'mo')
plot(area_medMhighVls(1,:), area_medMhighVls(2,:), 'b*')
plot(area_medMlowVls(1,:), area_medMlowVls(2,:), 'co')
plot(Y, Y, 'k-')
title('Area: linear sum vs multisensory')
xlabel('linear sum mean')
ylabel('multisensory mean')
legend('highMhighV', 'highMlowV', 'medMhighV', 'medMlowV', 'Y=X', 'location', 'best')
hold off

%%%%%%%%%%%%%%%%% by peak
% subplot 3: peak - lgst uni vs MSenh all
% get data to plot, by stimulus type
peak_highMhighV = [max(tadpole.peak_avg(2:3,:),[],1); tadpole.peak_avg(1,:)]; 
peak_highMlowV = [max(tadpole.peak_avg([3, 6],:),[],1); tadpole.peak_avg(5,:)];
peak_medMhighV = [max(tadpole.peak_avg([2, 9],:),[],1); tadpole.peak_avg(7,:)];
peak_medMlowV = [max(tadpole.peak_avg([6, 9],:),[],1); tadpole.peak_avg(8,:)];
Y_2 = -0.05:0.01:0.45;
% get values for plotting Y=X line

subplot(2,2,3)
hold on
plot(peak_highMhighV(1,:), peak_highMhighV(2,:), 'r*')
plot(peak_highMlowV(1,:), peak_highMlowV(2,:), 'mo')
plot(peak_medMhighV(1,:), peak_medMhighV(2,:), 'b*')
plot(peak_medMlowV(1,:), peak_medMlowV(2,:), 'co')

plot(Y_2, Y_2, 'k-')
title('Peak: largest uni vs MS enh')
xlabel('largest unisensory mean')
ylabel('MS enhancement mean')
legend('highMhighV', 'highMlowV', 'medMhighV', 'medMlowV', 'Y=X', 'location', 'best')
hold off

% subplot 4: peak - linear sum vs MS response
% get data to plot 
peak_highMhighVls = [(tadpole.peak_avg(2,:)+tadpole.peak_avg(3,:)); tadpole.peak_avg(1,:)]; 
peak_highMlowVls = [(tadpole.peak_avg(3,:)+tadpole.peak_avg(6,:)); tadpole.peak_avg(5,:)];
peak_medMhighVls = [(tadpole.peak_avg(2,:)+tadpole.peak_avg(9,:)); tadpole.peak_avg(7,:)];
peak_medMlowVls = [(tadpole.peak_avg(6,:)+tadpole.peak_avg(9,:)); tadpole.peak_avg(8,:)];


subplot(2,2,4)
hold on
plot(peak_highMhighVls(1,:), peak_highMhighVls(2,:), 'r*')
plot(peak_highMlowVls(1,:), peak_highMlowVls(2,:), 'mo')
plot(peak_medMhighVls(1,:), peak_medMhighVls(2,:), 'b*')
plot(peak_medMlowVls(1,:), peak_medMlowVls(2,:), 'co')
plot(Y_2, Y_2, 'k-')
title('Peak: linear sum vs multisensory')
xlabel('linear sum mean')
ylabel('multisensory mean')
legend('highMhighV', 'highMlowV', 'medMhighV', 'medMlowV', 'Y=X', 'location', 'best')
hold off


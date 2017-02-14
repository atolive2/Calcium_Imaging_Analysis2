
% TUTORIAL: Dimensionality reduction using PCA and classification using KNN
%
% Carlos Vargas-Irwin, Brown University, Jan 18 2017


% DATA

%Let's make some fake data representing 100 trials with 160 data points each 
data = rand(100,160);
% lets make some of the trials a bit different
% (representing 3 separate classes)
data(1:25,1:50) = data(1:25,1:50) + 0.2;
data(26:50,50:100) = data(26:50,50:100) - 0.2;
% and keep track of which are which
trial_type = [zeros(1,25)+1 zeros(1,25)+2 zeros(1,50)+3];

% Here's what the raw data looks like
figure
subplot(311)
plot(data(1:25,:)','r')
axis([1 160 -0.5 1.5 ])
title('Raw Data')
subplot(312)
plot(data(26:50,:)','b')
axis([1 160 -0.5 1.5 ])
subplot(313)
plot(data(51:100,:)','k')
axis([1 160 -0.5 1.5 ])


% DIMENSIONALITY REDUCTION

% use principal component analysis
[coeff, score] = pca(data);
% Isolate the projections from the first 3 PCs
data3D = score(:,1:3);


% Plot the 3D data, colored by trial type
figure
f1 = find(trial_type == 1)
h = plot3(data3D(f1,1),data3D(f1,2),data3D(f1,3),'ro'); set(h,'MarkerFaceColor','r')
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
hold on
f2 = find(trial_type == 2)
h = plot3(data3D(f2,1),data3D(f2,2),data3D(f2,3),'bo'); set(h,'MarkerFaceColor','b')
f3 = find(trial_type == 3)
h = plot3(data3D(f3,1),data3D(f3,2),data3D(f3,3),'ko'); set(h,'MarkerFaceColor','k')
legend([{'class 1'} {'class 2'} {'class 3'}])


% CLASSIFICATION

% Build a KNN classifier and test using 10-fold cross validation
% (We are giving it 3 different class labels, 
% so it will be a 3-way classification)
KNNClass =fitcknn(data3D,trial_type) 
CVKNN = crossval(KNNClass,'kfold',10);

% and report the mean and sdev for classification accuracy
klossKNN = kfoldLoss(CVKNN,'mode','individual');
accuracyM = mean(1-klossKNN) *100
accuracySTD = std(1-klossKNN) *100

title(['Decoding Accuracy: ' num2str(accuracyM) '% (+/-' num2str(accuracySTD)  ')'])

% What should we expect by chance? 
% Lets' shuffle the trial labels using 'randperm'
CHKNNClass =fitcknn(data3D,trial_type(randperm(numel(trial_type))));
CHRCVKNN = crossval(CHKNNClass,'kfold',10);
CHklossKNN = kfoldLoss(CHRCVKNN,'mode','individual');

CHaccuracyM = mean(1-CHklossKNN) *100
CHaccuracySTD = std(1-CHklossKNN) *100
% NOTE - for real data we would repeat the label shuffling step about 1000
% times or so, to get a more accurate estimate of the chance distribution


% Bonus Question: 
% what happens if we run the classifier using the original 160 dimensions?
% (using 'data' instead of 'data3D')


%% How much information do individual channels have?
klossKNNIC = zeros(size(data,2),1);
accuracyMIC = klossKNNIC;
accuracySTDIC = klossKNNIC;
for n = 1:size(data,2)
KNNClass =fitcknn(data(:,n),trial_type);
CVKNN = crossval(KNNClass,'kfold',10);
% and report the mean and sdev for classification accuracy
klossKNNIC = kfoldLoss(CVKNN,'mode','individual');
accuracyMIC(n) = mean(1-klossKNNIC) *100;
accuracySTDIC(n) = std(1-klossKNNIC) *100;
end

figure
hist(accuracyMIC,15)
xlabel('classification accuracy')
ylabel('# of channels')
title(['Single channel classification accuracy (ensemble =' num2str(accuracyM) '%)' ] )



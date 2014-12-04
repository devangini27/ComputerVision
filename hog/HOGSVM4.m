clear all
close all

imageextension = '.bmp';
imageextension2 = '.ppm';
% pedestrianspath = 'images//pos';
pedestrianspath = 'pedestrians128x64';
pedestriansimageprefix = [pedestrianspath '//per'];
personsimageprefix = 'persons//person_';
% bikespath = 'bikes';
bikespath = 'bikes';
bikesimageprefix = [ bikespath '//bike_'];

count3 = 15;
count4 = 15;
correctclasses= [ones(1 , count3 ) zeros(1, count4)]';

%check if there are any samples
listing = dir(pedestrianspath);
[samplesSize1, ~] = size(listing);
samplesSize1 = samplesSize1 - 2;
count1 = samplesSize1 - count3;

listing = dir(bikespath);
[samplesSize2, ~] = size(listing);
samplesSize2 = samplesSize2 - 2;
count2 = samplesSize2 - count4;


hogimagesfeatures = zeros(count1 + count2, 3780 ); %extra one for storing the class
classnumber = zeros(count1+count2,1);

avgmagnitude = zeros(128, 64);


%get all the hog features
formatSpec1 = '%05d';
for i = 1 : count1
    trainingimage = sprintf(formatSpec1,i);
    imagename = [pedestriansimageprefix trainingimage imageextension2];
    disp(imagename);
    [hogfeatures, logicalcellsindex, image, magnitude] = HOGFeature(imagename, 0, 1);
     avgmagnitude = avgmagnitude + magnitude;
    %     disp(hog);
    hogimagesfeatures(i,:) = hogfeatures;
    classnumber(i,1) = 1;
    
    %HOGVisualize(image, hogfeatures, logicalcellsindex );
    
end
avgmagnitude = avgmagnitude / count1;
figure;
imshow(uint8(avgmagnitude));

formatSpec2 = '%03d';
for i = 1 : count2
    trainingimage = sprintf(formatSpec2,i);
    imagename = [bikesimageprefix trainingimage imageextension];
    disp(imagename);
    [hogfeatures, logicalcellsindex, image, ~] = HOGFeature(imagename, 1, 0 );
    %     disp(hog);
    hogimagesfeatures(count1 + i,:) = hogfeatures;
    classnumber(count1 + i,1) = 0;
end

testhogimagesfeatures = zeros(count3 + count4, 3780 ); %extra one for storing the class
testclassnumber = zeros(count3 + count4,1);

for i = 1 : count3 + count4
    ppmmode = 0;
    if correctclasses(i,1) == 1
        %         testingimage = sprintf(formatSpec1,i);
        testingimage = sprintf(formatSpec1,count1+i);
        imagename = [pedestriansimageprefix testingimage imageextension2];
        ppmmode = 1;
    else
        %          testingimage = sprintf(formatSpec2, i-count3);
        testingimage = sprintf(formatSpec2,count2+i-count3);
        imagename = [bikesimageprefix testingimage imageextension];
    end
    disp(imagename);
    [hogfeatures, logicalcellsindex, image, ~] = HOGFeature(imagename, 1, ppmmode);
    %         disp(hogfeatures(1,1:5));
    testhogimagesfeatures(i,:) = hogfeatures;
end


%Train, and optionally cross validate, an SVM classifier using fitcsvm. The most common syntax is:
data =  hogimagesfeatures;
groups = classnumber;

%# number of cross-validation folds:
%# If you have 50 samples, divide them into 10 groups of 5 samples each,
%# then train with 9 groups (45 samples) and test with 1 group (5 samples).
%# This is repeated ten times, with each group used exactly once as a test set.
%# Finally the 10 results from the folds are averaged to produce a single 
%# performance estimation.
k=10;

%Randomly select training and test sets.
% cvFolds = crossvalind('Kfold', groups, k);   %# get indices of 10-fold CV
[train, test] = crossvalind('holdOut',groups);
cp = classperf(groups); %# init performance tracker

%5. Train an SVM classifier using a linear kernel function and plot the grouped data.
% svmStruct = svmtrain(data(train,:),groups(train),'kernel_function','linear');
svmStruct = svmtrain(data(train,:),groups(train),'kernel_function','linear','boxconstraint',0.01);
% svmStruct = svmtrain(data(train,:),groups(train),'kernel_function','rbf','boxconstraint',0.01);

%  Use the svmclassify function to classify the test set.
classes = svmclassify(svmStruct,data(test,:));


% for i = 1:k                                  %# for each fold
%     testIdx = (cvFolds == i);                %# get indices of test instances
%     trainIdx = ~testIdx;                     %# get indices training instances
% 
%     %# train an SVM model over training instances
%     svmStruct = svmtrain(data(trainIdx,:), groups(trainIdx), 'kernel_function','linear','boxconstraint',0.01 );
% 
%     %# test using test instances
%     pred = svmclassify(svmStruct, data(testIdx,:));
% 
%     %# evaluate and update performance object
%     cp = classperf(cp, pred, testIdx);
% end



% Evaluate the performance of the classifier.
classperf(cp,classes,test);
disp(['performance = ' num2str(cp.CorrectRate)])

cp.CountingMatrix

performance = 0;

for i = 1: count3 + count4
    
    testclassnumber(i,1) = svmclassify(svmStruct,testhogimagesfeatures(i,:));
    disp(['obtained classs = ' num2str(testclassnumber(i,1)) ' and correct class = ' num2str(correctclasses(i,1))]);
    
    if testclassnumber(i,1) == correctclasses(i,1)
        performance = performance + 1;
    end
    
end
performance = performance / (count3 +  count4);

disp(['testing performance on unseen data = ' num2str(performance) ]);


weights = 1000 * svmStruct.Alpha;
[count, ~] = size(weights);
% svmStruct.SupportVectors

positiveweights = zeros(1,3780);
negativeweights = zeros(1,3780);

for index = 1 : count
    weight = weights(index, 1);
    if weight > 0
        positiveweights = positiveweights + weight * svmStruct.SupportVectors(index, :);
    else
        negativeweights = negativeweights - weight * svmStruct.SupportVectors(index, :);
    end
end

rows = 15;
cols = 7;
positiveimage = zeros(rows,cols);
negativeimage = zeros(rows,cols);
blocksize = (2 * 2) * 9;

for i  = 1 : rows
    for j = 1 : cols
        index = (i - 1) * cols + j;
        
        blockindex = (index - 1) * blocksize + 1;
%         disp(['i = ' num2str(i) ' , j = ' num2str(j) ' , index = ' num2str(index) ' , blockindex = ' num2str(blockindex)]);
        %split that part of the positiveweights hogfeatures
%         disp(['extracting from index ' num2str(blockindex) ' to ' num2str(blockindex + blocksize - 1)]);
        blockhistograms = positiveweights(blockindex : blockindex + blocksize - 1 );
        maxvalue = max(blockhistograms);
        positiveimage(i,j) = maxvalue;  
        blockhistograms = negativeweights(blockindex : blockindex + blocksize - 1 );
        maxvalue = max(blockhistograms);
        negativeimage(i,j) = maxvalue;  
    end
end


figure;
imshow(uint8(255*positiveimage));

figure;
imshow(uint8(255*negativeimage));





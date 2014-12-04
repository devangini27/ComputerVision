clear all
close all


imageextension = '.bmp';
imageextension2 = '.ppm';
pedestriansimageprefix = 'pedestrians128x64//per';
personsimageprefix = 'persons//person_';
bikesimageprefix = 'bikes//bike_';

trainingimages1 = ['00001'
    '00002'
    '00003'
    '00004'
    '00005'
    '00006'
    '00007'
    '00008'
    '00009'
    '00010'];
% trainingimages1 = ['029'
%     '034'
%     '028'
%     '032'
%     '036'
%     '037'
%     '038'
%     '091'
%     '092'
%     '096'];
trainingimages2 = ['011'
    '012'
    '013'
    '015'
    '016'
    '017'
    '020'
    '021'
    '022'
    '023'
    ];

testingimages1 = [ '00011'
    '00012'
    '00013'
    '00014'
    '00015'];
testingimages2 = [
    '083'
    '085'
    '086'
    '087'
    '088'
    
    ];
correctclasses= [1,1,1,1,1,0,0,0,0,0]';

[count1, ~] = size(trainingimages1);
[count2, ~] = size(trainingimages2);
[count3, ~] = size(testingimages1);
[count4, ~] = size(testingimages2);

hogimagesfeatures = zeros(count1 + count2, 3780 ); %extra one for storing the class
classnumber = zeros(count1+count2,1);
avgmagnitude = zeros(128, 64);


% hogimagesfeatures = zeros(count1 , 3780 );
% classnumber = zeros(count1 ,1);

%get all the hog features
for i = 1 : count1
    imagename = [pedestriansimageprefix trainingimages1(i,:) imageextension2];
    %      imagename = [personsimageprefix trainingimages1(i,:) imageextension];
    disp(imagename);
    [hogfeatures, logicalcellsindex, image, magnitude] = HOGFeature(imagename, 0, 1);
    avgmagnitude = avgmagnitude + magnitude;
    %     disp(hog);
    hogimagesfeatures(i,1:end) = hogfeatures;
    classnumber(i,1) = 1;
    
    HOGVisualize(image, hogfeatures, logicalcellsindex );
    
end
avgmagnitude = avgmagnitude / count1;

figure;
imshow(uint8(avgmagnitude));

for i = 1 : count2
    imagename = [bikesimageprefix trainingimages2(i,:) imageextension];
    disp(imagename);
    [hogfeatures, logicalcellsindex, image, ~] = HOGFeature(imagename, 1, 0 );
    %     disp(hog);
    hogimagesfeatures(count1 + i,1:end) = hogfeatures;
    classnumber(count1 + i,1) = 0;
end

testhogimagesfeatures = zeros(count2 + count3, 3780 ); %extra one for storing the class
testclassnumber = zeros(count3 + count4,1);

for i = 1 : count3 + count4
    ppmmode = 0;
    if correctclasses(i,1) == 1
        imagename = [pedestriansimageprefix testingimages1(i,:) imageextension2];
        ppmmode = 1;
    else
        imagename = [bikesimageprefix testingimages2(i - count3,:) imageextension];
    end
    disp(imagename);
    [hogfeatures, logicalcellsindex, image, ~] = HOGFeature(imagename, 1, ppmmode);
    %     disp(hog);
    testhogimagesfeatures(i,1:end) = hogfeatures;
end


% train the svm
% http://www.egr.msu.edu/classes/ece480/capstone/spring11/group04/applicati
% on_Kan.pdf

%Create data, a two-column matrix containing sepal length and sepal width
%measurements for 150 irises.
data =  hogimagesfeatures;
% From the species vector, create a new column vector, groups, to classify data
%into two groups: data and non-data.
groups = classnumber;
%4. Randomly select training and test sets.
[train, test] = crossvalind('holdOut',groups);
cp = classperf(groups);
%5. Train an SVM classifier using a linear kernel function and plot the grouped data.
% svmStruct = svmtrain(data(train,:),groups(train),'kernel_function','linear');
svmStruct = svmtrain(data(train,:),groups(train),'kernel_function','linear','boxconstraint',0.01);
% svmStruct = svmtrain(data(train,:),groups(train),'kernel_function','rbf','boxconstraint',0.01);
%
% title(sprintf('Kernel Function: %s',...
%     func2str(svmStruct.KernelFunction)),...
%     'interpreter','none');
disp(['Kernel Function: ' func2str(svmStruct.KernelFunction)]);
%7. Use the svmclassify function to classify the test set.
classes = svmclassify(svmStruct,data(test,:));


%8. Evaluate the performance of the classifier.
classperf(cp,classes,test);
disp(['performance = ' num2str(cp.CorrectRate)])

%9. Use a one-norm, hard margin support vector machine classifier by changing the
%boxconstraint property.
% figure
% svmStruct = svmtrain(data(train,:),groups(train),...
%     'boxconstraint',1e6);
% classes = svmclassify(svmStruct,data(test,:));

% % try the proper test data
% classes1 = svmclassify(svmStruct,data(test,:),'showplot',true);
%
% %10. Evaluate the performance of the classifier.
% classperf(cp,classes1,correctclasses);
% cp.CorrectRate


for i = 1: count3 + count4
    
    classnumber = svmclassify(svmStruct,testhogimagesfeatures(i,:));
    disp(['obtained classs = ' num2str(classnumber) ' and correct class = ' num2str(correctclasses(i,1))]);
    
end

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

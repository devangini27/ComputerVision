clear all
close all
 addpath('svm_mex601//matlab');

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

hogimagesfeatures = zeros(count1 + count2, 3780 + 1 ); %extra one for storing the class
classnumber = zeros(count1+count2,1);

%get all the hog features
for i = 1 : count1
    imagename = [pedestriansimageprefix trainingimages1(i,:) imageextension2];
    disp(imagename);
    [hogfeatures, logicalcellsindex, image, magnitude] = HOGFeature(imagename, 0, 1);
    %     disp(hog);
    hogimagesfeatures(i,2:end ) = hogfeatures;
    hogimagesfeatures(i,1) = 1;
    
    %HOGVisualize(image, hogfeatures, logicalcellsindex );
    
end
for i = 1 : count2
    imagename = [bikesimageprefix trainingimages2(i,:) imageextension];
    disp(imagename);
    [hogfeatures, logicalcellsindex, image, ~] = HOGFeature(imagename, 1, 0 );
    %     disp(hog);
    hogimagesfeatures(count1 + i,2:end) = hogfeatures;
    hogimagesfeatures(count1 + i,1) = -1;
end

testhogimagesfeatures = zeros(count3 + count4, 3780 ); %extra one for storing the class
testclassnumber = zeros(count3 + count4,1);

for i = 1 : count3 + count4
    ppmmode = 0;
    if correctclasses(i,1) == 1
        imagename = [pedestriansimageprefix testingimages1(i,:) imageextension2];
        ppmmode = 1;
    else
        imagename = [bikesimageprefix testingimages2(i - count3 ,:) imageextension];
    end
    disp(imagename);
    [hogfeatures, logicalcellsindex, image, ~] = HOGFeature(imagename, 1, ppmmode);
%         disp(hogfeatures(1,1:5));
    testhogimagesfeatures(i,1:end) = hogfeatures;
end


%Train, and optionally cross validate, an SVM classifier using fitcsvm. The most common syntax is:
x = hogimagesfeatures;
y = classnumber;

%create a file

fileID = fopen('training_file','w');
formatSpec1 = '%d ';
formatSpec2 = '%d:%d ';
formatSpec3 = '\n';

[objectcount, featurecount] = size(hogimagesfeatures);
for i = 1 : objectcount
fprintf(fileID,formatSpec1,hogimagesfeatures(i,1));
for j = 2 : featurecount
fprintf(fileID,formatSpec2,[(j-1) ; hogimagesfeatures(i,j)]);
end
fprintf(fileID,formatSpec3);
end
fclose(fileID);

%create test file

fileID = fopen('test_file','w');

[objectcount, featurecount] = size(testhogimagesfeatures);
for i = 1 : objectcount
fprintf(fileID,formatSpec1,0);
for j = 1 : featurecount
fprintf(fileID,formatSpec2,[(j-1) ; testhogimagesfeatures(i,j)]);
end
fprintf(fileID,formatSpec3);
end
fclose(fileID);



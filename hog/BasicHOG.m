clear all
close all

imageextension2 = '.ppm';
pedestriansimageprefix = 'pedestrians128x64//per';

trainingimages1 = '00001';

%get all the hog features
imagename = [pedestriansimageprefix trainingimages1 imageextension2];
disp(imagename);
[hogfeatures, logicalcellsindex, image, magnitude] = HOGFeature(imagename, 0, 1);

HOGVisualize(image, hogfeatures, logicalcellsindex );



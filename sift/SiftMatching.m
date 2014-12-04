close all
clear all

imagename = 'daisy-photo.png';

image = imread(imagename);
[~,cols,channels] = size(image);
if channels == 3
    image = rgb2gray(image);
end

%split into two parts
image1 = image(:, 1 : cols/2);
image2 = image(:, cols/2 + 1: cols);

% figure(1);
% imshow(image1);
% 
% figure(2);
% imshow(image2);

% features1 = SIFTDescriptor('Yosemite1.jpg');
% features2 = SIFTDescriptor('Yosemite2.jpg');

features1 = SIFTDescriptor(image1);
% features2 = SIFTDescriptor(image2);
% 
% %find closest matching match
% 
% [features1size,~] = imsize(features1);
% [features2size,~] = imsize(features2);
% 
% 
% thresholdDistance = 5;
% matchingfeatures= [];
% 
% for i = 1:features1size
%     
%     minindex = -1;
%     minvalue = 100000;
%     
%     minindex2 = -1;
%     minvalue2 = 100000;
%     
%     for j = 1:features2size
%         %find distance between the two features
%         
%         distance = 0;
%         for k = 1:128
%             diff = features1(i,k) - features2(j,k);
%             distance = distance + (diff)^2;
%         end
%         distance = sqrt(distance);
%         
%         if distance < minvalue
%             minindex = j;
%             minvalue = distance;
%         end
%         
%     end
%     
%     
%     if minvalue <= thresholdDistance
%         matchingfeatures = [matchingfeatures ; [i minindex]];
%     end
% end
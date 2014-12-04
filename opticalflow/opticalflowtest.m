close all
clear all

width = 10;
image1 = zeros(width, width);
startx = 3;
starty = 3;
image1(startx,starty) = 255;
image1(startx,starty+1) = 50;
image1(startx+1,starty) = 60;
image1(startx+1,starty+1) = 100;
% image1(5,5) =255;


motionx = 0;
motiony = 2;
image2 = zeros(width, width);
% image2(6,4) = 255;
% image2(6,5) = 255;
% image2(7,4) = 255;
% image2(7,5) = 255;
image2(startx+motionx,starty+motiony) = 255;
image2(startx+motionx,starty+motiony+1) = 50;
image2(startx+motionx+1,starty+motiony) = 60;
image2(startx+motionx+1,starty+motiony+1) = 100;


% image2(5,6) =255;
%
% figure;
% imshow(image1);
%
% figure;
% imshow(image2);

% imagename1 = 'sphere1.jpg';
% imagename2 = 'sphere2.jpg';

% imagename1 = 'anim.00.tif';
% imagename2 = 'anim.01.tif';

% imagename1 = 'fg001.pgm';
% imagename2 = 'fg003.pgm';
% 
% image1 = imread(imagename1, 'pgm');
% image2 = imread(imagename2, 'pgm');
% 
% [~,~,channels] = size(image1);
% if channels == 3
%     image1 = rgb2gray(image1);
% end
% [~,~,channels] = size(image2);
% if channels == 3
%     image2 = rgb2gray(image2);
% end

figure;
imshow(image1);
figure;
imshow(image2);

% [u,v] = hornschunck1(image1, image2);

% [u,v] = lucaskanade(image1, image2);

[u,v] = lucaskanadepyramid(image1, image2);

[rows, cols] = size(image1);

[x,y] = meshgrid(1:cols, 1:rows);

figure;
%     imshow(image1);

%     [x2, y2] = meshgrid(1 : 5 : rows, 1 : 5 : cols);
%     u2 = u(x2,y2);
%     v2 = v(x2,y2);


hold on;
interval = 1;
x2 = x(1:interval:end,1:interval:end);
y2 = y(1:interval:end,1:interval:end);
u2 = u(1:interval:end,1:interval:end);
v2 = v(1:interval:end,1:interval:end);


% quiver(x2,y2,u2,v2);
    quiver(u,v);



% lucaskanadeiterative(image1, image2);

% 
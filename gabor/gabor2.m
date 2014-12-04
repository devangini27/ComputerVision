% image = imread('building2.jpg');
image = imread('checkboard.png');
% image = imread('eyeball.jpg');
image = rgb2gray(image);


%function [gb] = gabor1(sigma, gamma, theta, lambda, psi)
gabor = gabor1(6, 1, 0, 5, 0);

% newimage = conv2(double(image), gabor, 'same');
% newimage = newimage / norm(newimage,1);

newimage = imfilter(image, gabor, 'symmetric');   

figure(1);
imshow(uint8(newimage));
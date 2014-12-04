image = imread('checkboard.png');
image = rgb2gray(image);
figure(1);
imshow(image);
[rows, cols] = size(image);

%define some constants
k = 0.2;
threshold = 20;

% step 1 - computer derivative of image

filterdx = [-1 0 1; -1 0 1; -1 0 1];   % Derivative masks
filterdy = filterdx';

dx = conv2(double(image), filterdx, 'same');
dy = conv2(double(image), filterdy, 'same');
figure(2);
imshow(uint8(dx));
title('dx');
figure(3);
imshow(uint8(dy));
title('dy');

% step 2 - compute IxIx, IxIy, IyIy
Ix2 = dx.*dx;
Ixy = dx.*dy;
Iy2 = dy.*dy;

% step 3 - compute gaussian
% Generate Gaussian filter of size 6*sigma (+/- 3sigma) and of
% minimum size 1x1.
sigma=1;
gaussian = fspecial('gaussian',max(1,fix(6*sigma)), sigma);
%gaussian = fspecial('gaussian', [rows cols], size/6);

Ix2 = conv2(Ix2, gaussian, 'same'); % Smoothed squared image derivatives
Iy2 = conv2(Iy2, gaussian, 'same');
Ixy = conv2(Ixy, gaussian, 'same');


% for every search window calculate scalar cornerness value
windowSize = 12;
%calculate scalar cornerness
%k = 0.04;
%harris = (Ix2.*Iy2 - Ixy.^2) - k*(Ix2 + Iy2).^2; % Original Harris measure.
harris = (Ix2.*Iy2 - Ixy.^2)./(Ix2 + Iy2 + eps); % My preferred  measure.

% for i = 1:rows - windowSize
%     for j = 1:cols - windowSize
%step 4 - calculate M
% M = w(x,y) * [IxIx IxIy; IxIy IyIy]

%         m = zeros(2,2);
%         for x = i : i + windowSize - 1
%             for y = j : j +windowSize - 1
%
%                m = m + gaussian(i,j) * [ Ix2(i,j) Ixy(i,j) ; Ixy(i,j) Iy2(i,j) ];
%
%             end
%         end
%         disp(m);

%         % step 5 - calculate eigenvalues of m
%         eigenvalues = eig(m);
%         disp(eigenvalues);
%
%         %step 6 - calculate scalar measure R
%         sum = eigenvalues(1) + eigenvalues(2) ;
%         r  = eigenvalues(1)*eigenvalues(2) - k * (sum) * (sum);
%         disp(r);
%
%         if r > threshold
%             disp('found corner');
%         end
%     end
% end


% find local minima above threshold
% Extract local maxima by performing a grey scale morphological
% dilation and then finding points in the corner strength image that
% match the dilated image and are also greater than the threshold.\
radius =5;
thresh = 5000;
size1 = 2*radius+1;                   % Size of mask.
mx = ordfilt2(harris,size1^2,ones(size1)); % Grey-scale dilate.
maxima = (harris==mx)&(harris>thresh);       % Find maxima.

[r,c] = find(maxima);                  % Find row,col coords.

% overlay corners on original image
figure(4), imagesc(image), axis image, colormap(gray), hold on
plot(c,r,'ys'), title('corners detected');

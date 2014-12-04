image1 = imread('butterfly2.jpg');
image = rgb2gray(image1);
figure(1);
imshow(image);
[rows, cols] = size(image);

sigmavalue=1.6;
kvalue=sqrt(2.0);

minlevel=1;
maxlevel=5;
length = maxlevel-minlevel;
logvalues = zeros(rows,cols,length);

for i = minlevel: maxlevel
    sigma = sigmavalue ^ i;
    width = floor(6 * sigma + 1);
    log = fspecial('log', width, sigma);
    smoothedvalues = conv2(double(image), log, 'same');
    
    figure(i+1);
    imshow(smoothedvalues);
    
    logvalues(:,:,i) = smoothedvalues(:,:);
    %logvalues = [logvalues; [smoothedvalues]];
end

centerlevel = (minlevel+maxlevel)/2;
locations = zeros(rows,cols);
count=zeros(rows,cols);

for centerlevel = minlevel + 1 : maxlevel - 1
    for i = 2: rows - 1
        for j = 2:cols -1
            %find if the centerlevel pixel is the local maxima or minima
            centerpixel  = logvalues(i,j,centerlevel);
            neighbours = logvalues(i-1:i+1, j-1:j+1, centerlevel-1:centerlevel+1);
            neighbours= neighbours(:);
            neighbours(14) = [];
            
            %find if greater than all values or smaller than all values
            maxvalue = max(neighbours);
            minvalue = min(neighbours);
            
            if (centerpixel > maxvalue)
                locations(i,j) = 1;
                count(i,j) = count(i,j)+1;
            elseif (centerpixel < minvalue)
                locations(i,j) = 1;
                count(i,j) = count(i,j)+1;
            end
        end
    end
    
end

requiredlevels = maxlevel - minlevel - 2;

[r,c] = find(locations);                  % Find row,col coords.

% overlay corners on original image
figure(maxlevel+1 + centerlevel),
imagesc(image),
axis image,
colormap(gray),
hold on;
plot(c,r,'ys'),
title('corners detected');



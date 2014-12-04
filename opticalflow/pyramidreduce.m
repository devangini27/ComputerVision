function [reducedimage2] = pyramidreduce(image)

% image = imread('Lenna.png');
% image = imread(imagename);

[rows, cols,channels] = size(image);
if channels == 3
image = rgb2gray(image);
end

figure;
imshow(image);



trianglemask = [0 0.25 0.5 0.25 0];
gaussianmask = trianglemask;%[0.05 0.25 0.4 0.25 0.05];

reducedimage1 = zeros(rows, floor(cols/2));

for i = 1 : rows
    for j = 2 : 2 : cols
        sum = 0;
        if j - 2 >= 1
            sum = sum + gaussianmask(1) * image(i,j-2);
        end
        if j - 1 >= 1
            sum = sum + gaussianmask(2) * image(i,j-1);
        end
        sum = sum + gaussianmask(3) * image(i,j);
        if j + 1 <= cols
            sum = sum + gaussianmask(4) * image(i,j+1);
        end
        if j + 2 <= cols
            sum = sum + gaussianmask(5) * image(i,j+2);
        end
        
        reducedimage1(i,j/2) = sum;
    end
end

figure;
imshow(uint8(reducedimage1));

reducedimage2 = zeros(floor(rows/2), floor(cols/2));

for j = 1 : cols/2
    for i = 2 : 2 : rows
        sum = 0;
        if i - 2 >= 1
            sum = sum + gaussianmask(1) * reducedimage1(i-2,j);
        end
        if i - 1 >= 1
            sum = sum + gaussianmask(2) * reducedimage1(i-1,j);
        end
        sum = sum + gaussianmask(3) * reducedimage1(i,j);
        if i + 1 <= rows
            sum = sum + gaussianmask(4) * reducedimage1(i+1,j);
        end
        if i + 2 <= rows
            sum = sum + gaussianmask(5) * reducedimage1(i+2,j);
        end
        
        reducedimage2(i/2,j) = sum;
    end
end


figure;
imshow(uint8(reducedimage2));


end
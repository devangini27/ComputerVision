image1 = imread('butterfly2.jpg');
image = rgb2gray(image1);
figure(1);
imshow(image);
[rows, cols] = size(image);

sigmavalue=1.6;
kvalue=sqrt(2.0);

minlevel=1;
maxlevel=8;
length = maxlevel-minlevel;
gaussianvalues = zeros(rows,cols,length+1);
dogvalues = zeros(rows, cols, length);



sigma = sigmavalue;
width = floor(6 * sigma + 1);
gfilter = fspecial('gaussian', width, sigma);
gaussianvalues(:,:,minlevel) = conv2(double(image), gfilter, 'same');


for i = minlevel: maxlevel
    sigma = sigmavalue * (kvalue ^ i);
    width = floor(6 * sigma + 1);
    gfilter = fspecial('log', width, sigma);
    gaussianvalues(:,:,i+1) = conv2(double(image), gfilter, 'same');
    
    figure(i+1);
    imshow(gaussianvalues(:,:,i+1));
    
    dogvalues(:,:,i) =  gaussianvalues(:,:,i+1) -  gaussianvalues(:,:,i);
    
end

locations = zeros(rows,cols, length);
count=zeros(rows,cols);

for centerlevel = minlevel + 1 : maxlevel - 2
    for i = 2: rows - 1
        for j = 2:cols -1
            %find if the centerlevel pixel is the local maxima or minima
            centerpixel  = dogvalues(i,j,centerlevel);
            neighbours = dogvalues(i-1:i+1, j-1:j+1, centerlevel-1:centerlevel+1);
            neighbours= neighbours(:);
            neighbours(14) = [];
            
            %find if greater than all values or smaller than all values
            maxvalue = max(neighbours);
            minvalue = min(neighbours);
            
            if (centerpixel > maxvalue)
                locations(i,j, centerlevel) = 1;
                count(i,j) = count(i,j)+1;
            elseif (centerpixel < minvalue)
                locations(i,j, centerlevel) = 1;
                count(i,j) = count(i,j)+1;
            end
        end
    end
    
end


for centerlevel = minlevel + 1 : maxlevel - 2
    [r,c] = find(locations(:,:,centerlevel));                  % Find row,col coords.
    
    % overlay corners on original image
    figure(maxlevel+1 + centerlevel),
    imagesc(image),
    axis image,
    colormap(gray),
    hold on;
    %plot(c,r,'ys'),
    title('corners detected');
    
    [length,~] = size(r);
    sigma = sigmavalue * (kvalue ^ centerlevel);
    radius = sigma;
    for i = 1 : length
        y = r(i);
        x = c(i);
        ang=0:0.01:2*pi;
        xp=radius*cos(ang);
        yp=radius*sin(ang);
        plot(x+xp,y+yp);
    end
end


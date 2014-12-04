function [histogram, logicalcellsindex, image, magnitude] = HOGFeature(imagename, rotate,ppmmode)


if ppmmode
    image = imread(imagename, 'ppm');
else
image = imread(imagename);%'persons//person_015.bmp');
end
image = rgb2gray(image);

% figure;
% imshow(image);

if rotate
    image = imrotate(image,90);
end

image = imresize(image, [128,64]);
%image = image();
[rows, cols] = size(image);


% figure;
% imshow(image);

%compute gradients without smoothing
filterdx= [-1 0 1];
filterdy = filterdx';

dx = conv2(double(image), filterdx, 'same');

% figure;
% imshow(dx);


dy = conv2(double(image), filterdy, 'same');

% figure;
% imshow(dy);

%compute gradient magnitude and orienation
magnitude = sqrt(dx.*dx+ dy.*dy);


% figure;
% imshow(magnitude);


dx(find(dx == 0)) = 1;
dy(find(dx == 0)) = 1000000;
orientation = atand(dy./dx);
%convert then angles from -90 to 90  -> 0 to 180
% -90 -- > 90, -45 -> 135,
for i = 1: rows
    for j = 1: cols
        if orientation(i,j) < 0
            orientation(i,j) = 180 + orientation(i,j);
        end
    end
end

% cellindexx = zeros(128,64);
% cellindexy = zeros(128,64);
% blockindexx = zeros(128,64);
% blockindexy = zeros(128,64);

% cellindexx = [];% zeros(3780);
% cellindexy = [];% zeros(3780);
% blockindexx = [];% zeros(3780);
% blockindexy = [];% zeros(3780);
logicalcellsindex = [];

% indexx = [];
% indexy = [];
imagevalues = zeros(rows, cols);

sigma = 1;
width = 6 * sigma + 1;
gfilter = fspecial('gaussian', [width width], sigma );

histogram = [];
anglevalues = 10:20:170;
%divide the image into 16x16 blocks with overlap
for i = 1 : 15
    for j = 1: 7
        % 1-16,8-24,16-32,24-32
        startxb = (i-1) * 8 + 1;
        startyb = (j-1) * 8 + 1;
        %get 16x16 pixels block
        blockm = magnitude(startxb : startxb + 15, startyb : startyb +15);
        blocko = orientation(startxb : startxb + 15, startyb : startyb +15);
        
        %apply gaussian of sigma = 8 on block magntiude
%         blockm = conv2(blockm, gfilter,'same');
        
        histogramb = [];
        %divide each block into 2x2 cells
        for a = 1 : 2
            for b= 1 : 2
                %1-8,9-16
                startxc = 8 * (a - 1) + 1;
                startyc = 8 * (b-1) + 1;
                cellm = blockm(startxc:startxc+7, startyc: startyc+7);
                cello = blocko(startxc:startxc+7, startyc: startyc+7);
                
                histogramc = zeros(1,9);
                for c = 1:8
                    for d = 1:8
                        m = cellm(c,d);
                        o = cello(c,d);
                        
                        for bin = 1: 9
                            if anglevalues(bin) > o
                                break
                            end
                        end
                        leftneighbour = bin - 1;
                        rightneighbour = bin;
                        if leftneighbour == 0
                            leftneighbour = 1;
                        elseif rightneighbour == 10
                            rightneighbour = 9;
                        end
                        
                        leftdistance = abs(o - anglevalues(leftneighbour));
                        rightdistance = abs(anglevalues(rightneighbour) - o);
                        total = leftdistance + rightdistance;
                        
%                         disp(['angle ' num2str(o) ' left neighbour = ' num2str(anglevalues(leftneighbour)) ' right distance = ' num2str(anglevalues(rightneighbour)) ' left distance = ' num2str(leftdistance) ' right distance = ' num2str(rightdistance)]);
                        
                        histogramc(1,leftneighbour) = histogramc(1,leftneighbour) + m * rightdistance / total;
                        histogramc(1,rightneighbour) = histogramc(1,rightneighbour) + m * leftdistance / total; 
                        
                        %histogramc(1,bin) = histogramc(1,bin) + m;
                        
                        %add the index in cells and blocks array
                        diff = 2;
                        x = startxb + startxc + c - diff;
%                         disp(x);
                        y = startyb + startyc + d - diff;
%                         disp(y);
                        %                         cellindexx(x, y) = startxc;
                        %                         cellindexy(x, y) = startyc;
                        %                         blockindexx(x, y) = startxb;
                        %                         blockindexy(x, y) = startyb;
                        
%                         cellindexx = [cellindexx  ];
%                         blockindexx = [blockindexx ((startxb-1) * 7 + 1)];
%                         blockindexy = [blockindexy ((startyb-1) * 7 + 1)];
                        
%                         indexx = [indexx x];
%                         indexy = [indexy y];
                        imagevalues(x,y) = m;
                    end
                end
                
                %add logical cell index
                %shift from end of the cell to the start
                x = startxb + startxc - 1;
                y = startyb + startyc - 1;
                x2 = (x-1) / 8 + 1;
                y2 = (y-1) / 8 + 1;
                logicalindex = (x2 - 1) * 8 + (y2); 
                logicalcellsindex = [logicalcellsindex logicalindex];
                
%                 disp(['point (' num2str(x) ' , ' num2str(y) ') -> new cell index = (' num2str(x2) ' , ' num2str(y2) ') concatenated index = ' num2str(logicalindex)]);
                

                %l2-hys normalization of block histogram
%                 l2norm = norm(histogramb);
%                 histogramb = histogramb / l2norm;
%                 %clip max to 0.2
%                 histogramb(find (histogramb > 0.2) ) = 0.2;
%                 l2norm = norm(histogramb, 1);
%                 histogramb = histogramb / l2norm;
                
                histogramb = [histogramb histogramc];
                
            end
        end
        
        histogram = [histogram histogramb];
        
    end
    
end



% figure;
% imshow(uint8(imagevalues));


% blockindexx(1:300)
% blockindexy

end


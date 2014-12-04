close all
clear all

%image1 = imread('butterfly2.jpg');
image = imread('Yosemite1.jpg');
[~,~,channels] = size(image);
if channels == 3
    image = rgb2gray(image);
end
figure(1);
imshow(image);
[rows, cols] = size(image);

%% step 1 - scale space peak selection

%parameters of the scale space
octavecount = 1;
sigmavalue=1.6;%sqrt(2.0)%1.6;07;
kvalue=sqrt(2.0);

%number of scales per octave
minlevel=1;
maxlevel=5;
length = maxlevel - minlevel;

gaussianvalues1 = zeros(rows, cols, length + 1);
gaussianvalues2 = zeros(ceil(rows/2), ceil(cols/2), length + 1);
dogvalues1 = zeros(rows,cols,length);
dogvalues2 = zeros(ceil(rows/2), ceil(cols/2), length );

%stores the potential key points
locations = zeros(rows, cols, length, octavecount);
count=zeros(length, octavecount);

%consider all the octaves
for octave = 1:octavecount
    
    
    %do calculation for first octave
    if octave == 1
        %find n+1 gaussians and n dog by taking their difference
        %first level - initilize the gaussian
        sigma = sigmavalue;
        width = floor(6 * sigma + 1);
        gfilter = fspecial('gaussian', width, sigma);
        gaussianvalues1(:,:,minlevel,octave) = conv2(double(image), gfilter, 'same');
        
        %for every next level, calculate the gaussian and dog values
        for level = minlevel+1: maxlevel
            sigma = sigmavalue * (kvalue ^ (level));
            width = floor(6 * sigma + 1);
            gfilter = fspecial('gaussian', width, sigma);
            gaussianvalues1(:,:,level,octave) = conv2(double(image), gfilter, 'same');
            
            %             figure(i+1);
            dogvalues1(:,:,level-1,octave) =  gaussianvalues1(:,:,level,octave) -  gaussianvalues1(:,:,level-1,octave);
            %             imshow( dogvalues1(:,:,level-1,octave));
        end
        
        % find the extremum points in the center levels (all levels exclusing the first and last)
        for centerlevel = minlevel+1 : maxlevel-2
            
            for i = 2: rows - 1
                for j = 2: cols -1
                    % find if the centerlevel pixel is the local maxima or minima
                    centerpixel  = dogvalues1(i,j,centerlevel, octave);
                    %find all the neighbouring 26 points.
                    neighbours = dogvalues1(i-1:i+1, j-1:j+1, centerlevel-1:centerlevel+1, octave);
                    neighbours= neighbours(:);
                    neighbours(14) = [];
                    
                    %find if greater than all values or smaller than all values
                    maxvalue = max(neighbours);
                    minvalue = min(neighbours);
                    
                    % store all the points in location so that they can be
                    % used later
                    if (centerpixel > maxvalue)
                        locations(i,j, centerlevel, octave) = 1;
                        count(centerlevel, octave) = count(centerlevel, octave) + 1;
                    elseif (centerpixel < minvalue)
                        locations(i,j, centerlevel, octave) = 1;
                        count(centerlevel, octave) = count(centerlevel, octave) + 1;
                    end
                end
            end
        end
        
        %do calculation for second octave
    else
        %find gaussian, downsample it for first scale
        sigma = sigmavalue * (kvalue ^ (level+2));
        width = floor(6 * sigma + 1);
        gfilter = fspecial('gaussian', width, sigma);
        gimage= conv2(double(image), gfilter, 'same');
        gaussianvalues2(:,:,minlevel,octave) = gimage(1:2:end, 1:2:end);
        
        %for next levels, find gaussian, down sample it then find dog
        for level = minlevel+1: maxlevel
            lastgaussianvalues = gaussianvalues1(:,:,level,octave-1);
            gaussianvalues2(:,:,level,octave) = lastgaussianvalues(1:2:end, 1:2:end);
            
            figure(i+1);
            dogvalues2(:,:,level-1, octave) =  gaussianvalues2(:,:,level,octave) -  gaussianvalues2(:,:,level-1,octave);
            imshow( dogvalues2(:,:,level-1, octave));
        end
        
        
        %for every next level, calculate the gaussian and dog values
        for centerlevel = minlevel+1 : maxlevel-2
            
            for i = 2: ceil(rows/2) - 1
                for j = 2: ceil(cols/2) -1
                    
                    %find if the centerlevel pixel is the local maxima or minima
                    centerpixel  = dogvalues2(i,j,centerlevel, octave);
                    %find all the neighbouring 26 points.
                    neighbours = dogvalues2(i-1:i+1, j-1:j+1, centerlevel-1:centerlevel+1, octave);
                    neighbours= neighbours(:);
                    neighbours(14) = [];
                    
                    %find if greater than all values or smaller than all values
                    maxvalue = max(neighbours);
                    minvalue = min(neighbours);
                    
                    % store all the points in location so that they can be
                    % used later
                    if (centerpixel > maxvalue)
                        locations(2*i,2*j, centerlevel, octave) = 1;
                        count(centerlevel, octave) = count(centerlevel, octave) + 1;
                    elseif (centerpixel < minvalue)
                        locations(2*i,2*j, centerlevel, octave) = 1;
                        count(centerlevel, octave) = count(centerlevel, octave) + 1;
                    end
                end
            end
        end
        
        
    end
    
    
end


%% step - display all the potential points
for octave = 1 : octavecount
    for centerlevel = minlevel + 1 : maxlevel - 2
        
        %find the points at that octave and scale
        [r,c] = find(locations(:,:,centerlevel, octave));                  % Find row,col coords.
        
        
        %display the found points
        disp(['total points found at ' num2str(centerlevel) ]);
        disp(size(r))
        
        % overlay corners on original image
        figure,
        imagesc(image),
        axis image,
        colormap(gray),
        hold on;
        %plot(c,r,'ys'),
        title(['potential points detected', num2str(octave), ' - ', num2str(centerlevel) ]);
        
        [length,~] = size(r);
        sigma = sigmavalue * (kvalue ^ centerlevel);
        radius = sigma;
        for pointindex = 1 : length
            y = r(pointindex);
            x = c(pointindex);
            ang=0:0.01:2*pi;
            xp=radius*cos(ang);
            yp=radius*sin(ang);
            plot(x+xp,y+yp);
        end
        
    end
end

%% step 2 - key point localization

filterdx = [-1 0 1; -1 0 1; -1 0 1];   % Derivative masks
filterdy = filterdx';
finalpositions1 = zeros(rows, cols, length, octavecount);
finalpositions = zeros(rows, cols, length, octavecount);

for octave = 1 : octavecount
    for centerlevel = minlevel + 1 : maxlevel - 2
        
        %find the points at that octave and scale
        [r,c] = find(locations(:,:,centerlevel, octave));                  % Find row,col coords.
        
        
        % 1 . initial outlier rejection
        %taylor series to find the proper location
        
        
        if octave == 1
            %find first derivative and second derivative at those points
                        Dx = conv2(dogvalues1(:,:,centerlevel,octave), filterdx, 'same');
                        Dxx = conv2(Dx, filterdx, 'same');
                        Dxy = conv2(Dx, filterdy, 'same');
            %             dx1 =  conv2(dogvalues1(:,:,centerlevel+1,octave), filterdx, 'same');
            %             dxsigma = dx1 - dx;
                        Dy = conv2(dogvalues1(:,:,centerlevel,octave), filterdy, 'same');
                        Dyx = conv2(Dy, filterdx, 'same');
                        Dyy = conv2(Dy, filterdy, 'same');
            %             dy1 = conv2(dogvalues1(:,:,centerlevel+1,octave), filterdy, 'same');
            %             dysigma = dy1 - dy;
            %             dsigma = dogvalues1(:,:,centerlevel+1,octave) - dogvalues1(:,:,centerlevel,octave);
            %             dsigmax = conv2(dsigma, filterdx, 'same');
            %             dsigmay = conv2(dsigma, filterdy, 'same');
            %             dsigma1 = dogvalues1(:,:,centerlevel+1,octave) - dogvalues1(:,:,centerlevel,octave);
            %             dsigmasigma = dsigma1 - dsigma;
            
            
            % find the taylor series values for the potential points
            taylorseries = zeros(rows, cols);
            
            [length1,~] = size(r);
            for pointindex = 1 : length1
                locx = r(pointindex,1);
                locy = c(pointindex,1);
                propercenterlevel = centerlevel;
                
                %find the extremum offset
                
                %find first derivative
                %                 dX = [dx(i,j)
                %                     dy(i,j)
                %                     dsigma(i,j)
                %                     ];
                dx = (dogvalues1(locx+1,locy,centerlevel,octave) - dogvalues1(locx-1,locy,centerlevel,octave))/2;
                dy =  (dogvalues1(locx,locy+1,centerlevel,octave) - dogvalues1(locx,locy-1,centerlevel,octave))/2;
                dsigma = (dogvalues1(locx,locy,centerlevel+1,octave) - dogvalues1(locx,locy,centerlevel-1,octave))/2;
                dX = [dx
                    dy
                    dsigma
                    ];
                
                %find hessian
                dx1 = (dogvalues1(locx+1,locy,centerlevel,octave) - dogvalues1(locx,locy,centerlevel,octave));
                dx2 = (dogvalues1(locx,locy,centerlevel,octave) - dogvalues1(locx-1,locy,centerlevel,octave));
                dx3 = (dogvalues1(locx+1,locy+1,centerlevel,octave) - dogvalues1(locx,locy+1,centerlevel,octave));
                dx4 = (dogvalues1(locx+1,locy,centerlevel+1,octave) - dogvalues1(locx,locy,centerlevel+1,octave));
                dxx = dx1-dx2;
                dxy = dx3 - dx1;
                dxsigma = dx4 - dx1;
                dy1 = (dogvalues1(locx,locy+1,centerlevel,octave) - dogvalues1(locx,locy,centerlevel,octave));
                dy2 = (dogvalues1(locx,locy,centerlevel,octave) - dogvalues1(locx,locy-1,centerlevel,octave));
                dy3 = (dogvalues1(locx+1,locy+1,centerlevel,octave) - dogvalues1(locx+1,locy,centerlevel,octave));
                dy4 = (dogvalues1(locx,locy+1,centerlevel+1,octave) - dogvalues1(locx,locy,centerlevel+1,octave));
                dyy = dy1-dy2;
                dyx = dy3 - dy1;
                dysigma = dy4 - dy2;
                dsigma1 = (dogvalues1(locx,locy,centerlevel+1,octave) - dogvalues1(locx,locy,centerlevel,octave));
                dsigma2 = (dogvalues1(locx,locy,centerlevel,octave) - dogvalues1(locx,locy,centerlevel-1,octave));
                dsigma3 = (dogvalues1(locx+1,locy,centerlevel+1,octave) - dogvalues1(locx+1,locy,centerlevel,octave));
                dsigma4 = (dogvalues1(locx,locy+1,centerlevel+1,octave) - dogvalues1(locx,locy+1,centerlevel,octave));
                dsigmasigma = dsigma1 - dsigma2;
                dsigmax = dsigma3 - dsigma1;
                dsigmay = dsigma4 - dsigma1;
                
                %                 hessian = [dxx(locx,j) dxy(locx,j) dxsigma(i,j)
                %                     dyx(i,j) dyy(i,j) dysigma(i,j)
                %                     dsigmax(i,j) dsigmay(i,j) dsigmasigma(i,j)];
                
                hessian = [dxx dxy dxsigma
                    dyx dyy dysigma
                    dsigmax dsigmay dsigmasigma];
                
                %find extremum offset
                extremum = - hessian'* dX;
                
                %add the extremum if the values are less than 1
                limit = 2;
                
%                 finalpositions1(locx,locy,centerlevel,octave) = 1;
                
                if abs(extremum(1)) <= limit && abs(extremum(2)) <= limit && abs(extremum(3)) <= 1
%                     disp('found extremum');
%                     disp(extremum);
                    
                    % add the offset by rounding
                    locx = locx + round(extremum(1));
                    locy = locy + round(extremum(2));
                    propercenterlevel = propercenterlevel + round(extremum(3));
                    
%                     disp('new extremum');
%                     disp(locx);
%                     disp(locy);
%                     disp(propercenterlevel);
                    
                    %                     if abs(extremum(1)) > 0.5 || abs(extremum(2)) > 0.5 || abs(extremum(3)) > 0.5
                    %                         %disp(['hello: ', num2str(extremum(1)), ',' , num2str(extremum(2)) , ',' num2str(extremum(3))]);
                    %                     end
                    %
                    %                     if abs(extremum(1)) < 0.5 && abs(extremum(2)) < 0.5 && abs(extremum(3)) < 0.5
                    %
                    %                         disp('hello');
                    %                     end
                    
                    taylorseries(locx,locy) = dogvalues1(locx,locy,propercenterlevel,octave)+  0.5 * ( extremum(1)*dX(1) + extremum(2)*dX(2) + extremum(3)* dX(3)) ;
                    
%                     disp(['tailor series ' num2str(taylorseries(locx,locy))] );
                    finalpositions1(locx,locy,propercenterlevel,octave) = 1;
                    %remove all the points with less contrast
                    if abs( taylorseries(locx,locy)) > 0.03;
%                          finalpositions(locx,locy,propercenterlevel,octave) = 1;
                    end
                end
            end
            
            
            
            
            % 2. further outlier rejection
            %remove all the edges
            %use harris corner detector
            [r1,c1] = find(finalpositions1(:,:,centerlevel, octave));
            [length1,~] = size(r1);
            
            for k = 1 : length1
                locx = r1(k,1);
                locy = c1(k,1);
                
                dx1 = (dogvalues1(locx+1,locy,centerlevel,octave) - dogvalues1(locx,locy,centerlevel,octave));
                dx2 = (dogvalues1(locx,locy,centerlevel,octave) - dogvalues1(locx-1,locy,centerlevel,octave));
                dx3 = (dogvalues1(locx+1,locy+1,centerlevel,octave) - dogvalues1(locx,locy+1,centerlevel,octave));
                dxx = dx1-dx2;
                dxy = dx3 - dx1;
                dy1 = (dogvalues1(locx,locy+1,centerlevel,octave) - dogvalues1(locx,locy,centerlevel,octave));
                dy2 = (dogvalues1(locx,locy,centerlevel,octave) - dogvalues1(locx,locy-1,centerlevel,octave));
                dyy = dy1-dy2;
                
%                 hessian = [dxx dxy
%                     dyx dyy];
                hessian = [Dxx(locx,locy) Dxy(locx,locy)
                    Dyx(locx,locy) Dyy(locx,locy)];
                eigenvalues = eig(hessian);
                if eigenvalues(2) ~= 0
                    ratio = eigenvalues(1)/eigenvalues(2);
                    if ratio < 0
                        ratio = -ratio;
                    end
                    if abs(ratio) <= 10
                        %disp('hello2');
                        finalpositions(locx,locy, centerlevel ,octave)=1;
                    end
                end
            end
            
        end
    end
    
end

%% step - display all the potential points
for octave = 1 : octavecount
    for centerlevel = minlevel + 1 : maxlevel - 2
        
        %find the points at that octave and scale
        [r,c] = find(finalpositions(:,:,centerlevel, octave));                  % Find row,col coords.
        
        
        %display the found points
        disp(['total points found at ' num2str(centerlevel) ]);
        disp(size(r))
        
        % overlay corners on original image
        figure,
        imagesc(image),
        axis image,
        colormap(gray),
        hold on;
        %plot(c,r,'ys'),
        title(['potential points detected', num2str(octave), ' - ', num2str(centerlevel) ]);
        
        [length,~] = size(r);
        sigma = sigmavalue * (kvalue ^ centerlevel);
        radius = sigma;
        for pointindex = 1 : length
            y = r(pointindex);
            x = c(pointindex);
            ang=0:0.01:2*pi;
            xp=radius*cos(ang);
            yp=radius*sin(ang);
            plot(x+xp,y+yp);
        end
        
    end
end

% %% step show all the finalized key points
% for octave = 1 : octavecount
%     for centerlevel = minlevel + 1 : maxlevel - 2
%         
%         %find the points at that octave and scale
%         [r,c] = find(finalpositions(:,:,centerlevel, octave));                  % Find row,col coords.
%         disp(['total points found at '  num2str(centerlevel) ' = ']);
%         disp(size(r))
%         
%         % overlay corners on original image
%         figure;
%         imagesc(image),
%         axis image,
%         colormap(gray),
%         hold on;
%         %plot(c,r,'ys'),
%         title('key points detected');
%         
%         [length,~] = size(r);
%         sigma = sigmavalue * (kvalue ^ centerlevel);
%         radius = sigma;
%         for pointindex = 1 : length
%             y = r(pointindex);
%             x = c(pointindex);
%             ang=0:0.01:2*pi;
%             xp=radius*cos(ang);
%             yp=radius*sin(ang);
%             plot(x+xp,y+yp);
%         end
%         
%     end
%     
% end

%
%
% %% step 3 - orientation assignment
%
% featurelocations = [];
% featuredescriptions = [];
%
% for octave = 1 : octavecount
%     for centerlevel = minlevel + 1 : maxlevel - 2
%
%         if octave == 1
%             [r,c] = find(locations(:,:,centerlevel, octave));                  % Find row,col coords.
%             [length1,~] = size(r);
%             %do edge outliers
%             for pointindex = 1 : length1
%                 locx = r(pointindex);
%                 locy = c(pointindex);
%
%                 %store all the locations
%                 featurelocations = [featurelocations ; [locx, locy, centerlevel, octave] ];
%
%                 magnitudes = [];
%                 orientations = [];
%
%                 %find orientation and magnitude in the 4x4 neighbourhood
%                 for m = -1 : 2
%
%                     magnitudes1 = [];
%                     orientations1 = [];
%                     for n = -1 : 2
%
%                         i = locx + m;
%                         j = locy + n;
%
%                         %calculate magnitude and orientation
%                         dx = gaussianvalues1(i+1,j,centerlevel)-gaussianvalues1(i-1,j,centerlevel);
%                         dy = gaussianvalues1(i,j+1,centerlevel)-gaussianvalues1(i,j-1,centerlevel);
%                         magnitude = sqrt(dx*dx + dy*dy);
%                         magnitudes1 =[magnitudes1 ;magnitude];
%
%                         orientation = atand(dy/dx);
%                         %correct the orientation from 0-180 to 0-360 range
%                         if orientation > 0 && orientation < 90
%                             if gaussianvalues1(i-1,j,centerlevel) > gaussianvalues1(i+1,j,centerlevel)
%                                 orientation = orientation - 180;
%                             end
%                         elseif orientation <= 0 && orientation > -90
%                             if gaussianvalues1(i-1,j,centerlevel) > gaussianvalues1(i+1,j,centerlevel)
%                                 orientation = orientation + 180;
%                             end
%                         end
%                         orientations1 = [orientations1 orientation];
%
%                     end
%
%                     magnitudes =[magnitudes ;magnitudes1];
%                     orientations = [orientations orientations1];
%                 end
%
%                 %smooth the magnitude values with gaussian of sigma
%                 %with 1.5 sigma
%                 sigma = 1.5* sigmavalue * (kvalue ^ locations(a,3));
%                 width = ceil(6 * sigma + 1);
%                 gfilter = fspecial('gaussian', [width width], sigma);
%                 magnitudes = conv2(magnitudes, gfilter, 'same');
%
%                 %find histogram of orientation in neighbourhood
%                 histogram = zeros(1,36);
%                 bins = 1:36;
%                 [length7,length2] = size(magnitudes);
%
%                 for a = 1 : length7
%                     for b = 1 : length2
%                         bin = ceil(orientations(a,b)/10);
%                         histogram(1,bin ) =  histogram(1,hindex) + magnitudes(a,b);
%                     end
%                 end
%                 % Example Visualise:
%                 %figure;
%                 %bar(bins, histogram)
%
%                 %find the dominant directions in the histogram
%                 combined = [histogram; 1:36]';
%                 sortedvalues = sortrows(combined, 1);
%                 maxvalue = [sortedvalues(36,1)];
%                 domiantdirections = [sortedvalues(36,2)];
%                 for m = 35 : -1 : 1
%                     if sortedvalues(m) > 0.8 * maxvalue(1)
%                         % maxvalue = [maxvalue sortedvalues(m,1)];
%                         domiantdirections = [domiantdirections sortedvalues(m,2)];
%                     end
%                 end
%
%                 disp(domiantdirections);
%
%                 %% step 4 - keypoint descriptor
%
%                 %find the feature descriptor for 16x16 block
%                 featuredescriptor = [];
%                 %first calculate the magnitude and relative orientation for
%                 %all 16x16 values
%
%                 % determine the boundary, correct if it goes over the image
%                 % boundaries
%                 m1 = locx - 8;
%                 m2 = locx + 7;
%                 n1 = locy - 8;
%                 n2 = locy + 7;
%                 if m1 < 2
%                     delta = 2 - m1;
%                     m2 = m2 + delta;
%                     m1 = 2;
%                 elseif m2 > rows - 2
%                     delta = m2 - rows + 1;
%                     m1 = m1 - delta;
%                     m2 = rows - 1;
%                 end
%                 if n1 < 2
%                     delta = 2 - n1;
%                     n2 = n2 + delta;
%                     n1 = 2;
%                 elseif n2 > rows - 2
%                     delta = n2 - rows + 1;
%                     n1 = n1 - delta;
%                     n2 = rows - 1;
%                 end
%
%                 %disp(['boundaries ' , num2str(m1), ' , ', num2str(m2), ' , ' , num2str(n1) , ' , ' num2str(n2)]);
%
%                 %calculate for 4x4 blocks
%                 for m = 1: 4
%                     for n = 1:4
%
%                         mvalue = m1 + (4 * (m-1));
%                         nvalue = n1 + (4 * (n-1));
%
%                         magnitudes = [];
%                         orientations = [];
%
%                         %calculate for 4x4 pixels
%                         for a = 1: 4
%
%                             magnitudes1 =[];
%                             orientations1= [];
%
%                             for b = 1:4
%
%                                 i = mvalue + a - 1;
%                                 j = nvalue + b - 1;
%
%                                 %calculate magnitude and orientation
%                                 dx = gaussianvalues1(i+1,j,centerlevel)-gaussianvalues1(i-1,j,centerlevel);
%                                 dy = gaussianvalues1(i,j+1,centerlevel)-gaussianvalues1(i,j-1,centerlevel);
%                                 magnitude = sqrt(dx*dx + dy*dy);
%                                 magnitudes1 =[magnitudes1 ;magnitude];
%
%                                 orientation = atand(dy/dx);
%                                 %correct the orientation from 0-180 to 0-360 range
%                                 if orientation > 0 && orientation < 90
%                                     if gaussianvalues1(i-1,j,centerlevel) > gaussianvalues1(i+1,j,centerlevel)
%                                         orientation = orientation - 180;
%                                     end
%                                 elseif orientation <= 0 && orientation > -90
%                                     if gaussianvalues1(i-1,j,centerlevel) > gaussianvalues1(i+1,j,centerlevel)
%                                         orientation = orientation + 180;
%                                     end
%                                 end
%                                 orientation = orientation - domiantdirections(1);
%
%                                 orientations = [orientations orientation];
%                                 orientations1 =[orientations1 ;orientations];
%                             end
%
%                             magnitudes =[magnitudes ;magnitudes1];
%                             orientations =[orientations ;orientations1];
%
%                         end
%
%
%                         %smooth the magnitude values with gaussian of sigma
%                         %with 1.5 sigma
%                         sigma = 1.5* sigmavalue * (kvalue ^ locations(a,3));
%                         width = ceil(6 * sigma + 1);
%                         gfilter = fspecial('gaussian', [width width], sigma);
%                         magnitudes = conv2(magnitudes, gfilter, 'same');
%
%                         %form 8 bin histogram
%                         histogram = zeros(1,8);
%                         bins = 1:8;
%                         [length7,length2] = size(magnitudes);
%
%                         for a = 1 : length7
%                             for b = 1 : length2
%                                 bin = ceil(orientations(a,b)/45);
%                                 histogram(1,bin ) =  histogram(1,hindex) + magnitudes(a,b);
%                             end
%                         end
%
%
%                         %append to feature descriptor, to get the
%                         %descriptor of one point
%                         featuredescriptor = [featuredescriptor histogram];
%
%
%                     end
%                 end
%
%                 %append to the list of descriptors for all the points
%                 %disp(featuredescriptor);
%                 featuredescriptions= [featuredescriptions ; [featuredescriptor]];
%
%
%             end
%
%
%         end
%     end
% end
%
% %% step 4 - (a) - post processing
% %normalize to unit vector
% % remove any element greater than 0.2
% % again renormalize
%
% [length6, ~] = size(featuredescriptions);
% for i = 1 : length6
%     featuredescriptions(i,:) =  featuredescriptions(i,:)/norm(featuredescriptions(i,:));
%
%     featuredescriptions(i,find(featuredescriptions(i,:)) > 0.2) = 0;
%
%     featuredescriptions(i,:) =  featuredescriptions(i,:)/norm(featuredescriptions(i,:));
% end




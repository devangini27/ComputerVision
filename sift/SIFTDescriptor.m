function [ output_args ] = SIFTDescriptor( image )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%image1 = imread('butterfly2.jpg');

figure(1);
imshow(image);
[rows, cols] = size(image);

octavecount = 2;
sigmavalue=1.6;%sqrt(2.0)%1.6;
kvalue=sqrt(2.0);

minlevel=1;
maxlevel=5;
length = maxlevel - minlevel;
gaussianvalues1 = zeros(rows, cols, length + 1);
gaussianvalues2 = zeros(ceil(rows/2), ceil(cols/2), length + 1);
dogvalues1 = zeros(rows,cols,length);
dogvalues2 = zeros(ceil(rows/2), ceil(cols/2), length );

centerlevel = (minlevel+maxlevel)/2;
locations = zeros(rows, cols, length, octavecount);
count=zeros(length, octavecount);


for octave = 1:octavecount
    
    
    % find all the gaussian and dog for the octave
    
    if octave == 1
        %perform gaussian for first level and find dog
        sigma = sigmavalue;
        width = floor(6 * sigma + 1);
        log = fspecial('gaussian', width, sigma);
        gaussianvalues1(:,:,minlevel,octave) = conv2(double(image), log, 'same');
        
        for level = minlevel+1: maxlevel
            sigma = sigmavalue * (kvalue ^ (level));
            width = floor(6 * sigma + 1);
            log = fspecial('gaussian', width, sigma);
            gaussianvalues1(:,:,level,octave) = conv2(double(image), log, 'same');
            
%             figure(i+1);
            dogvalues1(:,:,level-1,octave) =  gaussianvalues1(:,:,level,octave) -  gaussianvalues1(:,:,level-1,octave);
%             imshow( dogvalues1(:,:,level-1,octave));
            
            %logvalues = [logvalues; [smoothedvalues]];
        end
        
        for centerlevel = minlevel+1 : maxlevel-2
            % find the extremas
            for i = 2: rows - 1
                for j = 2: cols -1
                    
                    %find if the centerlevel pixel is the local maxima or minima
                    centerpixel  = dogvalues1(i,j,centerlevel, octave);
                    neighbours = dogvalues1(i-1:i+1, j-1:j+1, centerlevel-1:centerlevel+1, octave);
                    neighbours= neighbours(:);
                    
                    %find if greater than all values or smaller than all values
                    maxvalue = max(neighbours);
                    minvalue = min(neighbours);
                    
                    if (centerpixel == maxvalue)
                        locations(i,j, centerlevel, octave) = 1;
                        count(centerlevel, octave) = count(centerlevel, octave) + 1;
                    elseif (centerpixel == minvalue)
                        locations(i,j, centerlevel, octave) = 1;
                        count(centerlevel, octave) = count(centerlevel, octave) + 1;
                    end
                end
            end
        end
        
    else
        %else downsample last level gaussian values and find dog
        %down sample the image
        %Xdown = X(1:2:end,1:2:end);
        
        sigma = sigmavalue * (kvalue ^ (level+2));
        width = floor(6 * sigma + 1);
        gfilter = fspecial('gaussian', width, sigma);
        gimage= conv2(double(image), log, 'same');
        
        %lastgaussianvalues = gaussianvalues(:,:,minlevel,octave-1);
        gaussianvalues2(:,:,minlevel,octave) = gimage(1:2:end, 1:2:end);
        
        for level = minlevel+1: maxlevel
            lastgaussianvalues = gaussianvalues1(:,:,level,octave-1);
            gaussianvalues2(:,:,level,octave) = lastgaussianvalues(1:2:end, 1:2:end);
            
            figure(i+1);
            dogvalues2(:,:,level-1, octave) =  gaussianvalues2(:,:,level,octave) -  gaussianvalues2(:,:,level-1,octave);
            imshow( dogvalues2(:,:,level-1, octave));
            
            %logvalues = [logvalues; [smoothedvalues]];
        end
        
        
        
        for centerlevel = minlevel+1 : maxlevel-2
            % find the extremas
            for i = 2: ceil(rows/2) - 1
                for j = 2: ceil(cols/2) -1
                    
                    %find if the centerlevel pixel is the local maxima or minima
                    centerpixel  = dogvalues2(i,j,centerlevel, octave);
                    neighbours = dogvalues2(i-1:i+1, j-1:j+1, centerlevel-1:centerlevel+1, octave);
                    neighbours= neighbours(:);
                    
                    %find if greater than all values or smaller than all values
                    maxvalue = max(neighbours);
                    minvalue = min(neighbours);
                    
                    if (centerpixel == maxvalue)
                        locations(2*i,2*j, centerlevel, octave) = 1;
                        count(centerlevel, octave) = count(centerlevel, octave) + 1;
                    elseif (centerpixel == minvalue)
                        locations(2*i,2*j, centerlevel, octave) = 1;
                        count(centerlevel, octave) = count(centerlevel, octave) + 1;
                    end
                end
            end
        end
        
        
    end
    
    
    
    
end

%taylorseries = zeros(rows, cols, maxlevel-minlevel-2);
filterdx = [-1 0 1; -1 0 1; -1 0 1];   % Derivative masks
filterdy = filterdx';
finalpositions = zeros(rows, cols, length, octave);


for octave = 1 : octavecount
    for centerlevel = minlevel + 1 : maxlevel - 2
        [r,c] = find(locations(:,:,centerlevel, octave));                  % Find row,col coords.
        
        
        disp(['total points found at ' num2str(centerlevel) ]);
        disp(size(r))
        
        % overlay corners on original image
        figure(maxlevel+1 + centerlevel + octave * 2),
        imagesc(image),
        axis image,
        colormap(gray),
        hold on;
        %plot(c,r,'ys'),
        title(['corners detected', num2str(octave), ' - ', num2str(centerlevel) ]);
        
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
        
        %         [length1,~] = size(r)
        %         for k = 1 : length1
        %             centx = r(k);
        %             centy = c(k);
        %             sigma = sigmavalue * (kvalue ^ (level+2));
        %             width = sigma%floor(3* sigma);
        %             radius = sigma;
        %             theta = 0 : (2 * pi / 10000) : (2 * pi);
        %             pline_x = radius * cos(theta) + centx;
        %             pline_y = radius * sin(theta) + centy;
        %             %k = ishold;
        %             plot(pline_x, pline_y, 'r-', 'LineWidth', 1);
        %         end
        
        
        
        %do outlier detection
        %find first derivative and second derivative at those points
        
        if octave == 1
            dx = conv2(dogvalues1(:,:,centerlevel-1,octave), filterdx, 'same');
            dxx = conv2(dx, filterdx, 'same');
            dxy = conv2(dx, filterdy, 'same');
            dx1 =  conv2(dogvalues1(:,:,centerlevel+1,octave), filterdx, 'same');
            dxsigma = dx1 - dx;
            dy = conv2(dogvalues1(:,:,centerlevel-1,octave), filterdy, 'same');
            dyx = conv2(dy, filterdx, 'same');
            dyy = conv2(dy, filterdy, 'same');
            dy1 = conv2(dogvalues1(:,:,centerlevel+1,octave), filterdy, 'same');
            dysigma = dy1 - dy;
            dsigma = dogvalues1(:,:,centerlevel+1,octave) - dogvalues1(:,:,centerlevel,octave);
            dsigmax = conv2(dsigma, filterdx, 'same');
            dsigmay = conv2(dsigma, filterdy, 'same');
            dsigma1 = dogvalues1(:,:,centerlevel+1,octave) - dogvalues1(:,:,centerlevel,octave);
            dsigmasigma = dsigma1 - dsigma;
            
            
            
            taylorseries = zeros(rows, cols);
            
            [length1,~] = size(r);
            for k = 1 : length1
                i = r(k,1);
                j = c(k,1);
                
                %find the extremum offset
                %find first derivative
                dX = [dx(i,j)
                    dy(i,j)
                    dsigma(i,j)
                    ];
                
                %                 dX = [(dogvalues1(i+1,j,centerlevel,octave) - dogvalues1(i-1,j,centerlevel,octave))/2
                %                     (dogvalues1(i,j+1,centerlevel,octave) - dogvalues1(i,j-1,centerlevel,octave))/2
                %                     (dogvalues1(i,j,centerlevel+1,octave) - dogvalues1(i,j,centerlevel-1,octave))/2];
                
                %find hessian
                hessian = [dxx(i,j) dxy(i,j) dxsigma(i,j)
                    dyx(i,j) dyy(i,j) dysigma(i,j)
                    dsigmax(i,j) dsigmay(i,j) dsigmasigma(i,j)];
                
                %find extremum offset
                extremum = - hessian * dX;
                
                if abs(extremum(1)) > 0.5 || abs(extremum(2)) > 0.5 || abs(extremum(3)) > 0.5
                    %disp(['hello: ', num2str(extremum(1)), ',' , num2str(extremum(2)) , ',' num2str(extremum(3))]);
                end
                
                if abs(extremum(1)) < 0.5 && abs(extremum(2)) < 0.5 && abs(extremum(3)) < 0.5
                    
                    disp('hello');
                    taylorseries(i,j) = dogvalues1(i,j,centerlevel,octave)+  0.5 * ( extremum(1)*dX(1) + extremum(2)*dX(2) + extremum(3)* dX(3)) ;
                end
            end
            %finalpositions(find(taylorseries > 0.03)) = 1;
            
            
            [length1,~] = size(r)
            %do edge outliers
            for k = 1 : length1
                i = r(k,1);
                j = c(k,1);
                
                hessian = [dxx(i,j) dxy(i,j)
                    dyx(i,j) dyy(i,j)];
                eigenvalues = eig(hessian);
                if eigenvalues(2) ~= 0
                    ratio = eigenvalues(1)/eigenvalues(2);
                    if ratio < 0
                        ratio = -ratio;
                    end
                    if abs(ratio) <= 10
                        %disp('hello2');
                        finalpositions(i,j, centerlevel ,octave)=1;
                    end
                end
            end
            
            
            
            [r,c] = find(finalpositions);                  % Find row,col coords.
            
            
            disp(['total points found = ' ]);
            disp(size(r))
            
            % overlay corners on original image
            figure;
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
        
    end
    
    
    
end



%calculate orientation

featurelocations = [];
featuredescriptions = [];
% magnitudes = [];
% orientations = [];

for octave = 1 : octavecount
    for centerlevel = minlevel + 1 : maxlevel - 2
        
        if octave == 1
            [r,c] = find(locations(:,:,centerlevel, octave));                  % Find row,col coords.
            [length1,~] = size(r);
         
            for k = 1 : length1
                locx = r(k);
                locy = c(k);
                
                featurelocations = [featurelocations ; [locx, locy, centerlevel, octave] ];
                
                
                %featurelocations = [];
                magnitudes = [];
                orientations = [];
                
                %find orientation and magnitude in the 4x4 neighbourhood
                for m = -1 : 2
                    for n = -1 : 2
                        
                        i = locx + m;
                        j = locy + n;
                        
                        %featurelocations = [featurelocations; [i, j, centerlevel]];
                        
                        %calculate magnitude and orientation
                        dx = gaussianvalues1(i+1,j,centerlevel)-gaussianvalues1(i-1,j,centerlevel);
                        dy = gaussianvalues1(i,j+1,centerlevel)-gaussianvalues1(i,j-1,centerlevel);
                        magnitude = sqrt(dx*dx + dy*dy);
                        magnitudes =[magnitudes ;magnitude];
                        
                        orientation = atand(dy/dx);
                        %correct the orientation from 0-180 to 0-360 range
                        if orientation > 0 && orientation < 90
                            if gaussianvalues1(i-1,j,centerlevel) > gaussianvalues1(i+1,j,centerlevel)
                                orientation = orientation - 180;
                            end
                        elseif orientation <= 0 && orientation > -90
                            if gaussianvalues1(i-1,j,centerlevel) > gaussianvalues1(i+1,j,centerlevel)
                                orientation = orientation + 180;
                            end
                        end
                        orientations = [orientations orientation];
                        
                    end
                end
                
                %find histogram of orientation in neighbourhood
                histogram = zeros(1,36);
                bins = 1:36;
                [~,length2] = size(magnitudes);
                for a = 1 : length2
                    bin = floor(orientations(a)/10);
                    sigma = sigmavalue * (kvalue ^ locations(a,3));
                    gausswindow = gausswin(floor(1.5 * sigma));
                    gausswindow = gausswindow/norm(gausswindow);
                    gaussianweight = magnitudes(a) * gausswindow;
                    
                    [length3 ,~] = size(gaussianweight);
                    
                    width = floor(length3)/2;
                    for b = bin - width : bin - width + length3 - 1
                        %disp(j);
                        hindex = b;
                        if b <=0
                            hindex = hindex + 36;
                        end
                        histogram(1,hindex) =  histogram(1,hindex) + gaussianweight(b - bin + width + 1);
                    end
                    
                end
                % Example Visualise:
                figure;
                bar(bins, histogram)
                
                %find the dominant directions in the histogram
                combined = [histogram; 1:36]';
                sortedvalues = sortrows(combined, 1);
                maxvalue = [sortedvalues(36,1)];
                domiantdirections = [sortedvalues(36,2)];
                for m = 35 : -1 : 1
                    if sortedvalues(m) > 0.8 * maxvalue(1)
                        % maxvalue = [maxvalue sortedvalues(m,1)];
                        domiantdirections = [domiantdirections sortedvalues(m,2)];
                    end
                end
                
                disp(domiantdirections);
                
                
                %find the feature descriptor for 16x16 block
                featuredescriptor = [];
                %first calculate the magnitude and relative orientation for
                %all 16x16 values
                
                % determine the boundary
                m1 = locx - 8;
                m2 = locx + 7;
                n1 = locy - 8;
                n2 = locy + 7;
                if m1 < 2
                    delta = 2 - m1;
                    m2 = m2 + delta;
                    m1 = 2;
                elseif m2 > rows - 2
                    delta = m2 - rows + 1;
                    m1 = m1 - delta;
                    m2 = rows - 1;
                end
                if n1 < 2
                    delta = 2 - n1;
                    n2 = n2 + delta;
                    n1 = 2;
                elseif n2 > rows - 2
                    delta = n2 - rows + 1;
                    n1 = n1 - delta;
                    n2 = rows - 1;
                end
                
                %disp(['boundaries ' , num2str(m1), ' , ', num2str(m2), ' , ' , num2str(n1) , ' , ' num2str(n2)]);
                
                for m = 1: 4
                    for n = 1:4
                        
                        mvalue = m1 + (4 * (m-1));
                        nvalue = n1 + (4 * (n-1));
                        
                        magnitudes = [];
                        orientations = [];
                        
                        %calculate for 4x4 block
                        for a = 1: 4
                            for b = 1:4
                                
                                i = mvalue + a - 1;
                                j = nvalue + b - 1;
                                
                                %calculate magnitude and orientation
                                dx = gaussianvalues1(i+1,j,centerlevel)-gaussianvalues1(i-1,j,centerlevel);
                                dy = gaussianvalues1(i,j+1,centerlevel)-gaussianvalues1(i,j-1,centerlevel);
                                magnitude = sqrt(dx*dx + dy*dy);
                                magnitudes =[magnitudes ;magnitude];
                                
                                orientation = atand(dy/dx);
                                %correct the orientation from 0-180 to 0-360 range
                                if orientation > 0 && orientation < 90
                                    if gaussianvalues1(i-1,j,centerlevel) > gaussianvalues1(i+1,j,centerlevel)
                                        orientation = orientation - 180;
                                    end
                                elseif orientation <= 0 && orientation > -90
                                    if gaussianvalues1(i-1,j,centerlevel) > gaussianvalues1(i+1,j,centerlevel)
                                        orientation = orientation + 180;
                                    end
                                end
                                orientation = orientation - domiantdirections(1);
                                
                                orientations = [orientations orientation];
                            end
                        end
                        
                        
                        %form 8 bin histogram
                        histogram = zeros(1,8);
                        bins = 1:8;
                        [~,length2] = size(magnitudes);
                        for a = 1 : length2
                            bin = floor(orientations(a)/45);
                            sigma = sigmavalue * (kvalue ^ locations(a,3));
                            gausswindow = gausswin(floor(1.5 * sigma));
                            gausswindow = gausswindow/norm(gausswindow);
                            gaussianweight = magnitudes(a) * gausswindow;
                            
                            [length3 ,~] = size(gaussianweight);
                            
                            width = floor(length3)/2;
                            for b = bin - width : bin - width + length3 - 1
                                %disp(j);
                                hindex = b;
                                if b <=0
                                    hindex = hindex + 8;
                                end
                                histogram(1,hindex) =  histogram(1,hindex) + gaussianweight(b - bin + width + 1);
                            end
                            
                        end
                        
                        %append to feature descriptor
                        featuredescriptor = [featuredescriptor histogram];
                        
                        
                        
                        
                        
                        
                    end
                end
                
                
                %disp(featuredescriptor);
                featuredescriptions= [featuredescriptions ; [featuredescriptor]];
                
                
            end
            
            
            
            
        end
    end
end

%post processing
%normalize to unit vector
% remove any element greater than 0.2
% again renormalize

[length6, ~] = size(featuredescriptions);
for i = 1 : length6
    featuredescriptions(i,:) =  featuredescriptions(i,:)/norm(featuredescriptions(i,:));
    
    featuredescriptions(i,find(featuredescriptions(i,:)) > 0.2) = 0;
    
    featuredescriptions(i,:) =  featuredescriptions(i,:)/norm(featuredescriptions(i,:));
end





end


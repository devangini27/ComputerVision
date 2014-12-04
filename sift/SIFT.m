classdef SIFT<handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        
        %% parameters of the scale space
        octavecount = 1;
        sigmavalue=1.6;%sqrt(2.0)%1.6;07;
        kvalue=sqrt(2.0);
        
        %number of scales per octave
        minlevel=1;
        maxlevel=5;
        length1 = 4;%maxlevel - minlevel;
        
        gaussianvalues1 = [];
        gaussianvalues2 = [];
        dogvalues1 = [];
        dogvalues2 = [];
        
        %% stores the potential key points
        locations = [];
        count = [];
        finalpositions = [];
        
        %% image properties
        image = [];
        rows = 0;
        cols = 0;
        
    end
    
    methods
        function obj = SIFT()
        end
        
        
        %% main function that is called
        
        function [ descriptors, locations] = findSIFTFeatures(self, imagename)
            
            
            
            self.image = imread(imagename);
            [~,~,channels] = size(self.image);
            if channels == 3
                self.image = rgb2gray(self.image);
            end
            figure(1);
            imshow(self.image);
            [self.rows, self.cols] = size(self.image);
            
            %initialize all the matrices
            self.gaussianvalues1 = zeros(self.rows, self.cols, self.length1 + 1);
            self.gaussianvalues2 = zeros(ceil(self.rows/2), ceil(self.cols/2), self.length1 + 1);
            self.dogvalues1 = zeros(self.rows,self.cols,self.length1);
            self.dogvalues2 = zeros(ceil(self.rows/2), ceil(self.cols/2), self.length1 );
            
            %stores the potential key points and counts
            self.locations = zeros(self.rows, self.cols, self.length1, self.octavecount);
            self.count=zeros(self.length1, self.octavecount);
            
            % start actual sift process
            
            %1 - scale space peak selection
            disp('performing scale space peak selection');
            self.performScaleSpacePeakSelection();
            self.displayKeypoints(self.locations);
            
            %2 - outlier rejection
            %             self.outlierRejection();
            disp('performing outlier rejection');
            self.harrisCornerRejection();
            self.displayKeypoints(self.finalpositions);
            
            % 3 -find dominant direction
            dominantDirections = self.calculateDominantOrientation();
            self.displayKeypointsDirections(self.finalpositions, dominantDirections);
            
            %4- create sift descriptors
            [descriptors, locations] = self.calculateSIFTDescriptor(dominantDirections);
            
            %5 - post processing
            descriptors = postprocessDescriptors(self, descriptors);
        end
        
        
        %% step 1 - scale space peak selection
        
        
        function performScaleSpacePeakSelection(self)
            
            %consider all the octaves
            for octave = 1 : self.octavecount
                
                if octave == 1
                    gaussianarray = self.gaussianvalues1;
                    dogarray = self.dogvalues1;
                    rowvalue = self.rows;
                    colvalue = self.cols;
                else
                    gaussianarray = self.gaussianvalues2;
                    dogarray = self.dogvalues2;
                    rowvalue = ceil(self.rows/2);
                    colvalue = ceil(self.cols/2);
                end
                
                %find n+1 gaussians and n dog by taking their difference
                %first level - initilize the gaussian
                gaussianarray(:,:,self.minlevel,octave) = self.computeGaussian(octave, 0);
                
                %for every next level, calculate the gaussian and dog values
                for level = self.minlevel: self.maxlevel-1
                    gaussianarray(:,:,level+1,octave) = self.computeGaussian(octave, level);
                    %             figure(level+1);
                    dogarray(:,:,level,octave) =  gaussianarray(:,:,level+1,octave) -  gaussianarray(:,:,level,octave);
                    %             imshow( dogvalues1(:,:,level-1,octave));
                end
                
                % find the extremum points in the center levels (all levels exclusing the first and last)
                for centerlevel = self.minlevel+1 : self.maxlevel-2
                    for i = 2: rowvalue - 1
                        for j = 2: colvalue -1
                            if self.findextrema(dogarray, i, j, centerlevel, octave) == 1
                                self.locations(i,j, centerlevel, octave) = 1;
                                self.count(centerlevel, octave) = self.count(centerlevel, octave) + 1;
                            end
                        end
                    end
                end
                
                %write the values back
                if octave == 1
                    self.gaussianvalues1 = gaussianarray;
                    self.dogvalues1 = dogarray;
                else
                    self.gaussianvalues2 = gaussianarray;
                    self.dogvalues2 = dogarray;
                end
                
            end
        end
        
        
        
        %% gaussian related functions
        
        function [width] = computeGaussianRadius(self, sigma)
            width = floor(6 * sigma + 1);
        end
        
        
        function [sigma] = computeGaussianSigma(self, octave, level)
            if octave == 1
                sigma = self.sigmavalue * (self.kvalue ^ (level));
            else
                %increase level by 2 for the second octave
                sigma = self.sigmavalue * (self.kvalue ^ (level+2));
            end
        end
        
        
        function [gimage] = computeGaussian(self, octave, level)
            
            sigma = self.computeGaussianSigma(octave, level);
            width = self.computeGaussianRadius(sigma);
            gfilter = fspecial('gaussian', width, sigma);
            gimage = conv2(double(self.image), gfilter, 'same');
            
            if octave == 2
                % find gaussian, downsample it for first scale
                gimage = gimage(1:2:end, 1:2:end);
            end
        end
        
        
        
        %% function to find extrema in 3x3x3 neigbhourhood
        
        function [truthvalue] = findextrema(self, array, i, j, centerlevel, octave)
            truthvalue= 0;
            % find if the centerlevel pixel is the local maxima or minima
            centerpixel  = array(i,j,centerlevel, octave);
            %find all the neighbouring 26 points.
            neighbours = array(i-1:i+1, j-1:j+1, centerlevel-1:centerlevel+1, octave);
            neighbours= neighbours(:);
            neighbours(14) = [];
            
            %find if greater than all values or smaller than all values
            maxvalue = max(neighbours);
            minvalue = min(neighbours);
            
            % store all the points in location so that they can be
            % used later
            if (centerpixel > maxvalue)
                truthvalue = 1;
            elseif (centerpixel < minvalue)
                truthvalue = 1;
            end
        end
        
        %% step - display all the potential points
        
        
        function displayKeypoints(self, locations)
            
            for octave = 1 : self.octavecount
                %which figure to display according to octave
                image1 = self.image;
                if octave == 2
                    image1 = self.image(1:2:end, 1:2:end);
                end
                
                for centerlevel = self.minlevel + 1 : self.maxlevel - 2
                    
                    %find the points at that octave and scale
                    [r,c] = find(locations(:,:,centerlevel, octave));                  % Find row,col coords.
                    
                    %display the found points
                    disp(['total points found at ' num2str(centerlevel) ' in octave' num2str(octave) ]);
                    disp(size(r))
                    
                    % overlay corners on original image
                    figure,
                    imagesc(image1),
                    axis image,
                    colormap(gray),
                    hold on;
                    %plot(c,r,'ys'),
                    title(['potential points detected', num2str(octave), ' - ', num2str(centerlevel) ]);
                    
                    [self.length1,~] = size(r);
                    %                     sigma = self.sigmavalue * (self.kvalue ^ (centerlevel-1));
                    radius = self.computeGaussianSigma(octave, centerlevel);
                    for pointindex = 1 : self.length1
                        y = r(pointindex);
                        x = c(pointindex);
                        ang=0:0.01:2*pi;
                        xp=radius*cos(ang);
                        yp=radius*sin(ang);
                        plot(x+xp,y+yp);
                    end
                    
                end
            end
        end
        
        %% step 2 - key point localization
        
        
        function outlierRejection(self)
            
            filterdx = [-1 0 1; -1 0 1; -1 0 1];   % Derivative masks
            filterdy = filterdx';
            self.finalpositions = zeros(self.rows, self.cols, self.length1, self.octavecount);
            
            for octave = 1 : self.octavecount
                for centerlevel = self.minlevel + 1 : self.maxlevel - 2
                    
                    % 1 . initial outlier rejection
                    %taylor series to find the proper location
                    
                    
                    if octave == 1
                        %find first derivative and second derivative at those points
                        Dx = conv2(self.dogvalues1(:,:,centerlevel-1,octave), filterdx, 'same');
                        Dxx = conv2(Dx, filterdx, 'same');
                        Dxy = conv2(Dx, filterdy, 'same');
                        Dy = conv2(self.dogvalues1(:,:,centerlevel-1,octave), filterdy, 'same');
                        Dyx = conv2(Dy, filterdx, 'same');
                        Dyy = conv2(Dy, filterdy, 'same');
                        
                        %                         %             dx1 =  conv2(dogvalues1(:,:,centerlevel+1,octave), filterdx, 'same');
                        %                         %             dxsigma = dx1 - dx;
                        
                        %                         %             dy1 = conv2(dogvalues1(:,:,centerlevel+1,octave), filterdy, 'same');
                        %                         %             dysigma = dy1 - dy;
                        %                         %             dsigma = dogvalues1(:,:,centerlevel+1,octave) - dogvalues1(:,:,centerlevel,octave);
                        %                         %             dsigmax = conv2(dsigma, filterdx, 'same');
                        %                         %             dsigmay = conv2(dsigma, filterdy, 'same');
                        %                         %             dsigma1 = dogvalues1(:,:,centerlevel+1,octave) - dogvalues1(:,:,centerlevel,octave);
                        %                         %             dsigmasigma = dsigma1 - dsigma;
                        %
                        %
                        %                         % find the taylor series values for the potential points
                        %                         taylorseries = zeros(self.rows, self.cols);
                        %
                        %                         [length11,~] = size(r);
                        %                         for pointindex = 1 : length11
                        %                             locx = r(pointindex,1);
                        %                             locy = c(pointindex,1);
                        %                             propercenterlevel = centerlevel;
                        %
                        %                             %find the extremum offset
                        %
                        %                             %find first derivative
                        %                             %                 dX = [dx(i,j)
                        %                             %                     dy(i,j)
                        %                             %                     dsigma(i,j)
                        %                             %                     ];
                        %                             dx = (self.dogvalues1(locx+1,locy,centerlevel,octave) - self.dogvalues1(locx-1,locy,centerlevel,octave))/2;
                        %                             dy =  (self.dogvalues1(locx,locy+1,centerlevel,octave) - self.dogvalues1(locx,locy-1,centerlevel,octave))/2;
                        %                             dsigma = (self.dogvalues1(locx,locy,centerlevel+1,octave) - self.dogvalues1(locx,locy,centerlevel-1,octave))/2;
                        %                             dX = [dx
                        %                                 dy
                        %                                 dsigma
                        %                                 ];
                        %
                        %                             %find hessian
                        %                             dx1 = (self.dogvalues1(locx+1,locy,centerlevel,octave) - self.dogvalues1(locx,locy,centerlevel,octave));
                        %                             dx2 = (self.dogvalues1(locx,locy,centerlevel,octave) - self.dogvalues1(locx-1,locy,centerlevel,octave));
                        %                             dx3 = (self.dogvalues1(locx+1,locy+1,centerlevel,octave) - self.dogvalues1(locx,locy+1,centerlevel,octave));
                        %                             dx4 = (self.dogvalues1(locx+1,locy,centerlevel+1,octave) - self.dogvalues1(locx,locy,centerlevel+1,octave));
                        %                             dxx = dx1-dx2;
                        %                             dxy = dx3 - dx1;
                        %                             dxsigma = dx4 - dx1;
                        %                             dy1 = (self.dogvalues1(locx,locy+1,centerlevel,octave) - self.dogvalues1(locx,locy,centerlevel,octave));
                        %                             dy2 = (self.dogvalues1(locx,locy,centerlevel,octave) - self.dogvalues1(locx,locy-1,centerlevel,octave));
                        %                             dy3 = (self.dogvalues1(locx+1,locy+1,centerlevel,octave) - self.dogvalues1(locx+1,locy,centerlevel,octave));
                        %                             dy4 = (self.dogvalues1(locx,locy+1,centerlevel+1,octave) - self.dogvalues1(locx,locy,centerlevel+1,octave));
                        %                             dyy = dy1-dy2;
                        %                             dyx = dy3 - dy1;
                        %                             dysigma = dy4 - dy2;
                        %                             dsigma1 = (self.dogvalues1(locx,locy,centerlevel+1,octave) - self.dogvalues1(locx,locy,centerlevel,octave));
                        %                             dsigma2 = (self.dogvalues1(locx,locy,centerlevel,octave) - self.dogvalues1(locx,locy,centerlevel-1,octave));
                        %                             dsigma3 = (self.dogvalues1(locx+1,locy,centerlevel+1,octave) - self.dogvalues1(locx+1,locy,centerlevel,octave));
                        %                             dsigma4 = (self.dogvalues1(locx,locy+1,centerlevel+1,octave) - self.dogvalues1(locx,locy+1,centerlevel,octave));
                        %                             dsigmasigma = dsigma1 - dsigma2;
                        %                             dsigmax = dsigma3 - dsigma1;
                        %                             dsigmay = dsigma4 - dsigma1;
                        %
                        %                             %                 hessian = [dxx(locx,j) dxy(locx,j) dxsigma(i,j)
                        %                             %                     dyx(i,j) dyy(i,j) dysigma(i,j)
                        %                             %                     dsigmax(i,j) dsigmay(i,j) dsigmasigma(i,j)];
                        %
                        %                             hessian = [dxx dxy dxsigma
                        %                                 dyx dyy dysigma
                        %                                 dsigmax dsigmay dsigmasigma];
                        %
                        %                             %find extremum offset
                        %                             extremum = - hessian'* dX;
                        %
                        %                             %add the extremum if the values are less than 1
                        %                             limit = 2;
                        %                             if abs(extremum(1)) <= limit && abs(extremum(2)) <= limit && abs(extremum(3)) <= 1
                        %                                 %                     disp('found extremum');
                        %                                 %                     disp(extremum);
                        %
                        %                                 % add the offset by rounding
                        %                                 locx = locx + round(extremum(1));
                        %                                 locy = locy + round(extremum(2));
                        %                                 propercenterlevel = propercenterlevel + round(extremum(3));
                        %
                        %                                 %                     disp('new extremum');
                        %                                 %                     disp(locx);
                        %                                 %                     disp(locy);
                        %                                 %                     disp(propercenterlevel);
                        %
                        %                                 %                     if abs(extremum(1)) > 0.5 || abs(extremum(2)) > 0.5 || abs(extremum(3)) > 0.5
                        %                                 %                         %disp(['hello: ', num2str(extremum(1)), ',' , num2str(extremum(2)) , ',' num2str(extremum(3))]);
                        %                                 %                     end
                        %                                 %
                        %                                 %                     if abs(extremum(1)) < 0.5 && abs(extremum(2)) < 0.5 && abs(extremum(3)) < 0.5
                        %                                 %
                        %                                 %                         disp('hello');
                        %                                 %                     end
                        %
                        %                                 taylorseries(locx,locy) = self.dogvalues1(locx,locy,propercenterlevel,octave)+  0.5 * ( extremum(1)*dX(1) + extremum(2)*dX(2) + extremum(3)* dX(3)) ;
                        %
                        %                                 %                     disp(['tailor series ' num2str(taylorseries(locx,locy))] );
                        %                                 finalpositions(locx,locy,propercenterlevel,octave) = 1;
                        %                                 %remove all the points with less contrast
                        %                                 if abs( taylorseries(locx,locy)) > 0.03;
                        %                                     %                         finalpositions(locx,locy,propercenterlevel,octave) = 1;
                        %                                 end
                        %                             end
                        %                         end
                        
                    end
                end
                
            end
        end
        
        
        %% perform harris corner and reject all edges or patches
        
        function harrisCornerRejection(self)
            
            self.finalpositions = zeros(self.rows, self.cols, self.length1, self.octavecount);
            
            for octave = 1 : self.octavecount
                %which dog value to use for derivative computation
                if octave == 1
                    array = self.dogvalues1;
                else
                    array = self.dogvalues2;
                end
                
                for centerlevel = self.minlevel + 1 : self.maxlevel - 2
                    
                    disp(['at level ' num2str(centerlevel) ' of octave ' num2str(octave)]);
                    
                    %find the points at that octave and scale, Find row,col coords.
                    [r,c] = find(self.locations(:,:,centerlevel, octave));
                    
                    [Dxx, Dxy, Dyy] = self.findDerivatives(array, centerlevel);
                    %                     disp(Dxx);
                    % 2. further outlier rejection
                    [lengtha,~] = size(r);
                    disp(lengtha);
                    for k = 1 : lengtha
                        i = r(k,1);
                        j = c(k,1);
                        %use harris corner detector
                        hessian = [Dxx(i,j) Dxy(i,j)
                            Dxy(i,j) Dyy(i,j)];
                        eigenvalues = eig(hessian);
                        if eigenvalues(2) ~= 0
                            ratio = eigenvalues(1)/eigenvalues(2);
                            if ratio < 0
                                ratio = -ratio;
                            end
                            %remove all the edges
                            if abs(ratio) <= 10
                                self.finalpositions(i,j, centerlevel ,octave)=1;
                            end
                        end
                    end
                    
                end
            end
            
        end
        
        function [Dxx, Dxy, Dyy] =  findDerivatives(self, array, level)
            % Derivative masks
            filterdx = [-1 0 1; -1 0 1; -1 0 1];
            filterdy = filterdx';
            
            %find first derivative and second derivative at those points
            Dx = conv2(array(:,:,level), filterdx, 'same');
            Dxx = conv2(Dx, filterdx, 'same');
            Dxy = conv2(Dx, filterdy, 'same');
            Dy = conv2(array(:,:,level), filterdy, 'same');
            Dyy = conv2(Dy, filterdy, 'same');
        end
        
        
        %% compute the dominant direction
        function [magnitude] = calculateMagnitude(self, array, i, j, centerlevel)
            dx = array(i+1,j,centerlevel)-array(i-1,j,centerlevel);
            dy = array(i,j+1,centerlevel)-array(i,j-1,centerlevel);
            magnitude = sqrt(dx*dx + dy*dy);
        end
        
        function [orientation] = calculateOrientation(self, array, i, j, centerlevel)
            dx = array(i+1,j,centerlevel)-array(i-1,j,centerlevel);
            dy = array(i,j+1,centerlevel)-array(i,j-1,centerlevel);
            if dx == 0
                if dy == 0
                    orientation = 1;
                    return
                else
                    dx = 1;
                    dy = 100000;
                end
            end
            orientation = atand(dy/dx);
            
%             disp(['angle = ' num2str(orientation) ' dx = ' num2str(dx) ' dy = ' num2str(dy)]);
            %correct the orientation from -90<->90 to 0<->360 range
            if orientation > 0 && orientation < 90
                if array(i-1,j,centerlevel) > array(i+1,j,centerlevel)
                    orientation = orientation + 180;
                end
            elseif orientation <= 0 && orientation > -90
                if array(i-1,j,centerlevel) > array(i+1,j,centerlevel)
                    orientation = orientation + 360;
                else
                    orientation = orientation + 180;
                end
            end
        end
        
        function [magnitudes] = smoothMagnitudeGaussian(self, magnitudes, centerlevel, octave)
            %smooth the magnitude values with gaussian of sigma
            %with 1.5 sigma
            sigma = 1.5* self.computeGaussianSigma(octave, centerlevel-1);
            width = self.computeGaussianRadius(sigma);
            gfilter = fspecial('gaussian', width, sigma);
            magnitudes = conv2(magnitudes, gfilter, 'same');
        end
        
        function [histogram] = computeHistogram(self, binsize, magnitudes, orientations)
            %find histogram of orientation in neighbourhood
            histogram = zeros(1,binsize);
            %                         bins = 1:binsize;
            [lengtha,lengthb] = size(magnitudes);
            angle = 360/binsize;
            for a = 1 : lengtha
                for b = 1 : lengthb
%                     disp(['angle2 = ' num2str(orientations(a,b))]);
                    bin = ceil(orientations(a,b)/angle);
                    histogram(1,bin ) =  histogram(1,bin) + magnitudes(a,b);
                end
            end
            
            % Example Visualise:
            %figure;
            %bar(bins, histogram)
        end
        
        %% step - display all the potential points and their dominant
        %% directions
        
        
        function displayKeypointsDirections(self, locations, directions)
            
            index = 1;
            
            for octave = 1 : self.octavecount
                %which figure to display according to octave
                image1 = self.image;
                if octave == 2
                    image1 = self.image(1:2:end, 1:2:end);
                end
                
                for centerlevel = self.minlevel + 1 : self.maxlevel - 2
                    
                    %find the points at that octave and scale
                    [r,c] = find(locations(:,:,centerlevel, octave));                  % Find row,col coords.
                    
                    %display the found points
                    disp(['total points found at ' num2str(centerlevel) ' in octave' num2str(octave) ]);
                    disp(size(r))
                    
                    % overlay corners on original image
                    figure,
                    imagesc(image1),
                    axis image,
                    colormap(gray),
                    hold on;
                    %plot(c,r,'ys'),
                    title(['potential points detected', num2str(octave), ' - ', num2str(centerlevel) ]);
                    
                    [self.length1,~] = size(r);
                    %                     sigma = self.sigmavalue * (self.kvalue ^ (centerlevel-1));
                    radius = self.computeGaussianSigma(octave, centerlevel);
                    
                    for pointindex = 1 : self.length1
                        y = r(pointindex);
                        x = c(pointindex);
                        ang=0:0.01:2*pi;
                        xp=radius*cos(ang);
                        yp=radius*sin(ang);
                        plot(x+xp,y+yp);
                        
                        % display the dominant direction
                        dominantdirection = directions(index);
                        x1 = x;%-radius*cos(dominantdirection)*0.5;
                        x2 = x+radius*cos(dominantdirection);
                        y1 = y;%-radius*sin(dominantdirection)*0.5;
                        y2 = y+radius*sin(dominantdirection);
                        
                        x = linspace(x1,x2,100);
                        y = linspace(y1,y2,100);
                        plot(x,y);
                        
                        index=index+1;
                    end
                    
                    
                end
            end
        end
        
        
        function [dominantDirections] = calculateDominantOrientation(self)
            
            %calculate orientation
            %             featurelocations = [];
            dominantDirections = [];
            
            for octave = 1 : self.octavecount
                
                if octave == 1
                    array = self.gaussianvalues1;
                    rowvalue = self.rows;
                    colvalue=self.cols;
                else
                    array = self.gaussianvalues2;
                    rowvalue = ceil(self.rows/2);
                    colvalue = ceil(self.cols/2);
                end
                
                for centerlevel = self.minlevel + 1 : self.maxlevel - 2
                    
                    [r,c] = find(self.finalpositions(:,:,centerlevel, octave));                  % Find row,col coords.
                    [lengtha,~] = size(r);
                    %do edge outliers
                    for k = 1 : lengtha
                        locx = r(k);
                        locy = c(k);
                        %   featurelocations = [featurelocations ; [locx, locy, centerlevel, octave] ];
                        
                        magnitudes = [];
                        orientations = [];
                        
                        %find orientation and magnitude in the 4x4 neighbourhood
                        m1 = locx -1;
                        m2 = locx + 2;
                        n1 = locy - 1;
                        n2 = locy + 2;
                        % correct the boundary
                        if m1 < 2
                            delta = 2 - m1;
                            m2 = m2 + delta;
                            m1 = 2;
                        elseif m2 > rowvalue - 2
                            delta = m2 - rowvalue + 1;
                            m1 = m1 - delta;
                            m2 = rowvalue - 1;
                        end
                        if n1 < 2
                            delta = 2 - n1;
                            n2 = n2 + delta;
                            n1 = 2;
                        elseif n2 > colvalue - 2
                            delta = n2 - colvalue + 1;
                            n1 = n1 - delta;
                            n2 = colvalue - 1;
                        end
                        for i = m1 : m2
                            
                            magnitudes1 = [];
                            orientations1 = [];
                            
                            for j = n1 : n2
                                
                                %calculate magnitude and orientation
                                magnitude = self.calculateMagnitude(array, i, j, centerlevel);
                                magnitudes1 =[magnitudes1 magnitude];
                                orientation = self.calculateOrientation(array, i, j, centerlevel);
                                orientations1 = [orientations1 orientation];
                                
                            end
                            
                            magnitudes =[magnitudes ;magnitudes1];
                            orientations = [orientations; orientations1];
                        end
                        
                        %smooth magnitude with gaussian
                        magnitudes = self.smoothMagnitudeGaussian(magnitudes, centerlevel, octave);
                        
                        histogram = self.computeHistogram(36, magnitudes, orientations);
                        
                        %find the dominant directions in the histogram
                        combined = [histogram; 1:36]';
                        sortedvalues = sortrows(combined, 1);
                        %                         maxvalue = [sortedvalues(36,1)];
                        domiantdirections = [sortedvalues(36,2)];
                        %                         for m = 35 : -1 : 1
                        %                             if sortedvalues(m) > 0.8 * maxvalue(1)
                        %                                 % maxvalue = [maxvalue sortedvalues(m,1)];
                        %                                 domiantdirections = [domiantdirections sortedvalues(m,2)];
                        %                             end
                        %                         end
                        
                        dominantDirections = [dominantDirections domiantdirections(1)];
                        
                    end
                end
            end
        end
        
        
        %% compute sift feature descriptor
        function [featuredescriptions, featurelocations] =calculateSIFTDescriptor(self, dominantDirections)
            
            first = 1;
            featuredescriptions = [];
            featurelocations = [];
            index= 1;
            
            for octave = 1 : self.octavecount
                for centerlevel = self.minlevel + 1 : self.maxlevel - 2
                    
                    if octave == 1
                        array = self.gaussianvalues1;
                        rowvalue = self.rows;
                        colvalue=self.cols;
                    else
                        array = self.gaussianvalues2;
                        rowvalue = ceil(self.rows/2);
                        colvalue = ceil(self.cols/2);
                    end
                    
                    [r,c] = find(self.finalpositions(:,:,centerlevel, octave));                  % Find row,col coords.
                    [lengtha,~] = size(r);
                    %do edge outliers
                    for k = 1 : lengtha
                        locx = r(k);
                        locy = c(k);
                        
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
                        elseif m2 > rowvalue - 2
                            delta = m2 - rowvalue + 1;
                            m1 = m1 - delta;
                            m2 = rowvalue - 1;
                        end
                        if n1 < 2
                            delta = 2 - n1;
                            n2 = n2 + delta;
                            n1 = 2;
                        elseif n2 > colvalue - 2
                            delta = n2 - colvalue + 1;
                            n1 = n1 - delta;
                            n2 = colvalue - 1;
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
                                    magnitudes1 = [];
                                    orientations1 = [];
                                    
                                    for b = 1:4
                                        
                                        i = mvalue + a - 1;
                                        j = nvalue + b - 1;
                                        
                                        %calculate magnitude and orientation
                                        magnitude = self.calculateMagnitude(array, i, j, centerlevel);
                                        magnitudes1 =[magnitudes1 magnitude];
                                        orientation = self.calculateOrientation(array, i, j, centerlevel) - dominantDirections(index) ;
                                        if orientation <= 0
                                            orientation= orientation + 360;
                                        end
                                        orientations1 = [orientations1 orientation];
                                        
                                    end
                                    
                                    magnitudes =[magnitudes ;magnitudes1];
                                    orientations = [orientations; orientations1];
                                end
                                
                                %smooth magnitude with gaussian
                                magnitudes = self.smoothMagnitudeGaussian(magnitudes, centerlevel, octave);
                                
                                histogram = self.computeHistogram(8, magnitudes, orientations);
                                
                                %append to feature descriptor
                                featuredescriptor = [featuredescriptor histogram];
                                
                                
                                
                            end
                        end
                        
%                         if first
%                             
%                             figure;
%                             subimage = array(m1:m2, n1:n2, centerlevel);
%                             imshow(uint8(subimage));
%                             
%                             figure;
%                             [~,lengthc] = size(featuredescriptor);
%                             bins = 1 : lengthc;
%                             bar(bins, featuredescriptor)
%                             
%                             first = 0;
%                         end
                        
                        %disp(featuredescriptor);
                        featuredescriptions= [featuredescriptions ; [featuredescriptor]];
                        
                        featurelocations = [featurelocations ; [locx, locy, octave]];
                        
                        index = index+1;
                    end
                    
                end
                
            end
            
        end
        
        %% post processing
        function [featuredescriptions] = postprocessDescriptors(self, featuredescriptions)
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
        
        %% end of all the methods
    end
    
    
    
    
end


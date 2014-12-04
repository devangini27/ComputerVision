image = imread('railway.jpg');
figure(1);
imshow(image);
title('image');


% step 1 - smooth image with gaussian
gaussian= fspecial('gaussian', [7 7], 1);
smooth = conv2(double(image), gaussian, 'same');
figure(2);
imshow(uint8(smooth));
title('smoothed');

% step 2 - computer derivative of image
filterdx = [1 0 -1
    2 0 -2
    1 0 -1];
filterdy = [1 2 1
    0 0 0
    -1 -2 -1];
dx = conv2(double(smooth), filterdx, 'same');
dy = conv2(double(smooth), filterdy, 'same');
figure(3);
imshow(uint8(dx));
title('dx');
figure(4);
imshow(uint8(dy));
title('dy');


% step 3 - find magnitude and orientation
magnitude = dx.*dx + dy.*dy;
magnitude = sqrt(double(magnitude));
orientation = atand(double(dx./dy));
figure(5);
imshow(uint8(magnitude));
title('magnitude');


%step 4 - non maximum suppression
[rows,cols] = size(magnitude);
for i = 2 : rows - 1
    for  j = 2 : cols - 1
        disp(['i = ' num2str(i) ' j = ' num2str(j) ]);
        angle = orientation(i,j);
        neighbour1 = 0;
        neighbour2 = 0;
        if  angle >=-22.5 && angle <22.5
            neighbour1 =magnitude(i, j+1);
            neighbour2 = magnitude(i, j - 1);
        elseif angle >=22.5 && angle < 67.5
            neighbour1 = magnitude(i+1, j+1);
            neighbour2 = magnitude(i-1, j - 1);
        elseif angle >=67.5 && angle <=90
            neighbour1 = magnitude(i+1, j);
            neighbour2 = magnitude(i-1, j);
        elseif angle <= -22.5 && angle > -67.5
            neighbour1 = magnitude(i-1, j+1);
            neighbour2 = magnitude(i+1, j - 1);
        elseif angle <=- 67.5 && angle >= -90
            neighbour1 = magnitude(i+1, j);
            neighbour2 = magnitude(i-1, j);
        else
            disp(['unknown angle '  num2str(angle)] );
        end
        
        if magnitude(i,j) > neighbour1 && magnitude(i,j) > neighbour2
        else
            magnitude(i,j) = 0;
        end
    end
end
figure(6);
imshow(uint8(magnitude));
title('magnitude');

disp('starting hystersis thresholding');

% step 5 - hystersis thresholding
threshold1 = 200;
threshold2 = 220;
neighbours = [];
edges = zeros(rows,cols);
for i = 1 : rows
    for j = 1 : cols
        disp(['i = ' num2str(i) ' j = ' num2str(j) ]);
        if magnitude(i,j) >= threshold2
            edges(i,j) = 255;
            % look at the four neighbours
            if i > 1 % top neighbour
                neighbour = [ [i-1; j] ];
                neighbours = [neighbours neighbour];
            elseif i < rows % bottom neighbour
                neighbour = [ [i+1; j] ];
                neighbours = [neighbours neighbour];
            elseif j > 1 % leftneighbour
                neighbour = [ [i; j-1] ];
                neighbours = [neighbours neighbour];
            elseif j < cols
                neighbour = [ [i; j+1] ];
                neighbours = [neighbours neighbour];
            end
            
            %consume all the nighbours
            while ~isempty(neighbours)
                
                neighbour = neighbours(:,1);
                if length(neighbours) > 1
                       neighbours = neighbours(3 : end);
                end
                if magnitude(neighbour(1), neighbour(2)) > threshold1
                    edges( neighbour(1), neighbour(2) ) = 255;
                    
                    disp('recursive neighbours');
                    k = neighbour(1);
                    l = neighbour(2);
                    % add all its neighbours
                    if k > 1 % top neighbour
                        neighbour = [ [k-1; l] ];
                        neighbours = [neighbours neighbour];
                    elseif k < rows % bottom neighbour
                        neighbour = [ [k+1; l] ];
                        neighbours = [neighbours neighbour];
                    elseif l > 1 % leftneighbour
                        neighbour = [ [k; l-1] ];
                        neighbours = [neighbours neighbour];
                    elseif l < cols
                        neighbour = [ [k; l+1] ];
                        neighbours = [neighbours neighbour];
                    end
                    
                end
                
            end;
            
            
            
        end
    end
end
figure(7);
imshow(edges);
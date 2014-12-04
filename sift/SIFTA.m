matchingfeatures = [];
count = 0;

for i = 1:features1size
    
    minindex1 = -1;
    minvalue1 = 100000;
    
    minindex2 = -1;
    minvalue2 = 100000;
    
    for j = 1:features2size
        %find distance between the two features
        
        distance = 0;
        for k = 1:128
            diff = features1(i,k) - features2(j,k);
            distance = distance + (diff)^2;
        end
        distance = sqrt(distance);
        
        
        %initialize min values
        if minindex2 == -1
            minindex2 = j;
            minvalue2 = distance;
        elseif minindex1 == -1
            if distance < minvalue2
                minindex1 = j;
                minvalue1 = distance;
            else
                minindex1 = minindex2;
                minvalue1 = minvalue2;
                minindex2 = j;
                minvalue2 = distance;
            end
        else
            %now compare properly
            if distance < minvalue1
                %shift 1 to 2
                minindex2 = minindex1;
                minvalue2 = minvalue1;
                minindex1 = j;
                minvalue1 = distance;
            elseif distance >= minvalue1 && distance < minvalue2
                %replace 2
                minindex2 = j;
                minvalue2 = distance;
            end
        end
        
    end
    
%     disp(['distances = ' num2str(minvalue1) ' , ' num2str(minvalue2)]);
    
    ratio = minvalue1/minvalue2;
%     disp(ratio);
    if ratio <= 0.5
        matchingfeatures = [matchingfeatures ; [i minindex1]];
        count = count +1;
    end
end

disp(['number of matching features ' num2str(count)]);


image1 = imread(imagename1);
[~,cols,channels] = size(image1);
if channels == 3
    image1 = rgb2gray(image1);
end
image2 = imread(imagename2);
[~,~,channels] = size(image2);
if channels == 3
    image2 = rgb2gray(image2);
end
image = [image1(:,:) image2(:,:)];

figure;
imshow(image);
hold on;

for k = 1 : count
    
    i = matchingfeatures(k,1);
    j = matchingfeatures(k,2);
    
    x1 = locations1(i, 2);
    y1 = locations1(i, 1);
    if locations1(i,3) == 2
        x1 = x1 * 3;
        y1 = y1 * 2;
    end
    
    x2 = locations2(j, 2);
    y2 = locations2(j, 1);
    if locations2(j,3) == 2
        x2 = x2 * 3;
        y2 = y2 * 2;
    end
    
    disp(['point1 = (' num2str(x1) ',' num2str(y1) ') point2 = ' num2str(x2) ',' num2str(y2) ')']);
    
    x2 = x2 + cols;
    
    x = linspace(x1,x2,100);
    y = linspace(y1,y2,100);
    plot(x,y);
end
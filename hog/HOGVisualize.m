function HOGVisualize(image, hogfeatures, logicalcellsindex)
% Visualize HOG features

% The code below takes a HOG descriptor and computes the gradient magnitude per cell
% by looping through all blocks and adding the gradient strengths for each of the cells per block to each individual cell
% and finally averages the gradient strenghts.
% The gradient strengths are then visualized by line lengths in the direction of the corresponding gradient bin.

figure;

image(8:8:end,:,:) = 0;       %# Change every tenth row to black
image(:,8:8:end,:) = 0;       %# Change every tenth column to black

[rows, cols]=size(image);
% blankimage= zeros(rows, cols);
imshow(image);
hold on;


%total number of cells (logical)
% 64x128 / 8x8 = 8 x 16 = 128
cellcount = 128;
binsize = 9;
for cellindex = 1 : cellcount
    
    %find the corresponding logical indexes
    indexes = find(logicalcellsindex==cellindex);
    %disp(['found cellindex = ' num2str(cellindex) ' at ' mat2str(indexes)]);
    
    [~,count] = size(indexes);
    avghistogram = zeros(1,binsize);
    for blocksindex = 1: count
        %find the corresponding histograms
        block = indexes(blocksindex);
        hogindex = (block-1)*binsize + 1;
        histogram = hogfeatures(1,hogindex: hogindex + binsize - 1);
        avghistogram = avghistogram + histogram;
        %         disp(histogram);
    end
    
    avghistogram = avghistogram/count;
    %divide by 255;
    avghistogram = avghistogram/255;
    
    maxvalue = avghistogram;
    domiantdirections = 1:binsize;
%     domiantdirections = domiantdirections';
    
%     disp(avghistogram);
    
%     dominantcond = 1;
%     
%     if dominantcond
%         %find dominant directions
%         combined = [avghistogram; 1:binsize]';
%         sortedvalues = sortrows(combined, 1);
%         maxvalue = [sortedvalues(binsize,1)];
%         domiantdirections = [sortedvalues(binsize,2)];
%         for m = binsize - 1 : -1 : 1
%             if sortedvalues(m,1) > 0.8 * maxvalue(1,1)
%                 % maxvalue = [maxvalue sortedvalues(m,1)];
%                 domiantdirections = [domiantdirections sortedvalues(m,2)];
%                 maxvalue = [maxvalue sortedvalues(m,1)];
%             end
%         end
%         domiantdirections = 20 * (1:binsize);
%         maxvalue = avghistogram;
%     else
%         
%         
%         
%     end
    
    %find middle of the cell in the image
    %get x, y in terms of cell first
    x = floor((cellindex - 1)/8) + 1;
    y = (cellindex - (x-1) * 8);
    %mutliply by 8
    x = (x-1)*8 + 1;
    y = (y-1)*8 + 1;
    %add 3 to get the center
    x = x + 3;
    y = y + 3;
    
    %draw the lines for the dominant directions
%     disp(domiantdirections);
%     disp(maxvalue);
    
    [~,length] = size(domiantdirections);
    for i =1 : length
%         disp(domiantdirections(1,i));
%         disp(maxvalue(1,i));
        
        angle = domiantdirections(1,i);
        length = maxvalue(1,i);
        viz_factor = 1;%5;%3;
        xcomponent = round(length/2 * cosd(angle * 20 - 10)* viz_factor);
        ycomponent = round(length/2 * sind(angle * 20 - 10)* viz_factor);
        
        %disp(['cell index ' num2str(cellindex) ' to pixel index = ( ' num2str(x) ' , ' num2str(y) ') and line created, xcomp =  ' num2str(xcomponent) ' and y comp = ' num2str(ycomponent) ]);
        x1 = x-xcomponent;
        y1 = y-ycomponent;
        x2 = x+xcomponent;
        y2 = y+ycomponent;
        %disp(['line drawn from ( ' num2str(x1) ' , ' num2str(y1) ' ) to ( ' num2str(x2) ' , ' num2str(y2) ' ) ']);
        
        % plot([cols1,cols2],[rows1,rows2],'Color','r','LineWidth',2)
        plot([y1,y2],[x1,x2],'Color','r','LineWidth',1)
        
    end
    
end

end

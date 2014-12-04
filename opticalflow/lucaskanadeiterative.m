function [u0,v0] = lucaskanadeiterative(image1, image2, u3, v3)

[rows, cols] = size(image1);

windowx= 2;
windowy =2;


[fx,fy,ft1] = calculatederivatives(image1, image2);


%calculate the value of ft for a block
%ft = fx*u  + fy*v
ft = zeros(rows, cols);
for i = 1 : rows
    for j = 1 : cols
        disp(['optical flow for (' num2str(i) ',' num2str(j) ') = ' num2str(u3(i,j)) ',' num2str(v3(i,j))]);
        
        % perform interpolation
        x = i+u3(i,j);
        y = j+v3(i,j);
        floorx = floor(x);
        floory = floor(y);
        ceilx = ceil(x);
        ceily = ceil(y);
        dfloorx = x - floorx;
        dfloory = y - floory;
        dceilx = 1 - dfloorx;
        dceily = 1 - dfloory;
        
        %change vales if at boundary
        if floorx == 0
            floorx = 1;
            %keep dfloorx as it is
        end
        if floory == 0
            floory = 1;
            %keep dfloory as it is
        end
        %perform flow interpolation
        image2value= dceilx*dceily * image2(floorx, floory) + dceilx*dfloory * image2(floorx, ceily) +  dfloorx*dceily * image2(ceilx, floory) + dfloorx*dfloory * image2(ceilx, ceily);
        
        
        
        ft(i,j) = image2value - image1(i,j);
        % ft(i,j) = image2(i+u3(i,j) ,j+v3(i,j)) - image1(i,j);
    end
end

%do lucas kanade
g = zeros(rows, cols,3);


%for all pixels
for i = windowx+1 : rows - windowx
    for j = windowy+1 : cols - windowy
        
        sumfx2 = 0;
        sumfy2 = 0;
        sumfxfy = 0;
        
        %construct a 3x3 window
        for a = -windowx : windowx
            for b = -windowy : windowy
                
                %calculate all the sums
                sumfx2 = sumfx2 + fx(i+a, j+b)^2;
                sumfy2 = sumfy2 + fy(i+a, j+b)^2;
                sumfxfy = sumfxfy + fx(i+a, j+b) * fy(i+a, j+b);
            end
        end
        
        %compose spatial gradient matrix
        g(i,j,1) = sumfx2;
        g(i,j,2) = sumfxfy;
        g(i,j,3) =  sumfy2;
    end
end

u0 = zeros(rows, cols);
v0 = zeros(rows, cols);
difference = zeros(rows, cols);

for i = 1 + windowx : rows - windowx
    for j = 1 + windowy : cols - windowy
        %calculate image difference
        x = i+u3(i,j) + u0(i,j);
        y = j+v3(i,j) + v0(i,j);
        
        disp(['i = ' num2str(i) ' , u3 = ' num2str(u3(i,j)) ' ,u0 = ' num2str(u0(i,j)) ' , x = ' num2str(x)]);
        disp(['j = ' num2str(i) ' , v3 = ' num2str(u3(i,j)) ' ,v0 = ' num2str(u0(i,j)) ' , y = ' num2str(x)]);
        
        floorx = floor(x);
        floory = floor(y);
        ceilx = ceil(x);
        ceily = ceil(y);
        dfloorx = x - floorx;
        dfloory = y - floory;
        dceilx = 1 - dfloorx;
        dceily = 1 - dfloory;
        
        %change vales if at boundary
        if floorx == 0
            floorx = 1;
            %keep dfloorx as it is
        end
        if floory == 0
            floory = 1;
            %keep dfloory as it is
        end
        %perform flow interpolation
        image2value= dceilx*dceily * image2(floorx, floory) + dceilx*dfloory * image2(floorx, ceily) +  dfloorx*dceily * image2(ceilx, floory) + dfloorx*dfloory * image2(ceilx, ceily);
        
        difference(i,j) =image1(i,j) - image2value;
    end
end

accuracythreshold = 0.5;

for i = 1 + windowx : rows - windowx
    for j = 1 + windowy : cols - windowy
        accuracy = 1;
        
        while accuracy < accuracythreshold
            sumdiffx = 0;
            sumdiffy = 0;
            for a = -windowx : windowx
                for b = -windowy : windowy
                    
                    %calculate all the sums
                    sumdiffx = sumdiffx + difference(i+a, j+b) * fx(i,j);
                    sumdiffy = sumdiffy + difference(i+a, j+b) * fy(i,j);
                end
            end
            mismatcherror = [sumdiffx
                sumdiffy];
            gmatrix = [g(i,j,1) g(i,j,2)
                g(i,j,2) g(i,j,3)];
            
            gmatrix
            mismatcherror
            
            if det(gmatrix) == 0
                error = [0 0];
            else
                error = inv(gmatrix ) * mismatcherror;
            end
            
            error
            
            u0(i,j) = u0(i,j) + error(1);
            v0(i,j) = v0(i,j) + error(2);
            
            accuracy = sqrt(error(1) * error(1) + error(2)*error(2));
            
        end
        
    end
end



figure;
% imshow(image1);
% hold on;
quiver(u0,v0);

end
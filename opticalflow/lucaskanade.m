function [u,v] = lucaskanade(image1, image2)

[rows, cols] = size(image1);
[fx,fy,ft] = calculatederivatives(image1, image2);


%do lucas kanade
u = zeros(rows, cols);
v = zeros(rows, cols);

%for all pixels
for i = 2 : rows - 1
    for j = 2 : cols -1
        
        sumfx2 = 0;
        sumfy2 = 0;
        sumfxfy = 0;
        sumfxft = 0;
        sumfyft = 0;
        
        %construct a 3x3 window
        for a = -1 : 1
            for b = -1 : 1
                
                %calculate all the sums
                sumfx2 = sumfx2 + fx(i+a, j+b)^2;
                sumfy2 = sumfy2 + fy(i+a, j+b)^2;
                sumfxfy = sumfxfy + fx(i+a, j+b) * fy(i+a, j+b);
                sumfxft = sumfxft + fx(i+a, j+b) * ft(i+a, j+b);
                sumfyft = sumfyft + fy(i+a, j+b) * ft(i+a, j+b);
                
            end
        end
        
        %calculate the constant demoninator
        denominator = sumfx2 * sumfy2 - sumfxfy ^ 2;

        if denominator == 0
            u(i,j) = 0;
            v(i,j) = 0;
        else
        %calculate the values of u and v
        u(i,j) = - sumfy2 * sumfxft + sumfxfy * sumfyft;
        u(i,j) = u(i,j) / denominator;
        v(i,j) = sumfxft * sumfxfy - sumfx2 * sumfyft;
        v(i,j) = v(i,j) / denominator;
        end
        
    end
end

figure;
% imshow(image1);
% hold on;
quiver(u,v);

end
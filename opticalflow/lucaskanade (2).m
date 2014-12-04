function [u,v] = lucaskanade(image1, image2)

[rows, cols] = size(image1);

%calculate fx, fy, ft
filterdx = [-0.5 0.5; -0.5 0.5];
filterdy  = filterdx';
filterdt = [-1 -1; -1 -1];

fx1 = conv2(double(image1) , filterdx, 'same');
fx2 = conv2(double(image2) , filterdx, 'same');
fx = fx1 + fx2;

figure;
imshow(fx);

fy1 = conv2(double(image1) , filterdy, 'same');
fy2 = conv2(double(image2) , filterdy, 'same');
fy = fy1 + fy2;

figure;
imshow(fy);

ft1 = conv2(double(image1) , filterdt, 'same');
ft2 = conv2(double(image2) , -filterdt, 'same');
ft = ft1 + ft2;

figure;
imshow(ft);

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
        
        %calculate the values of u and v
        u(i,j) = - sumfy2 * sumfxft + sumfxfy * sumfyft;
        u(i,j) = u(i,j) / denominator;
        v(i,j) = sumfxft * sumfxfy - sumfx2 * sumfyft;
        v(i,j) = v(i,j) / denominator;
        
    end
end

figure;
% imshow(image1);
% hold on;
quiver(u,v);

end
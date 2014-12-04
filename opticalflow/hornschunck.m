function [u,v] = hornschunck(image1, image2)

[rows, cols] = size(image1);



%calculate fx, fy, ft
filterdx = [-0.5 0.5; -0.5 0.5];
filterdy  = filterdx';
filterdt = [-1 -1; -1 -1];

fx1 = conv2(double(image1) , filterdx, 'same');
fx2 = conv2(double(image2) , filterdx, 'same');
fx = fx1 + fx2;

% figure;
% imshow(fx);

fy1 = conv2(double(image1) , filterdy, 'same');
fy2 = conv2(double(image2) , filterdy, 'same');
fy = fy1 + fy2;

% figure;
% imshow(fy);

ft1 = conv2(double(image1) , filterdt, 'same');
ft2 = conv2(double(image2) , -filterdt, 'same');
ft = ft1 + ft2;

% figure;
% imshow(ft);

%do iterative horn schunck
k=0;
lambda = 0.25;
u = zeros(rows, cols);
v = zeros(rows, cols);

threshold = 0.25;


for k = 0 : 100

    udiff = u;
    vdiff = v;

    %calculate average values for all the neighbours
    filteravg = [0 0.25 0; 0.25 0 0.25; 0 0.25 0];
    uavg = conv2(double(u), filteravg, 'same');
    vavg = conv2(double(v), filteravg, 'same');
    
    %calculate the parallel and normal flows
    p = fx.*uavg + fy.*vavg + ft;
    d = lambda + fx.^2 + fy.^2;
    ratio = p./d;
    
    %calculate the corrected values of u and v
    u = uavg - fx.*ratio;
    v = vavg - fy.*ratio;
    
    udiff = udiff - u;
    vdiff = vdiff - v;
    
    if (sign(udiff) < threshold) & (sign(vdiff) < threshold)
        disp(['k = ' num2str(k)])
        break
    end
    
end

end
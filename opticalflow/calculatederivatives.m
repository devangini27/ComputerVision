function [fx,fy,ft] = calculatederivatives(image1, image2)

%calculate fx, fy, ft
filterdx = 0.25*[-1 1; -1 1];
filterdy  = filterdx';
filterdt = 0.25*[-1 -1; -1 -1];

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


end
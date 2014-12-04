function [u,v] = hornschunck1(image1, image2)

% http://dspace.mit.edu/bitstream/handle/1721.1/6337/AIM-572.pdf?sequence=2
% http://www.journal.au.edu/au_techno/2011/july2011/journal151_article02.pdf

[rows, cols] = size(image1);



% smooth the images

% sigma=1;
% kSize=2*(sigma*3);
% x=-(kSize/2):(1+1/kSize):(kSize/2);
% G=(1/(sqrt(2*pi)*sigma)) * exp (-(x.^2)/(2*sigma^2));
% image1=conv2(image1,G,'same');
% image1=conv2(image1,G','same');
% image2=conv2(image2,G,'same');
% image2=conv2(image2,G','same');

[fx,fy,ft] = calculatederivatives(image1, image2);

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
    
    % Averaging kernel
    kernel_avg=[1/12 1/6 1/12;1/6 0 1/6;1/12 1/6 1/12];

    %filteravg = [0 0.25 0; 0.25 0 0.25; 0 0.25 0];
    uavg = conv2(double(u), kernel_avg, 'same');
    vavg = conv2(double(v), kernel_avg, 'same');
    
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
% sigma=1
% theta=0
% lambda=1
% psi=1
% gamma=1
%  
% sigma_x = sigma;
% sigma_y = sigma/gamma;
%  
% % Bounding box
% nstds = 3;
% xmax = max(abs(nstds*sigma_x*cos(theta)),abs(nstds*sigma_y*sin(theta)));
% xmax = ceil(max(1,xmax));
% ymax = max(abs(nstds*sigma_x*sin(theta)),abs(nstds*sigma_y*cos(theta)));
% ymax = ceil(max(1,ymax));
% xmin = -xmax; ymin = -ymax;
% [x,y] = meshgrid(xmin:xmax,ymin:ymax);
%  
% % Rotation 
% x_theta=x*cos(theta)+y*sin(theta);
% y_theta=-x*sin(theta)+y*cos(theta);
%  
% gb= exp(-.5*(x_theta.^2/sigma_x^2+y_theta.^2/sigma_y^2)).*cos(2*pi/lambda*x_theta+psi);

img = imread('cameraman.tif');
gaborArray = gaborFilterBank(5,8,39,39);  % Generates the Gabor filter bank
featureVector = gaborFeatures(img,gaborArray,4,4);   % Extracts Gabor feature vector, 'featureVector', from the image, 'img'.
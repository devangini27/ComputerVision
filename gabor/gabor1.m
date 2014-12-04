function [gb] = gabor1(sigma, gamma, theta, lambda, psi)

% sigma = 3;
% gamma = 2;
% theta = 0;
% lambda = 4 * pi/6;
% psi = 0;


sigma_x = sigma;
sigma_y = sigma/gamma;

nstds = 3;
xmax = max(abs(nstds*sigma_x*cos(theta)),abs(nstds*sigma_y*sin(theta)));
xmax = ceil(max(1,xmax));
ymax = max(abs(nstds*sigma_x*sin(theta)),abs(nstds*sigma_y*cos(theta)));
ymax = ceil(max(1,ymax));
xmin = -xmax; ymin = -ymax;
[x,y] = meshgrid(xmin:xmax,ymin:ymax);


x_theta=x*cos(theta)+y*sin(theta);
y_theta=-x*sin(theta)+y*cos(theta);

gb= exp(-.5*(x_theta.^2/sigma_x^2+y_theta.^2/sigma_y^2)).*cos(2*pi/lambda*x_theta+psi);

% figure(1);
% imshow(gb);
% 
% figure(2);
% [rows, cols] = size(gb);
% [x_image, y_image] = meshgrid(1:cols, 1:rows);
% 
% surf(x_image,y_image,gb)
% % surf(x_image,y_image,gb,gradient(gb))
% colorbar
% 
% figure(3);
% G2 = fft2(gb);
% imshow(log(abs(fftshift(G2)) + 1), [])

end
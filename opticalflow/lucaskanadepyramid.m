function [u3, v3] = lucaskanadepyramid(image1, image2)

[rows, cols] = size(image1);

%calculate lucas kanade for highest level
% [u,v] = lucaskanade(image, image2);

maxlevel = 2;
windowx = 2;
windowy = 2;


%% calculate pyramid for the images

% add boundary in the image
% image1new = [image1(1,1) image1(1,:) image1(1,end)
%     image1(:,1) image1 image1(:,end)
%     image1(end,1) image1(end,:) image1(end,end)];
% image2new = [image2(1,1) image2(1,:) image2(1,end)
%     image2(:,1) image2 image2(:,end)
%     image2(end,1) image2(end,:) image2(end,end)];
% 
% filterreduce = [0.0625 0.125 0.0625
%     0.125 0.25 0.125
%     0.0625 0.125 0.0625];
% height = floor((rows+1)/2);
% width = floor((cols+1)/2);
% imagea = zeros(height, width);
% imageb = zeros(height, width);
% 
% for i = 1 : height
%     for j = 1 : width
%         neighbourhood1 = image1new(2*i-1:2*i+1, 2*j-1:2*j+1);
%         answer1 = filterreduce.*neighbourhood1;
% %         disp(answer);
%         value1 = sum(sum(answer1)');
% %         disp(value);
%         imagea(i,j) = value1;
%         neighbourhood2 = image2new(2*i-1:2*i+1, 2*j-1:2*j+1);
%         answer2 = filterreduce.*neighbourhood2;
% %         disp(answer);
%         value2 = sum(sum(answer2)');
% %         disp(value);
%         imageb(i,j) = value2;
%     end
% end

imagea = pyramidreduce(image1);
imageb = pyramidreduce(image2);

%% calculate the optical flow at the deepest level and this wil be used as the initial guess
[u2,v2] = lucaskanade(imagea, imageb);

disp(u2);
disp(v2);

[rows2, cols2] = size(imagea);

%% now multiply by two so it belongs to the next level
u2 = 2 * u2;
v2 = 2 * v2;

% %% bilinearly interpolate the optical flow
% % get double resolution
u3 = zeros(rows, cols);
v3 = zeros(rows, cols);
for i = 1 : rows2 - 1 % - 1 since we are calculated the right as well
    for j = 1 : cols2 - 1
        %compute the values in all four directions for u and v
        u3(2*i-1,2*j-1) = u2(i,j);
        
        u3(2*i-1,2*j+1) = u2(i,j+1);
        u3(2*i+1,2*j-1) = u2(i+1,j);
        u3(2*i+1,2*j+1) = u2(i+1,j+1);
        
        v3(2*i-1,2*j-1) = v2(i,j);
        v3(2*i-1,2*j+1) = v2(i,j+1);
        v3(2*i+1,2*j-1) = v2(i+1,j);
        v3(2*i+1,2*j+1) = v2(i+1,j+1);
        
        %apply interpolation formula
        u3(2*i-1,2*j) = 0.5*u3(2*i-1,2*j-1) +  0.5*u3(2*i-1,2*j+1);
        u3(2*i,2*j-1) = 0.5*u3(2*i-1,2*j-1) +  0.5*u3(2*i+1,2*j-1);
        u3(2*i,2*j-1) = 0.25*u3(2*i-1,2*j-1) + 0.25*u3(2*i-1,2*j+1) + 0.25*u3(2*i+1,2*j-1) + 0.25*u3(2*i+1,2*j+1);
        u3(2*i,2*j+1) = 0.5*u3(2*i-1,2*j+1) +  0.5*u3(2*i+1,2*j+1);
        u3(2*i+1,2*j) = 0.5*u3(2*i+1,2*j-1) +  0.5*u3(2*i+1,2*j+1); 
    end
end




%% compute residual error
[uerror, verror] = lucaskanadeiterative(image1, image2, u3, v3);

%% add the corrections
u3 = u3 + uerror;
v3 = v3 + verror;

% figure;
% quiver(u3, v3);

end
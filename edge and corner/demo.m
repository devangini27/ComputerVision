img = imread('railway.jpg');
%greyim = img; %rgb2gray(img);
figure(1);
imshow(img);


g  = fspecial('gaussian',[7, 7], 1);

greyim = conv2(img, g, 'same');


odd = [255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0 , 0 , 0 , 0, 0, 0 , 0 , 0 , 0, 0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0 , 0 , 0 , 0, 0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255];
even = [0 , 0 , 0, 0, 0, 0 , 0 , 0 , 0, 0 , 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0 , 0 , 0 , 0, 0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0 , 0 , 0 , 0, 0];

greyim1 = [odd;
                 odd;
                 odd;
                 odd;
                 odd;
                 odd;
                 odd;
                 odd;
                 odd;
                 odd;
                 even;
                 even;
                 even;
                 even;
                 even;
                 even;
                 even;
                 even;
                 even;
                 even;
                 odd;
                 odd;
                 odd;
                 odd;
                 odd;
                 odd;
                 odd;
                 odd;
                 odd;
                 odd;
                 even;
                 even;
                 even;
                 even;
                 even;
                 even;
                 even;
                 even;
                 even;
                 even;
                 odd;
                 odd;
                 odd;
                 odd;
                 odd;
                 odd;
                 odd;
                 odd;
                 odd;
                 odd];

%greyim = 255 * imread('chessboard.gif');             

figure(2);
imshow(greyim);

[rows, cols] = size(greyim);
horEdges = zeros(rows, cols-1);
template = [1 ; -1];

%disp(template(1));
%disp(template(2));

% looping on the image
for i = 1: rows
    for j = 1 : cols-1
        % apply the template
        horEdges(i,j) = abs(greyim(i,j) * template(1) + greyim(i, j+1) * template(2));
    end
end

figure(3);
imshow(horEdges);
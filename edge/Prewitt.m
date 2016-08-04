img = imread('railway.jpg');
figure(1);
imshow(img);

g1 = fspecial('prewitt')

smoothim = conv2(double(img), g1, 'same');
figure(2);
imshow(smoothim);

g2 = g1'

smoothim2 = conv2(double(img), g2, 'same');
figure(3);
imshow(smoothim2);

magnitude = smoothim.*smoothim + smoothim2.*smoothim2;
magnitude = uint8(sqrt(magnitude));

figure(4);
imshow(magnitude);

threshold = 220;
binary = zeros(size(magnitude));
binary(find(magnitude>threshold)) = 1;

figure(5);
imshow(binary);

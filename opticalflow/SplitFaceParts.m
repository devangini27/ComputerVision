function [parts] = SplitFaceParts(imageIndex, imagePath, imageName)

% SplitFaceParts(3, 'SMIC//SMIC_all_cropped//NIR//s11//micro//surprise//s11_sur_01//', 'reg_64773.bmp')

image = imread([imagePath imageName]);
[rows, cols] = size(image);

%% FIND ALL THE VISUAL LANDMARKS OF THE FACE

facePoints = FacePoints(imageIndex, imagePath, imageName);


% facePoints = [-3.2716   52.5983
%    -1.9779   72.4269
%     0.4201   92.1464
%     3.9240  112.3180
%    10.1001  132.4715
%    21.3813  149.7296
%    36.5548  162.9320
%    53.9824  172.4256
%    73.0104  176.0611
%    92.2768  171.7534
%   109.5130  161.3252
%   124.3246  147.3520
%   135.5236  129.8400
%   142.1749  110.2683
%   145.9137   91.0715
%   148.1827   72.4357
%   149.1433   53.6391
%     8.7972   39.3157
%    17.9364   31.8611
%    30.0650   29.2373
%    42.6790   30.6556
%    54.8904   34.4591
%    91.8729   35.0963
%   103.4458   31.5965
%   115.4555   30.2973
%   127.0376   32.6144
%   135.8340   39.5636
%    73.8164   53.5849
%    73.6137   66.8170
%    73.3994   80.0122
%    73.1843   93.0537
%    60.8471  100.1645
%    67.0218  102.4791
%    73.5201  103.6706
%    79.8128  102.5788
%    85.7591  100.3518
%    26.2143   54.4752
%    34.0456   49.3136
%    43.6755   49.2791
%    51.4274   54.5721
%    43.0880   57.1679
%    34.2469   57.4060
%    94.8581   55.0591
%   102.4957   50.0220
%   111.9897   50.0391
%   119.5897   55.3769
%   111.6472   58.2083
%   102.9566   57.7959
%    46.3813  125.5736
%    54.5343  121.0871
%    63.5667  118.0074
%    73.0689  120.0332
%    82.5448  118.0798
%    91.3983  121.2242
%    99.2052  125.9878
%    91.7458  131.5489
%    82.6129  134.5950
%    73.1384  135.6845
%    63.5683  134.4807
%    54.2310  131.2853
%    63.6651  125.1433
%    73.1466  126.5207
%    82.6068  125.2458
%    82.5981  124.3921
%    73.2578  125.1124
%    63.8512  124.3276
% ];


%% FIND ALL THE PARTS OF THE DETECTED LANDMARKS E.G. EYEBROWS, JAW, EYES, NOSE, OUTER MOUTH, INNER MOUTH.

figure(1);
imshow(image);
hold on;

for i  = 1 : 66
    plot(facePoints(:,1), facePoints(:,2),'*');
end

%plot(facePoints(1,1), facePoints(1,2),'r*');

jawStart = 1; % 1-17
jawLength = 17;
jaw = facePoints(jawStart:jawStart + jawLength - 1,:);
plot(jaw(:,1), jaw(:,2));

eyebrow1Start = 18; % 18-22
eyebrow1Length = 5;
eyebrow1 = facePoints(eyebrow1Start:eyebrow1Start + eyebrow1Length - 1,:);
plot(eyebrow1(:,1), eyebrow1(:,2));

eyebrow2Start = 23; % 23-27
eyebrow2Length = 5;
eyebrow2 = facePoints(eyebrow2Start:eyebrow2Start + eyebrow2Length - 1,:);
plot(eyebrow2(:,1), eyebrow2(:,2));

noseStart = 28; % 28-36
noseLength = 9;
nose = facePoints(noseStart:noseStart + noseLength - 1,:);
plot(nose(:,1), nose(:,2));

eye1Start = 37; % 37-42
eye1Length = 6;
eye1 = facePoints(eye1Start:eye1Start + eye1Length - 1, :);
plot(eye1(:,1), eye1(:,2));

eye2Start = 43; % 43-48
eye2Length = 6;
eye2 = facePoints(eye2Start:eye2Start + eye2Length - 1, :);
plot(eye2(:,1), eye2(:,2));

outerMouthStart = 49; %49-60
outerMouthLength = 12;
outerMouth = facePoints(outerMouthStart:outerMouthStart + outerMouthLength - 1, :);
plot(outerMouth(:,1), outerMouth(:,2));

innerMouthStart = 61; %61-66
innerMouthLength = 6;
innerMouth = facePoints(innerMouthStart:innerMouthStart + innerMouthLength - 1, :);
plot(innerMouth(:,1), innerMouth(:,2));

%% FIND ALL THE SPECIAL POINTS ON THE FACE E.G. MIDDLE OF THE EYEBROWS, TEARDUCTS, NOSE TIP, CHIN, 

figure(2);
imshow(image);
hold on;

eyebrowMiddle = [ (eyebrow1(eyebrow1Length, 1) + eyebrow2(1,1))/2 , (eyebrow1(1, 2) + eyebrow2(eyebrow2Length,2))/2  ];
plot(eyebrowMiddle(:,1), eyebrowMiddle(:,2),'*');

tearduct1 = eye1(eye1Length - 2, :);
plot(tearduct1(:,1), tearduct1(:,2),'*');

tearduct2 = eye2(1, :);
plot(tearduct2(:,1), tearduct2(:,2),'*');

nosemiddle = nose(2,:);
plot(nosemiddle(:,1), nosemiddle(:,2),'*');

nosetip = nose(noseLength - 2, :);
plot(nosetip(:,1), nosetip(:,2),'*');

chin1 = jaw(9, :);
plot(chin1(:,1), chin1(:,2),'*');

chin2 = jaw(8, :);
plot(chin2(:,1), chin2(:,2),'*');

chin3 = jaw(10,:);
plot(chin3(:,1), chin3(:,2),'*');

chin4 = [  chin1(1),    (chin2(2) + chin3(2))/2];
plot(chin4(:,1), chin4(:,2),'*');

cheekp1 = jaw(6,:);
plot(cheekp1(:,1), cheekp1(:,2),'*');

cheekp2 = jaw(12,:);
plot(cheekp2(:,1), cheekp2(:,2),'*');

%% FIND ALL THE BOX COMPONENTS OF THE FACE E.G. EYEBROWS, EYES, NOSE, CHEEKS, MOUTH.

figure(3);
imshow(image);
hold on;

faceclipx = 10;
faceStart = 1 + faceclipx;
chinclipy = 10;

eyebrowbox1 = [faceStart,   10,     eyebrowMiddle(1) - faceStart,    eyebrowMiddle(2)-10];
rectangle('Position',eyebrowbox1);

eyebrowbox2 = [eyebrowMiddle(1),    10,     jaw(jawLength, 1) - eyebrowMiddle(1) - faceclipx,    eyebrowMiddle(2)-10];
rectangle('Position',eyebrowbox2);

eyebox1 = [faceStart,   eyebrowMiddle(2),   tearduct1(1) - faceStart,    nosemiddle(2) - eyebrowMiddle(2)];
rectangle('Position',eyebox1);

eyebox2 = [tearduct2(1),    eyebrowMiddle(2),  jaw(jawLength, 1) - tearduct2(1) - faceclipx,   nosemiddle(2) - eyebrowMiddle(2)];
rectangle('Position',eyebox2);

nose = [tearduct1(1),   eyebrowMiddle(2),   tearduct2(1) - tearduct1(1),    nosetip(2) - eyebrowMiddle(2)];
rectangle('Position',nose);

mouth = [chin2(1),   nosetip(2),    chin3(1) - chin2(1),    chin4(2) - nosetip(2) ];
rectangle('Position',mouth);

cheek1 = [cheekp1(1),   nosemiddle(2) + 5,    chin2(1) - cheekp1(1),    chin4(2) - nosemiddle(2) - 20 - chinclipy ];
rectangle('Position',cheek1);

cheek2 = [chin3(1),   nosemiddle(2) + 5,    cheekp2(1) - chin3(1),    chin4(2) - nosemiddle(2) - 20 - chinclipy ];
rectangle('Position',cheek2);


parts = uint8([ eyebrowbox1
    eyebrowbox2
    eyebox1
    eyebox2
    nose
    mouth
    cheek1
    cheek2]);

end
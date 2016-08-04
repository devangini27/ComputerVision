img = imread('railway.jpg');
figure(1);
imshow(img);

MH1 = edge(img,'log',0,1.0);
MH2 = edge(img,'log',0,2.0);
MH3 = edge(img,'log',0,3.0);
MH4 = edge(img,'log',0,4.0);

% form mosaic
EFGH = [ MH1 MH2; MH3 MH4];

%% show mosaic in Matlab Figure window
log = figure('Name','Marr/Hildreth: UL: s=1  UR: s=2  BL: s=3 BR: s=4');
iptsetpref('ImshowBorder','tight');
imshow(EFGH,'InitialMagnification',100);
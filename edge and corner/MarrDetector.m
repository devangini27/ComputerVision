img = imread('railway.jpg');
figure(1);
imshow(img);

g1 = fspecial('gaussian', [7 7], 1);
%surf(g1)

smoothim = uint8(conv2(double(img), g1, 'same'));
figure(2);
imshow(smoothim);

lx = fspecial('laplacian',0)
laplacianx = uint8(conv2(double(smoothim), lx, 'same'));
figure(3);
imshow(laplacianx);

ly = [1 1 1
    1 -8 1
    1 1 1]
laplaciany = uint8(conv2(double(smoothim), ly, 'same'));
figure(3);
imshow(laplaciany);

laplacian = laplacianx + laplaciany;


g3 = conv2(g1, lx,'same')


g4  = [0 0 1 0 0;
       0 1 2 1 0;
       1 2 -16 2 1;
       0 1 2 1 0;
       0 0 1 0 0];
   
log = fspecial('log', [7 7],1);
   
laplacian2 =  uint8(conv2(double(img), log, 'same'));
figure(4);
imshow(laplacian2);


[rows, cols] = size(img);

count = 0;
for i = 1:rows
    for j = 1:cols
        if laplacian2(i,j) < 0
            count=count+1
        end
    end
end

%{
threshold = 20;
slopem = zeros(rows,cols);
binary = zeros(rows,cols);

for i = 1:rows
    for j = 1:cols-1
        
        slope = 0;
        a = laplacian(i,j);
        b = laplacian(i,j+1);
        
        if (a>0 && b< 0 ) || (a< 0 && b > 0)
            slope = abs(a-b);
            if slope > threshold
                slopem(i,j) = slope;
                binary(i,j) = 255;
            end
        elseif b == 0 && j ~= cols-1
            c= smoothim2(i,j+2);
            if (a>0 && c < 0) || (a< 0 && c > 0)
                slope = abs(a-c);
                if slope > threshold
                    slopem(i,j+1) = slope;
                    binary(i,j+1) = 255;
                end
            end
        end
        
    end
end

figure(4);
imshow(slopem);
%}

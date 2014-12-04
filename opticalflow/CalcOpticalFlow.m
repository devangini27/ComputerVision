function CalcOpticalFlow(imagepath, imagearray, count, subjectString, microString, expressionString)

algorithmTypes = ['HS','LK','LKP'];
opticalflowOption = 1;

mkdir([algorithmTypes(opticalflowOption) '//images//' subjectString]);
mkdir([algorithmTypes(opticalflowOption) '//images//' subjectString '//' microString]);
mkdir([algorithmTypes(opticalflowOption) '//images//' subjectString '//' microString '//'  expressionString]);
mkdir([algorithmTypes(opticalflowOption) '//images//' subjectString '//' microString '//'  expressionString '//flow']);
mkdir([algorithmTypes(opticalflowOption) '//images//' subjectString '//' microString '//' expressionString '//filtered']);
mkdir([algorithmTypes(opticalflowOption) '//images//' subjectString '//' microString '//' expressionString '//resultant']);

im1 = imread([imagepath imagearray(1,:) ]);
% parts = SplitFaceParts(1, imagepath, imagearray(1,:));


for index = 2 : count
    
    im2 = imread([imagepath imagearray(index,:) ]);
    [u,v] = hornschunck(im1, im2);
    uv = [u,v];
    
    % Display estimated flow fields
    h=figure; subplot(1,2,1);imshow(uint8(flowToColor(uv))); title('Middlebury color coding');
    subplot(1,2,2); plotflow(uv);   title('Vector plot');
    
    saveas(h,[algorithmTypes(opticalflowOption) '//images//' subjectString '//' microString '//'  expressionString  '//flow//flow' num2str(index)],'jpg')
  
    % [flow2] = CalcResultantOpticalFlow(uv, parts, index, subjectString, microString, expressionString);
    im1 = im2;
    
    
    
end
end
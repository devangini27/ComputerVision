% DECLARE ALL THE FOLDER NAMES SO THAT THEY CAN BE EASILY USED FOR PATH CONSTRUCTION

addpath('.//flowColorCode');
addpath('.//blacklib');

PARENTFOLDER = '..//..//micro database//';
MAINIMAGEFOLDERS1 = [PARENTFOLDER 'SMIC//SMIC_all_cropped'; PARENTFOLDER  'SMIC//SMIC_all_raw    '];
MAINIMAGEFOLDERS = cellstr(MAINIMAGEFOLDERS1);
MAINIMAGEINDEX = 1;

LIGHTTYPESFOLDERS1 = ['HS '; 'NIR'; 'VIS']; % HS = High speed,  VIS = visible light, %NIR = near infrared
LIGHTTYPESFOLDERS = cellstr(LIGHTTYPESFOLDERS1);
LIGHTTYPEINDEX = 2;

SUBJECTSCOUNT = 20;
SUBJECTSFOLDERS1 = ['s1 '; 's2 '; 's3 '; 's4 '; 's5 '; 's6 '; 's7 '; 's8 '; 's9 '; 's10'; 's11'; 's12'; 's13'; 's14'; 's15'; 's16'; 's17'; 's18'; 's19'; 's20'];
SUBJECTSFOLDERS = cellstr(SUBJECTSFOLDERS1);

MICROCOUNT = 2;
MICROFOLDERS1 = ['micro    '; 'non_micro'];
MICROFOLDERS = cellstr(MICROFOLDERS1);
microIndex = 1;

EXPRESSIONCOUNT = 3;
EXPRESSIONTYPES1 = ['negative'; 'positive'; 'surprise'];
EXPRESSIONTYPES = cellstr(EXPRESSIONTYPES1);
expressionIndex = 3;


%read one sequence in all the expressions
breakCondition = false;

% read all the images for subjects
for subjectIndex = 1 : SUBJECTSCOUNT % 11%1: subjectsCount 18,11
    
    for microIndex = 1 : MICROCOUNT%1, 1: MICROCOUNT , 1
        
        for expressionIndex = 1 : EXPRESSIONCOUNT % 1 : EXPRESSIONCOUNT, 3
            
            
            
            
        
            path = '';
            if microIndex == 1
                path = char(strcat( MAINIMAGEFOLDERS{MAINIMAGEINDEX} , '//' , LIGHTTYPESFOLDERS{LIGHTTYPEINDEX}, '//', SUBJECTSFOLDERS{subjectIndex} , '//' , MICROFOLDERS{microIndex} , '//',  EXPRESSIONTYPES{expressionIndex}));
            else
                path = char(strcat( MAINIMAGEFOLDERS{MAINIMAGEINDEX} , '//' , LIGHTTYPESFOLDERS{LIGHTTYPEINDEX}, '//', SUBJECTSFOLDERS{subjectIndex} , '//' , MICROFOLDERS{microIndex}));
            end
            %disp(path);
            
            %check if this path exists
            if isequal(exist(path, 'dir'),7) % 7 = directory.
                %display('a folder');
                disp(path);
                %check if there are any samples
                listing = dir(path);
                [samplesSize, ~] = size(listing);
                
                if samplesSize > 2 % THE FIRST TWO ENTRIES ARE . AND ..
                    
                    folderIndex = 1;
                    folderName = char(strcat(path, '//', listing(2 + folderIndex).name)); % THE FIRST TWO ENTRIES ARE . AND ..
                    
                    imageListing = dir(folderName);
                    [imageCount, ~] = size(imageListing);
                    imagePath = char(strcat(folderName, '//')) % , imageName
                    
                    imagesList = '';
                    imagearray = [];
                    
                    for imageIndex = 3 : imageCount % THE FIRST TWO ENTRIES ARE . AND ..
                        imageName = imageListing(imageIndex).name
                        imagearray = [imagearray; imageName];
                    end
                    
                    
                    CalcOpticalFlow(imagePath, imagearray, imageCount-2, SUBJECTSFOLDERS{subjectIndex} , MICROFOLDERS{microIndex} ,  EXPRESSIONTYPES{expressionIndex});
                    
                    
                    breakCondition = true;
                    break;
                end
            end
            
            if breakCondition
                break
            end
            
            
        end
        
        
    end
end



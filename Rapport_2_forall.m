clear all, close all, clc
%%
% Loading the data
addpath(genpath('data'));
addpath(genpath('Matlab'));
%%
[multiIm, annotationIm] = loadMulti('multispectral_day13.mat','annotation_day13.png');

[fatPix, fatR, fatC] = getPix(multiIm, annotationIm(:,:,2));
[meatPix, meatR, meatC] = getPix(multiIm, annotationIm(:,:,3));

% Converting from int8 to double.
multiImDouble = double(multiIm);

% Combing annotations to get image of background. 
background = sum(annotationIm,3);

fatSigma = cov(fatPix);
meatSigma = cov(meatPix);

meanFat = mean(fatPix);
meanMeat = mean(meatPix);

lenMeat = length(meatPix);
lenFat = length(fatPix);

pooledSigma = (1/((lenMeat-1)+(lenFat-1))).*((lenMeat-1).*meatSigma+(lenFat-1).*fatSigma);
pooledSigmaInv = inv(pooledSigma);

% Defining the S-functions
Sf_fat = @(x) permute(x,[1,3,2])*pooledSigmaInv*permute(meanFat,[2,1])-(1/2)*meanFat*pooledSigmaInv*permute(meanFat,[2,1]);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
Sf_meat = @(x) permute(x,[1,3,2])*pooledSigmaInv*permute(meanMeat,[2,1])-(1/2)*meanMeat*pooledSigmaInv*permute(meanMeat,[2,1]);

%% Doing LDA and error on the other days

days = [1,6,13,20,28];
absErrors = zeros(1,5);
numPixels = zeros(1,5);

for k = 1:5
    [multiIm, annotationIm] = loadMulti(strcat('multispectral_day',sprintf('%02d',days(k)),'.mat'),strcat('annotation_day',sprintf('%02d',days(k)),'.png'));
    [fatPix, fatR, fatC] = getPix(multiIm, annotationIm(:,:,2));
    [meatPix, meatR, meatC] = getPix(multiIm, annotationIm(:,:,3));

    multiImDouble = double(multiIm);
 
    % Computing S-values for every pixel.
    
    Sfat = zeros(514);
    Smeat = zeros(514);

    for i = 1:514
        for j = 1:514
            Sfat(i,j) = Sf_fat(multiImDouble(i,j,:));
            Smeat(i,j) = Sf_meat(multiImDouble(i,j,:));
        end
    end

    Sdif = Smeat./Sfat;
    
    % Classifying pixels using advanced model. 

    sClass = zeros(514);

    for i = 1:514
        for j = 1:514
            if Sdif(i,j) > 1 && background(i,j) == 1
                sClass(i,j) = 1;
            end
        end
    end
    
    subplot(2,3,k);
    imshow(sClass);
    title(strcat('Day',{' '}, int2str(days(k))));
    sgtitle('Classifiation');
    
    % Counting errors using advanced model.
    
    for i = 1:514
        for j = 1:514
            if annotationIm(i,j,3) == 1 || annotationIm(i,j,2) == 1
                numPixels(k) = numPixels(k) + 1;
            end
            if (annotationIm(i,j,3) == 1 && sClass(i,j) == 0) || (annotationIm(i,j,2) == 1 && sClass(i,j) == 1)
                absErrors(k) = absErrors(k) + 1;
            end
        end
    end
end

%%

errorRate = absErrors./numPixels;
figure(2)
plot(days,errorRate);

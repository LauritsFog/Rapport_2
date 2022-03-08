addpath(genpath('data'));
addpath(genpath('Matlab'));

%%

[multiIm, annotationIm] = loadMulti('multispectral_day01.mat','annotation_day01.png');

[fatPix, fatR, fatC] = getPix(multiIm, annotationIm(:,:,2));
[meatPix, meatR, meatC] = getPix(multiIm, annotationIm(:,:,3));

% Finding spectral band with best discriminative properties for meat and
% fat. Using mean and median.

meanDif = abs(mean(fatPix)-mean(meatPix));

[max, idx] = max(meanDif);

% Calculating threshold values for day 1. 

meanThresholds = (mean(fatPix)+mean(meatPix))/2;

%%

figure(1)
plot(mean(meatPix), 'b');
hold on
plot(mean(fatPix), 'r');
title('Mean fatPix vs mean meatPix');

%%

errorRate = zeros(1,19);

for l = 1:18
    for p = 1:756
        if fatPix(p,l) < meanThresholds(l)
            errorRate(l) = errorRate(l) + 1;
        end
        if meatPix(p,l) > meanThresholds(l)
            errorRate(l) = errorRate(l) + 1;
        end
    end
end

%%

figure(2)
plot(errorRate);
title('Error rate');

%%

fatClass = zeros(514);

% Classifying every meat-pixel and giving it value 1. 

for i = 1:514
    for j = 1:514
        if multiIm(i,j,idx) < meanThresholds(idx) && background(i,j) == 1
            fatClass(i,j) = 1;
        end
    end
end

%%

figure(3);
imshow(fatClass);
title('Meat-pixels classified and colored white with simple model');

%%

fatSigma = cov(fatPix);
meatSigma = cov(meatPix);

meanFat = mean(fatPix);
meanMeat = mean(meatPix);

lenMeat = length(meatPix);
lenFat = length(fatPix);

pooledSigma = (1/((lenMeat-1)+(lenFat-1))).*((lenMeat-1).*meatSigma+(lenFat-1).*fatSigma);

pooledSigmaInv = inv(pooledSigma);

%%

% Converting from int8 to double.

multiImDouble = double(multiIm);

%%

% Combing annotations to get image of background. 

background = sum(annotationIm,3);

%%

% Defining the S-functions. 

Sf_fat = @(x) permute(x,[1,3,2])*pooledSigmaInv*permute(meanFat,[2,1])-(1/2)*meanFat*pooledSigmaInv*permute(meanFat,[2,1]);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
Sf_meat = @(x) permute(x,[1,3,2])*pooledSigmaInv*permute(meanMeat,[2,1])-(1/2)*meanMeat*pooledSigmaInv*permute(meanMeat,[2,1]);

%%

% Computing Sfat and Smeat for each pixel. 

Sfat = zeros(514);
Smeat = zeros(514);

for i = 1:514
    for j = 1:514
        Sfat(i,j) = Sf_fat(multiImDouble(i,j,:));
        Smeat(i,j) = Sf_meat(multiImDouble(i,j,:));
    end
end

%%

% Computing difference. If Sdif > 1 the probability for meat is largest. 

Sdif = Smeat./Sfat;

%%

% Coloring in every meat pixel. 

sClassDay1 = zeros(514);

for i = 1:514
    for j = 1:514
        if Sdif(i,j) > 1 && background(i,j) == 1
            sClassDay1(i,j) = 1;
        end
    end
end

%% 

figure(4);
imshow(sClassDay1);
title('Meat-pixels classified and colored white with advanced model');

%%

errorRateDay1 = 0;

for i = 1:514
    for j = 1:514
        if sClassDay1(i,j) ~= annotationIm(i,j,3) && annotationIm(i,j,3) == 1
            errorRateDay1 = errorRateDay1 + 1;
        end
    end
end

%%

days = [1,6,13,20,28];
absErrors = zeros(1,5);
numPixels = zeros(1,5);

for k = 1:5
    [multiIm, annotationIm] = loadMulti(strcat('multispectral_day',sprintf('%02d',days(k)),'.mat'),strcat('annotation_day',sprintf('%02d',days(k)),'.png'));
    [fatPix, fatR, fatC] = getPix(multiIm, annotationIm(:,:,2));
    [meatPix, meatR, meatC] = getPix(multiIm, annotationIm(:,:,3));

    multiImDouble = double(multiIm);

    background = sum(annotationIm,3);

    meanDif = abs(mean(fatPix)-mean(meatPix));

    [max, idx] = max(meanDif);

    meanThresholds = (mean(fatPix)+mean(meatPix))/2;
    
    absErrorSimple = zeros(1,19);
    
    % Counting errors using simple model.
    
    for l = 1:18
        for p = 1:756
            if fatPix(p,l) < meanThresholds(l) || meatPix(p,l) > meanThresholds(l)
                absErrorSimple(l) = absErrorSimple(l) + 1;
            end
        end
    end
    
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

%%

figure(5)
bar(days,absErrors);
figure(6)
bar(days,errorRate);

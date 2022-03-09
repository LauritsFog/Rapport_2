addpath(genpath('data'));
addpath(genpath('Matlab'));

%%

[multiIm, annotationIm] = loadMulti('multispectral_day01.mat','annotation_day01.png');

[fatPix, fatR, fatC] = getPix(multiIm, annotationIm(:,:,2));
[meatPix, meatR, meatC] = getPix(multiIm, annotationIm(:,:,3));

% Finding spectral band with best discriminative properties for meat and
% fat. Using mean and median.

meanDif = abs(mean(fatPix)-mean(meatPix));

[maxVal, idx] = max(meanDif);

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

% Combining annotations to get image of background. 

background = sum(annotationIm,3);

%%

% Converting from int8 to double.

multiImDouble = double(multiIm);

%%

meatClass = zeros(514);

% Classifying every meat-pixel and giving it value 1. 

for i = 1:514
    for j = 1:514
        if multiImDouble(i,j,idx) < meanThresholds(idx) && background(i,j) == 1
            meatClass(i,j) = 1;
        end
    end
end

%%

figure(3);
imshow(meatClass);
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

advClassDay1 = zeros(514);

for i = 1:514
    for j = 1:514
        if Sdif(i,j) > 1 && background(i,j) == 1
            advClassDay1(i,j) = 1;
        end
    end
end

%% 

figure(4);
imshow(advClassDay1);
title('Meat-pixels classified and colored white with advanced model');

%%

errorRateDay1 = 0;

for i = 1:514
    for j = 1:514
        if advClassDay1(i,j) ~= annotationIm(i,j,3) && annotationIm(i,j,3) == 1
            errorRateDay1 = errorRateDay1 + 1;
        end
    end
end

%%

days = [1,6,13,20,28];
absErrorsAdv = zeros(1,5);
absErrorsSimple = zeros(1,5);
numPixels = zeros(1,5);

for k = 1:5
    [multiIm, annotationIm] = loadMulti(strcat('multispectral_day',sprintf('%02d',days(k)),'.mat'),strcat('annotation_day',sprintf('%02d',days(k)),'.png'));
    
    multiImDouble = double(multiIm);

    background = sum(annotationIm,3);
    
    simpleClass = zeros(514);
    
    % Classifying using simple model. 1 = meat. 
    
    for i = 1:514
        for j = 1:514
            if multiImDouble(i,j,idx) < meanThresholds(idx) && background(i,j) == 1
                simpleClass(i,j) = 1;
            end
        end
    end
    
    figure(5)
    subplot(2,3,k);
    imshow(simpleClass);
    title(strcat('Day',{' '}, int2str(days(k))));
    sgtitle('Classification using simple model');
    
    % Computing S-values for every pixel.
    
    Sfat = zeros(514);
    Smeat = zeros(514);

    for i = 1:514
        for j = 1:514
            if background(i,j) == 1
                Sfat(i,j) = Sf_fat(multiImDouble(i,j,:));
                Smeat(i,j) = Sf_meat(multiImDouble(i,j,:));
            end
        end
    end

    Sdif = Smeat./Sfat;
    
    % Classifying pixels using advanced model. 1 = meat. 

    advClass = zeros(514);

    for i = 1:514
        for j = 1:514
            if Sdif(i,j) > 1 && background(i,j) == 1
                advClass(i,j) = 1;
            end
        end
    end
    
    figure(6)
    subplot(2,3,k);
    imshow(advClass);
    title(strcat('Day',{' '}, int2str(days(k))));
    sgtitle('Classification using advanced model');
    
    for i = 1:514
        for j = 1:514
            if annotationIm(i,j,3) == 1 || annotationIm(i,j,2) == 1
                numPixels(k) = numPixels(k) + 1;
            end
            % Counting errors using simple model.
            if (annotationIm(i,j,3) == 1 && simpleClass(i,j) == 0) || (annotationIm(i,j,2) == 1 && simpleClass(i,j) == 1)
                absErrorsSimple(k) = absErrorsSimple(k) + 1;
            end
            % Counting errors using advanced model.
            if (annotationIm(i,j,3) == 1 && advClass(i,j) == 0) || (annotationIm(i,j,2) == 1 && advClass(i,j) == 1)
                absErrorsAdv(k) = absErrorsAdv(k) + 1;
            end
        end
    end
end

%%

errorRateAdv = absErrorsAdv./numPixels;
errorRateSimple = absErrorsSimple./numPixels;

%%

figure(7)
plot(days,errorRateAdv,'b');
hold on
plot(days,errorRateSimple,'r');

%%


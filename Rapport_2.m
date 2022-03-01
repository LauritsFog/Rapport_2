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

figure
plot(mean(meatPix), 'b');
hold on
plot(mean(fatPix), 'r');

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

plot(errorRate);

%%

fatClass = zeros(514);

% Classifying every fat-pixel and giving it value 1. 

for i = 1:514
    for j = 1:514
        if multiIm(i,j,idx) < meanThresholds(idx)
            fatClass(i,j) = 1;
        end
    end
end

imagesc(fatClass);

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

multiImDouble = double(multiIm);

%%

Sf_fat = @(x) permute(x,[1,3,2])*pooledSigmaInv*permute(meanFat,[2,1])-(1/2)*meanFat*pooledSigmaInv*permute(meanFat,[2,1]);

Sf_meat = @(x) permute(x,[1,3,2])*pooledSigmaInv*permute(meanMeat,[2,1])-(1/2)*meanMeat*pooledSigmaInv*permute(meanMeat,[2,1]);

%%

Sfat = zeros(514);
Smeat = zeros(514);

for i = 1:514
    for j = 1:514
        Sfat(i,j) = Sf_fat(multiImDouble(i,j,:));
        Smeat(i,j) = Sf_meat(multiImDouble(i,j,:));
    end
end

%%

Sdif = Smeat./Sfat;

%%

sClass = zeros(514);

for i = 1:514
    for j = 1:514
        if Sdif(i,j) > 1
            sClass(i,j) = 1;
        end
    end
end

imagesc(fatClass);
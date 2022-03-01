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
        
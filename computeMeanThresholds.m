function [meanThresholds,idx] = computeMeanThresholds(multiIm,annotationIm)
    [fatPix, fatR, fatC] = getPix(multiIm, annotationIm(:,:,2));
    [meatPix, meatR, meatC] = getPix(multiIm, annotationIm(:,:,3));

    % Finding spectral band with best discriminative properties for meat and
    % fat. Using mean and median.

    meanDif = abs(mean(fatPix)-mean(meatPix));

    [maxVal, idx] = max(meanDif);

    % Calculating threshold values for day 1. 

    meanThresholds = (mean(fatPix)+mean(meatPix))/2;
end
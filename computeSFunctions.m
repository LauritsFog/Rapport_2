function[Sf_fat, Sf_meat] = computeSFunctions(multiIm,annotationIm,pFat,pMeat)
    
    pFat = pFat;
    pMeat = pMeat;

    [fatPix, fatR, fatC] = getPix(multiIm, annotationIm(:,:,2));
    [meatPix, meatR, meatC] = getPix(multiIm, annotationIm(:,:,3));

    fatSigma = cov(fatPix);
    meatSigma = cov(meatPix);

    meanFat = mean(fatPix);
    meanMeat = mean(meatPix);

    lenMeat = length(meatPix);
    lenFat = length(fatPix);

    pooledSigma = (1/((lenMeat-1)+(lenFat-1))).*((lenMeat-1).*meatSigma+(lenFat-1).*fatSigma);

    pooledSigmaInv = inv(pooledSigma);

    Sf_fat = @(x) permute(x,[1,3,2])*pooledSigmaInv*permute(meanFat,[2,1])-(1/2)*meanFat*pooledSigmaInv*permute(meanFat,[2,1])+log(pFat);

    Sf_meat = @(x) permute(x,[1,3,2])*pooledSigmaInv*permute(meanMeat,[2,1])-(1/2)*meanMeat*pooledSigmaInv*permute(meanMeat,[2,1])+log(pMeat);
end
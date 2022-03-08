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

figure (1)
hold on 
plot(mean(meatPix), 'b');
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

figure (2)
hold on 
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

figure (3) 
hold on 
imagesc(fatClass);

%% 

sigma_meat = cov(meatPix);
sigma_fat = cov(fatPix);

m=[969,756];

sum_under=0;

for i=1:2
    sum_under=sum_under+m(i)-1;
end

Sigma_pooled = ((m(1)-1)*sigma_meat+(m(2)-1)*sigma_fat)./sum_under ;


%%

s_meat=zeros(514);
s_fat=zeros(514);

mean2=[mean(meatPix)', mean(fatPix)'];

%lp1=3; 

M=im2double(multiIm);


for j=1:514
    for k=1:514
            s_meat(j,k)=squeeze(M(j,k,:))'*(Sigma_pooled^-1)*mean2(:,1)-1/2*mean2(:,1)'*(Sigma_pooled^-1)*mean2(:,1);
    end
end

for j=1:514
    for k=1:514
            s_fat(j,k)=squeeze(M(j,k,:))'*(Sigma_pooled^-1)*mean2(:,2)-1/2*mean2(:,2)'*(Sigma_pooled^-1)*mean2(:,2);
    end
end

 for i=1:2

   s(i)=squeeze(M(200,200,:))'*(Sigma_pooled^-1)*mean2(:,i)-((1/2)*mean2(:,i))'*(Sigma_pooled^-1)*mean2(:,i);
            
 end






        
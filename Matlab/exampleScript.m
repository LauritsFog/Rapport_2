% Dette script er et eksemple på brug af funktioner til første opgave i
% matematisk modllering - 02526
% 
% Anders Bjorholm Dahl
% Februar 2015
% 

% Folder hvor mine datafiler ligger
% dirIn = '/Users/abd/Documents/Teaching/02526_matematisk_modellering/Exercises/data/toCampusnet/';


%% Eksempel på at loade multispektralt billede

[multiIm, annotationIm] = loadMulti('multispectral_day01.mat' , 'annotation_day01.png');
 
% multiIm er et multispektralt billede - her ses dimensioner
size(multiIm)

%% her vises det billed ved spektralbånd 7

figure
imagesc(multiIm(:,:,7))

%% annotationIm er et binært billede med tre dimensioner

size(annotationIm)

%% I hver dimesion er et binært billede 1 - baggrund med pølse, 2 - fedtannotering, 
% 3 - kødannotering. Her vises kødannoteringen

figure
imagesc(annotationIm(:,:,1))

%% Funktionen getPix udtrækker multispektrale pixels for annoteringerne 

% Her er et eksemple på kød- og fedtannoteringerne
[fatPix, fatR, fatC] = getPix(multiIm, annotationIm(:,:,2));
[meatPix, meatR, meatC] = getPix(multiIm, annotationIm(:,:,3));

% Her plottes middelværdierne for pixels med kød hhv. fedt
figure
plot(mean(meatPix), 'b');
hold on
plot(mean(fatPix), 'r');


%% Funktionen showHistogram laver et histogram som den returnerer som en
% vektor

% Her er et eksempel - sidste argument gør at funktionen plotter
% histogrmmet for kød og fedt
h = showHistograms(multiIm, annotationIm(:,:,2:3), 3, 1);
axis square

%% histogrammet er også givet i h - h er altid 256 dimensioner, og da der
% kun er lave pixelværdier viser jeg kun de første 50

figure
plot(h(1:50,:))

%% Funktionen setImagePix giver et farvebillede, hvor pixel-koordinater
% gives som argument.

% RGB billede loades
imRGB = imread(['color_day20.png']);

% Pixel-koordinater for fedtannoteringerne fås med getPix
[fatPix, fatR, fatC] = getPix(multiIm, annotationIm(:,:,2));

% Pixel-koordinater sættes sammen til en matrix
pixId = [fatR, fatC];

% Det nye billede laves
rgbOut = setImagePix(imRGB, pixId);
figure
imagesc(rgbOut)



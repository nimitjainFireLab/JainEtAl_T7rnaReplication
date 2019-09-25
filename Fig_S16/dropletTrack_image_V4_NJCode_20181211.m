clc;
close all;
clear all;
numberPlots=4;
maindir='DD/';
picturedir_bright=[maindir,'bright/'];
picturedir_fluor=[maindir,'blue/'];
filename_bright = 'well 2_bright 4x nd1 p2s exp 7.tif'; %finds the drops
filename_fluor = 'well 2_blue 4x nd1 p2s exp 7.tif';  %gives the fluorescence 
imag = im2double(imread([picturedir_bright,filename_bright]));
lowerCircleSizeBright=6;
lowerCircleSizeFluor=8;
upperCircleSize=20;
sensitivityToUse=0.92;
[centersBright,radiiBright]=imfindcircles(imag,[lowerCircleSizeBright upperCircleSize],'Sensitivity',sensitivityToUse);
figure(1)
subplot(2,2,1)
imshow(imag,[])
subplot(2,2,2)
imshow(imag,[])
viscircles(centersBright,radiiBright);

imag_fluor = im2double(imread([picturedir_fluor,filename_fluor]));
a=imadjust(imag_fluor,[0;0.01]);
[centersFluor,radiiFluor]=imfindcircles(a,[lowerCircleSizeFluor upperCircleSize],'Sensitivity',sensitivityToUse,'EdgeThreshold',0.02);
subplot(2,2,3)
imshow(a,[])
subplot(2,2,4)
imshow(a,[])
viscircles(centersFluor,radiiFluor);

display(length(centersBright))
display(length(centersFluor))
display((length(centersFluor)/length(centersBright))*100)


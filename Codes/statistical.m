function [con, cor, ene, hom, meanx,stdx] = statistical(imaggex) %This function calculates the statistical properties of the image
statsc = graycoprops((imaggex),'Contrast Correlation Energy Homogeneity');
con = statsc.Contrast;
cor = statsc.Correlation;
ene = statsc.Energy;
hom = statsc.Homogeneity;
meanx = mean2(imaggex); % Indicates the brightness of the image..
stdx = std2(imaggex); 
%Ref: https://www.mathworks.com/help/images/ref/graycoprops.html
%This function contributes to a 5 dimensional feature vector. We
%observed that con and cor when added together and fed into the training and testing phase, increases the accuracy of the classifier. 
end

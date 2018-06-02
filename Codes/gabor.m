function [feat] = gabor(imagexx,gaborArray)
[q w] = size(imagexx);
feat = (gaborFeatures(imagexx,gaborArray,q/2 - 1,w/2 - 1))'; %Returns a 36 dimensional feature vector for each image. 
end
%% Reference:
%   M. Haghighat, S. Zonouz, M. Abdel-Mottaleb, "CloudID: Trustworthy 
%   cloud-based and cross-enterprise biometric identification," 
%   Expert Systems with Applications, vol. 42, no. 21, pp. 7905-7916, 2015.
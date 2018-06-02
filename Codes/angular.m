function [movingReg] = angular(fixed) %This function induces angular distortion and tries to rectifies it..
% A mix of angular, scaling and shear distortions can be used to simulate
% measurement noise/ outside noise.
theta = 45; %can be changed to test against different angles of rotation. 
S = 0.7; %scaling 
ShearY = 1; %shear distortion
tform = affine2d([S.*cosd(theta) -S.*ShearY*sind(theta) 0;...
    S.*sind(theta) S.*cosd(theta) 0; 0 0 1]);
moving = imwarp(fixed,tform);
moving = moving + uint16(10*rand(size(moving)));
tformEstimate = imregcorr(moving,fixed);
Rfixed = imref2d(size(fixed));
movingReg = imwarp(moving,tformEstimate,'OutputView',Rfixed);
%Reference: https://www.mathworks.com/help/images/use-phase-correlation-as-preprocessing-step-in-registration.html?requestedDomain=in.mathworks.com
end
  
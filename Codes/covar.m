function [feature] = covar(imagec)
I = double(imagec);
[m1 n1] = size(I);
[Ix, Iy] = gradient(I);
[Ixx, Ixy] = gradient(Ix);
[Iyx, Iyy] = gradient(Iy);
Ix = abs(Ix);Iy = abs(Iy);Ixx = abs(Ixx);Iyy = abs(Iyy);
rootsq = sqrt(Ix.^2 + Iy.^2);
for m=1:m1
for n=1:n1
if Iy(m,n) == 0
Iynew(m,n) = 0.0000001; %So that NaN does not result in the arctan below.
else
Iynew(m,n) = Iy(m,n);
end
end
end
tann = Ix./Iynew;
arctan = atan(tann);
fbar = [mean2(I) mean2(Ix) mean2(Iy) mean2(Ixx) mean2(Iyy) mean2(rootsq) mean2(arctan)]';
Covc = zeros(7); 
for m=1:m1
for n=1:n1
final(m,n,:) = [I(m,n) Ix(m,n) Iy(m,n) Ixx(m,n) Iyy(m,n) rootsq(m,n) arctan(m,n)]';
d(:,1) = final(m,n,:);
C  = (d(:,1) - fbar)*(d(:,1) - fbar)'; 
Covc = Covc + C;
end
end
NewCov = (1/((m1*n1)-1))*Covc;
[A B] = eig(NewCov);
[M,yy] = max(B(:));
xx = sqrt(yy);
feature = [A(:,xx)' A(:,xx-1)'];
% Returns a 14 dimensional feature vector.
% Extracting the eigen vectors corresponding to the largest 2 eigen values of the covariance matrix.
% Reference for method: 
% Fatih Porikli and Tekin Kocak. Robust license plate detection using covariance descriptor in a neural network framework. In Proceedings of the IEEE International Conference on
% Video and Signal Based Surveillance
% (AVSS), page 107, 2006.
end

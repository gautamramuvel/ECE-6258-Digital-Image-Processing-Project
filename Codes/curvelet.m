function [avgg] = curvelet(imagex)
%Curvelet transform - Extracting coefficients - Means of the matrices of
%different levels.
sum1 = 0; sum2 = 0; sum3 = 0;
curve = fdct_wrapping(imagex,1, 1, 4, 8);
for ii = 1:8
sum1 = sum1 + mean2(curve{1,2}{1,ii});
end
for ii = 1:16
sum2 = sum2 + mean2(curve{1,3}{1,ii});
sum3 = sum3 + mean2(curve{1,4}{1,ii});
end
avg1 = sum1/8; avg2 = sum2/16; avg3 = sum3/16;
avgg = [avg1 avg2 avg3]; %The final feature vector (3 dimensions)
%Reference: By Laurent Demanet, 2004: http://www.curvelet.org/papers/curvelab.pdf
end

function [optimalscale,r] = findoptimalnumscales( img )
%Scale Selection algorithm:
%Finds the optimal scale to be used with curvelets ensuring centered
%impulse region is not divided.

X = fftshift(fft2(ifftshift(img)))/sqrt(numel(img));
[N1,N2]=size(img);
center_point1=floor(N1/2);
center_point2=floor(N2/2);


X=abs(X);
maximum=max(max(X));
minimum=min(min(X));
midvalue=sqrt(maximum*minimum);


r=3;
while r<(1/2*sqrt((N1-center_point1)^2+(N2-center_point2)^2))

    square=X(center_point1-r:center_point1+r,center_point2-r:center_point2+r);
    
    I=(square>=midvalue);

    if sum(sum(I))< 0.99* numel(I)
        break;
    end
    
    r=r+1;

end

r=r-1;

if r<=4
    
    optimalscale = max(ceil((log2(min(N1,N2)) - 3)),3);
else

coarsest_length=2*r+1;

    optimalscale=log2(min(N1,N2)/coarsest_length);

    optimalscale=max(ceil(optimalscale),4);

end
end


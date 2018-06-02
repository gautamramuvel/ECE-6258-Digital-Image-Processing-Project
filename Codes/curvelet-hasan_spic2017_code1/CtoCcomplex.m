%This function combines the real and imaginary components of curvelet
%coefficients. 
function [Ccomplex]=CtoCcomplex(C,nbangles)
Ccomplex=C;

for i=1:length(C)
    for j=1:length(C{i})/2
        Ccomplex{i}{j} = (C{i}{j}+sqrt(-1)*C{i}{j+sum(nbangles(1,i)+nbangles(2,i))})/sqrt(2);
    end
    
    
    for j=length(C{i})/2+1:length(C{i})
        Ccomplex{i}{j} = real(Ccomplex{i}{j-sum(nbangles(1,i)+nbangles(2,i))})-sqrt(-1)*imag(Ccomplex{i}{j-sum(nbangles(1,i)+nbangles(2,i))});
    end
end

end

function [haarx] = haarf(immage) %Extracting the high frequency components 
haarim = imresize(immage,[64 64]);
for k1 = 1:6
    [a b c dd] = haart2(haarim,k1); 
end
for k1 = 1:6
    haarx(k1) = mean2(dd{k1}); %6 dimensional feature vector
end 
end

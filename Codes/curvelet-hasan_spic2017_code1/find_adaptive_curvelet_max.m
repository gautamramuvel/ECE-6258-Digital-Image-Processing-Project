function [scale_locations,nbangles,PSNR]=find_adaptive_curvelet_max(img,noisy_img,sigma,measure)
%Function to find the optimal adaptive curvelet tiling by maximizing a denoising
%performance measure
% Input arguemnts-
% sigma: noise standard deviation
% measure: denoising performance metric can be
% 'denoiseMSE','denoisePSNR','gini','entropy', 'coeffv'
% gini is the gini index, entropy is shannon's entropy, and coeffv is the
% coefficient of variation.
% numlevels: number of curvelet decomposition levels J.
%
% Output arguments:
% scale_locations: the optimal scale locations found.
% nbangles: optimal angular decomposition.
% PSNR: Optimal value of the metric specified by 'measure'.

numlevels=findoptimalnumscales(img);
options=optimset('MaxIter',500,'TolX',1);
disp('First iteration');

%First iteration:
down4_img=imresize(img,1/4);
down4_noisy_img=imresize(noisy_img,1/4);
[N1,N2]=size(down4_img);
hsize=round(0.10*[N1,N2]);

h = fspecial('gaussian', hsize);
down4_img = imfilter(down4_img,h);
down4_noisy_img = imfilter(down4_noisy_img,h);





%Default scale locations:
scale_locations=zeros(2,numlevels);
scale_locations(:,end)=[4/3*N1;4/3*N2];
scale_locations(:,end-1)=[(4/3*N1-4/6*N1)/2+4/6*N1;(4/3*N2-4/6*N2)/2+4/6*N2];
for i=size(scale_locations,2)-2:-1:1
   scale_locations(:,i)=[(scale_locations(1,i+1)-4/6*N1)/2+4/6*N1;(scale_locations(2,i+1)-4/6*N2)/2+4/6*N2];    
end

scale_locations=round(scale_locations);


[nbangles]=findoptimalangular(measure,down4_noisy_img,down4_img,scale_locations,sigma);
disp('Angular decomposition')
disp(nbangles)

[X] = fminsearch(@(scale_locations) -criteria_denoise(measure,down4_noisy_img,down4_img,scale_locations,nbangles,sigma),scale_locations,options);
disp('Scale locations');
disp(round(X))
disp('Second iteration')
%Second iteration:
scale_locations=floor(2*X);
down2_img=imresize(img,1/2);
down2_noisy_img=imresize(noisy_img,1/2);

[N1,N2]=size(down2_img);
hsize=round(0.05*[N1,N2]);

h = fspecial('gaussian', hsize);
down2_img = imfilter(down2_img,h);
down2_noisy_img = imfilter(down2_noisy_img,h);

[nbangles]=findoptimalangular(measure,down2_noisy_img,down2_img,scale_locations,sigma);
disp('Angular decomposition')
disp(nbangles)

X = fminsearch(@(scale_locations) -criteria_denoise(measure,down2_noisy_img,down2_img,scale_locations,nbangles,sigma),scale_locations,options);

disp('Scale locations');
disp(round(X))


%Last iteration:
disp('Last iteration')

scale_locations=floor(2*X);

[nbangles]=findoptimalangular(measure,noisy_img,img,scale_locations,sigma);
disp('Angular decomposition')
disp(nbangles)



scale_locations = fminsearch(@(scale_locations) -criteria_denoise(measure,noisy_img,img,scale_locations,nbangles,sigma),scale_locations,options);
disp('Scale locations');
disp(round(scale_locations))

[nbangles,PSNR]=findoptimalangular(measure,noisy_img,img,scale_locations,sigma);
disp('Angular decomposition')
disp(nbangles)

end

function output=criteria_denoise(measure,noisy_img,img,scale_locations,nbangles,sigma)

errornbscale=0;
mindist=0;
[N1,N2]=size(noisy_img);
if scale_locations(1,end)<N1+mindist || scale_locations(2,end)<N2+mindist
    errornbscale=1;
end


if scale_locations(1,end-1)>N1 || scale_locations(2,end-1)>N2
    errornbscale=1;
end

for i=1:size(scale_locations,2)-1
    if i==1
        if scale_locations(1,1)-floor(scale_locations(1,end)/2)<mindist
               errornbscale=1;
        end 
    end
    
    if scale_locations(1,i+1)-scale_locations(1,i)<mindist
            errornbscale=1;
        
    end    
    
    if i==1
        if scale_locations(2,1)-floor(scale_locations(2,end)/2)<mindist
               errornbscale=1;
        end 
    end
    
    if scale_locations(2,i+1)-scale_locations(2,i)<mindist
            errornbscale=1;
    end    

end

%Constraint violation check
if errornbscale==0
    output=denoise(measure,noisy_img,img,scale_locations,nbangles,sigma);
else
    output=0;
end

end





function [ finaloutput,PSNR] = findoptimalangular(measure,noisy_barbara,barbara,scale_locations,sigma)
%Testing different frequency divisions:
 
sequence=[4 8 12 16 20 24];
nbangles=4*ones(2,2*size(scale_locations,2));   
nbangles(:,1)=1;
nbangles=[nbangles(:,1) 2*nbangles(:,2:size(scale_locations,2))];
testingnbangles=nbangles;
halfvariables=size(scale_locations,2);
numvariables=2*(size(scale_locations,2)-1);
output=0;

for i=1:numvariables            
       if i<=halfvariables-1              
           for testingangle=sequence
               testingnbangles(1,i+1)=testingangle;
               PSNR=denoise(measure,noisy_barbara,barbara,scale_locations,testingnbangles,sigma);
               curroutput=PSNR;                                                                                                                                                                                                 
            if curroutput>output && isnan(curroutput)==0
              output=curroutput;
               nbangles(1,i+1)=testingangle;
            end
           end           
           testingnbangles=nbangles;
       else            
           for testingangle=sequence
                
               testingnbangles(2,i-(halfvariables)+2)=testingangle; 

               PSNR=denoise(measure,noisy_barbara,barbara,scale_locations,testingnbangles,sigma); 
               curroutput=PSNR;                                                                                                                                                                                                      
            if curroutput>output && isnan(curroutput)==0
              output=curroutput;
               nbangles(2,i-(halfvariables)+2)=testingangle;
            end
           end
           
          testingnbangles=nbangles;

       end 
       
end


finaloutput=nbangles;
PSNR=output;
end


function output=denoise(measure,noisy_img,img,scale_locations,nbangles,sigma)
 switch measure
         case 'coeffv'
            C = adaptive_curvelet(img,1,scale_locations,nbangles);
            y=Ctoy(C);
            y=abs(log(y));
            output=sqrt(std(y)/mean(y));
       
        case 'gini'
            C = adaptive_curvelet(img,1,scale_locations,nbangles);
            y=Ctoy(C);
            y=abs(y);
            output=ginicoeff(y);
        
        case 'entropy'
            C = adaptive_curvelet(img,1,scale_locations,nbangles);
            y=Ctoy(C);
            y=abs(y);
            output=entropy(y);            
            
        case 'denoisePSNR'
            [N1,N2]=size(img);
            F = ones(N1,N2);
            X = fftshift(ifft2(F)) * sqrt(numel(F));
            C = adaptive_curvelet(X,1,scale_locations,nbangles); 

            E = cell(size(C));
            for s=1:length(C)
              E{s} = cell(size(C{s}));
              for w=1:length(C{s})
                A = C{s}{w};
                E{s}{w} = sqrt(sum(sum(A.*conj(A))) / numel(A));
              end
            end

            C = adaptive_curvelet(noisy_img, 1, scale_locations, nbangles);


            Ct = C;
            for s = 2:length(C)
              thresh = 3*sigma + sigma*(s == length(C));
              for w = 1:length(C{s})
                Ct{s}{w} = C{s}{w}.* (abs(C{s}{w}) > thresh*E{s}{w});
              end
            end

            restored_img =  adaptive_icurvelet(Ct, 1, size(img,1),size(img,2),scale_locations, nbangles);
            [N1,N2]=size(restored_img);
            MSE = sum(sum((double(img)-double(restored_img)).^2))/(N1*N2);
            output = 20*log10(255/sqrt(MSE));
            if isnan(output)
                output=0;
            end
            
            
     case 'denoiseMSE'                  
            [N1,N2]=size(img);
            F = ones(N1,N2);
            X = fftshift(ifft2(F)) * sqrt(numel(F));
            C = adaptive_curvelet(X,1,1,scale_locations,nbangles); 
            E = cell(size(C));
            for s=1:length(C)
              E{s} = cell(size(C{s}));
              for w=1:length(C{s})
                A = C{s}{w};
                E{s}{w} = sqrt(sum(sum(A.*conj(A))) / numel(A));
              end
            end

            C = adaptive_curvelet(noisy_img, 1, scale_locations, nbangles);

            
            Ct = C;
            for s = 2:length(C)
              thresh = 3*sigma + sigma*(s == length(C));
              for w = 1:length(C{s})
                Ct{s}{w} = C{s}{w}.* (abs(C{s}{w}) > thresh*E{s}{w});
              end
            end

            restored_img =  adaptive_icurvelet(Ct, 1, size(img,1),size(img,2),scale_locations, nbangles);
            [N1,N2]=size(restored_img);
            output = log10(sum(sum((double(img)-double(restored_img)).^2))/(N1*N2));
           if isnan(output)
               output=0;
           end


     otherwise 
            disp('wrong method');
            return;
 end
    



end


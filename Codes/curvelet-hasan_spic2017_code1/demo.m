% This demo compares default and adaptive curvelet and denoising and partial 
% reconstruction error experiments. Default curvelets are assumed to use
% the 'real' option and periodic extensions.
% 
% [1] Al-Marzouqi, Hasan, and Ghassan AlRegib. "Curvelet transform with learning-based tiling." Signal Processing: Image Communication 53 (2017): 24-39.

%% Denoising example:
img=imread('img.jpg');
img=double(img);
sigma=20;
noisy_img=img+sigma*randn(size(img));
%Finding the adaptive curvelet that maximizes denoising performance as
%measured in PSNR
[scale_locations,nbangles,PSNR]=find_adaptive_curvelet_max(img,noisy_img,sigma,'denoisePSNR');
%Denoising using default and adaptive curvelets
[output_default,output_adapt,restored_img_default,restored_img_optimal]=img_denoise_adaptive_curvelet(img,noisy_img,sigma,'PSNR',scale_locations,nbangles,2);
disp(['Default curvelet denoising: ',num2str(output_default),' PSNR']);
disp(['Adaptive curvelet denoising: ',num2str(output_adapt),' PSNR']);
%% Partial reconstruction error: The image is reconstructed from the highest 3% curvelet coefficients

% A function that combines the real and complex components of
% curvelet coefficients 
C = adaptive_curvelet(img, 1,scale_locations, nbangles);
C=CtoCcomplex(C,nbangles);
[~,b]=size(C);

%Computing the total number of coefficients
numcoeff=0;
for i=1:b
    [~,bb]=size(C{i});
    for ii=1:bb/2
        numcoeff=numcoeff+numel(C{i}{ii});
    end
end

%Taking the magnitude of curvelet coefficients and 
%arranging them in a 1D vector
coeff=zeros(1,numcoeff);
counter=1;
for i=1:b
   [~,bb]=size(C{i});
   for ii=1:bb/2           
      tile=reshape(C{i}{ii},1,numel(C{i}{ii}));
      coeff(counter:counter-1+numel(C{i}{ii}))=abs(tile);
      counter=counter+numel(C{i}{ii});
   end  
end

sortedcoeff_adapt=fliplr(sort(coeff));



nbangles_def =[1     4     8     8    16    16;1     4     8     8    16    16];
D=fdct_wrapping(img,1,1);
% A function that combines the real and complex components of
% curvelet coefficients 
D=CtoCcomplex(D,nbangles_def);
[~,b]=size(D);

%Computing the total number of coefficients
numcoeff=0;
for i=1:b
    [~,bb]=size(D{i});
    for ii=1:bb/2
        numcoeff=numcoeff+numel(D{i}{ii});
    end
end

%Taking the magnitude of curvelet coefficients and 
%arranging them in a 1D vector
coeff=zeros(1,numcoeff);

counter=1;
for i=1:b
   [~,bb]=size(D{i});
   for ii=1:bb/2              
      tile=reshape(D{i}{ii},1,numel(D{i}{ii}));
      coeff(counter:counter-1+numel(D{i}{ii}))=abs(tile);
      counter=counter+numel(D{i}{ii});
   end    
end

sortedcoeff_default=sort(coeff,'descend');


partialC=C;
[~,b]=size(C);

% Zeroing all coefficients:
for i=1:b
   [~,bb]=size(C{i});
   for ii=1:bb           
       partialC{i}{ii}(:,:)=0;       
   end    
end

%Partial reconstruction error is computed at 3% 
rec_per=3;


    
Cpartial=C;
indexth=round(rec_per*numel(img)/100);

Th=sortedcoeff_adapt(indexth);

[~,b]=size(C);

for i=1:b
   [~,bb]=size(C{i});
   for ii=1:bb  
       indicator=abs(C{i}{ii})>=Th;
       Cpartial{i}{ii}=C{i}{ii}.*indicator;
   end
end

%Seperating real and imaginary components to take the inverse curvelet
%transform
Cpartial=CcomplextoC(Cpartial,nbangles);
restored_img_adapt=adaptive_icurvelet(Cpartial,1,size(img,1),size(img,2),scale_locations,nbangles);


MSE = sum(sum((img-restored_img_adapt).^2))/(numel(img));
PSNR_adapt= 20*log10(255/sqrt(MSE));


Dpartial=D;
indexth=round(rec_per*numel(img)/100);


Th=sortedcoeff_default(indexth);
[~,b]=size(D);   
for i=1:b
   [~,bb]=size(D{i});
   for ii=1:bb
       indicator=abs(D{i}{ii})>=Th;
       Dpartial{i}{ii}=D{i}{ii}.*indicator;
   end
end  

%Seperating real and imaginary components to take the inverse curvelet
%transform
Dpartial=CcomplextoC(Dpartial,nbangles_def);

restored_img_def=ifdct_wrapping(Dpartial,1,size(img,1),size(img,2));

MSE = sum(sum((img-restored_img_def).^2))/(numel(img));

PSNR_default = 20*log10(255/sqrt(MSE));       

 figure;  
 imshow(restored_img_def,[]);
 title('Default curvelet image reconstructed from the highest 3% coefficients');
 
 
 figure;  
 imshow(restored_img_adapt,[]);
 title('Adaptive curvelet image reconstructed from the highest 3% coefficients');

 
 

disp('Default curvelet partial reconstruction error: ');
disp(PSNR_default)

disp('Adaptive curvelet partial reconstruction error: ');
disp(PSNR_adapt)


%%


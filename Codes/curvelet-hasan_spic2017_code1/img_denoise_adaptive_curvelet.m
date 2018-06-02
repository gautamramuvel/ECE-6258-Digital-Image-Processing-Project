%Denoising the input image using either default curvelets, adaptive
%curvelets or both.
%if method =0: default result.
%if method=1: adaptive result.
%if method=2: both
function  [output_default,output_optimal,restored_img_default,restored_img_optimal]= img_denoise_adaptive_curvelet(img,noisy_img,sigma,measure,nbscales,nbangles,method)
    
if method==0 ||method==2
    [output_default,restored_img_default]=img_denoise_default_curvelet(img,noisy_img,sigma,measure);   
    output_optimal=0;
    restored_img_optimal=0;
else
    output_default=0;
    restored_img_default=0;
end
    
if method==2 || method==1
    [N1,N2]=size(img);    

    F = ones(N1,N2);
    X = fftshift(ifft2(F)) * sqrt(numel(F));
    C = adaptive_curvelet(X,1,nbscales,nbangles); 
    E = cell(size(C));
    for s=1:length(C)
      E{s} = cell(size(C{s}));
      for w=1:length(C{s})
        A = C{s}{w};
        E{s}{w} = sqrt(sum(sum(A.*conj(A))) / numel(A));
      end
    end

    % Take curvelet transform
    C = adaptive_curvelet(noisy_img, 1, nbscales, nbangles);


    % Apply thresholding
    Ct = C;
    for s = 2:length(C)
      thresh = 3*sigma + sigma*(s == length(C));
      for w = 1:length(C{s})
          Ct{s}{w} = C{s}{w}.* (abs(C{s}{w}) > thresh*E{s}{w});
      end
    end

    
    restored_img_optimal =  adaptive_icurvelet(Ct, 1, size(img,1),size(img,2),nbscales, nbangles);
    

    switch measure
        case 'SSE'
            output_optimal = sum(sum((img-restored_img_optimal).^2));
        case 'MSE'
            output_optimal = sum(sum((img-restored_img_optimal).^2))/numel(img);
        case 'PSNR'
            output_optimal = sum(sum((img-restored_img_optimal).^2))/numel(img);
            output_optimal=20*log10(255/sqrt(output_optimal));        
        case 'SSIM'
            output_optimal=ssim(restored_img_optimal,img);
        otherwise 
            display('Wrong measure. Choose from MSE,PSNR,SSE and SSIM'); 
    end

end






end
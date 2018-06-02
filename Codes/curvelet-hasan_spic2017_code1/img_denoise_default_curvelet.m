function [output,restored_img]=img_denoise_default_curvelet(img,noisy_img,sigma,measure,numlevels)
[N1,N2]=size(img);

F = ones(N1,N2);
X = fftshift(ifft2(F)) * sqrt(numel(F));
if nargin>4
 C = fdct_wrapping(X,1,1,numlevels); 
else
    C = fdct_wrapping(X,1,1); 
end

E = cell(size(C));
for s=1:length(C)
  E{s} = cell(size(C{s}));
  for w=1:length(C{s})
    A = C{s}{w};
    E{s}{w} = sqrt(sum(sum(A.*conj(A))) / numel(A));
  end
end

if nargin>4
 C = fdct_wrapping(noisy_img,1,1,numlevels); 
else
 C = fdct_wrapping(noisy_img,1,1);     
end

% Apply thresholding
Ct = C;
for s = 2:length(C)
  thresh = 3*sigma + sigma*(s == length(C));
  for w = 1:length(C{s})
    Ct{s}{w} = C{s}{w}.* (abs(C{s}{w}) > thresh*E{s}{w});
  end
end


restored_img = real(ifdct_wrapping(Ct,1,N1,N2)); 


switch measure
    case 'SSE'
      output = sum(sum((img-restored_img).^2));
    case 'MSE'
      output = sum(sum((img-restored_img).^2))/numel(img);
    case 'PSNR'
      output = sum(sum((img-restored_img).^2))/numel(img); 
      output=20*log10(255/sqrt(output));
    case 'SSIM'
      output=ssim(restored_img,img);
    otherwise 
            display('Wrong measure. Choose from MSE,PSNR,SSE and SSIM'); 

end




end
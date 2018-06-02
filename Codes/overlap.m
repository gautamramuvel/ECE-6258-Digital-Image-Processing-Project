%Image size = [M,N]
%Overlapping blocks of size [rr,cc]
%Each block is shifted by yy and xx pixels. 

for i = 1:9
  Filename = sprintf('Image_0_0%d.tif', i);
  fullFileName = fullfile('F:\Fall 2017\Digital Image Processing\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class0 - control', Filename);
  control(:,:,i) = (imread(fullFileName));
end   
for i = 10:58
  Filename = sprintf('Image_0_%d.tif', i);
  fullFileName = fullfile('F:\Fall 2017\Digital Image Processing\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class0 - control', Filename);
  control(:,:,i) = (imread(fullFileName));
end

for kk = 1:40
Im = control(:,:,i);
[M,N,~] = size(Im);
rr = 100; cc = 100; xx = 1; yy = 1;
% rr = 2; cc = 2; xx = 1; yy = 1;

numBlocksYY = numel(1:rr-xx:(M-(rr-1)));
numBlocksXX = numel(1:cc-yy:(N-(cc-1)));
[numBlocksYY, numBlocksXX];
C = cell(numBlocksYY*numBlocksXX,1);
counter = 1;
for ii=1:rr-xx:(M-(rr-1))
    for jj=1:cc-yy:(N-(cc-1))
        fprintf('[%d:%d, %d:%d]\n',ii,ii+rr-1,jj,jj+cc-1);
        C{counter} =  Im(ii:(ii+rr-1), jj:(jj+cc-1), : );
        counter = counter + 1;
    end
    fprintf('\n');
end
% figure;
% for ii=1:numBlocksYY*numBlocksXX
%     subplot(numBlocksYY,numBlocksYY,ii), imagesc( C{ii} ); axis image; colormap gray;
% end

for i=1:numBlocksYY*numBlocksXX
   %imwrite(C{i},sprintf('%d.jpg',i))
  Filename = sprintf('Image_%d.tif', (i + ((kk-1)*numBlocksYY*numBlocksXX)));
  fullFileName = fullfile('F:\Fall 2017\Digital Image Processing\DIP Project\Partitioned Images', Filename);
  imwrite(C{i},fullFileName);
end   
end


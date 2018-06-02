function [control, osteo] = getimages(x) %The file names must be changed according to the folder where the train and test images are to be stored..
%imadjust is used to enhance the input images (contrast enhancement)
%imnoise is used to induce Gaussian noise with mean 0 and variance 0.01 (default values in MATLAB)
%weiner2 is used to denoise the image 
%angular() is used to induce rotation and to rectify distortion due to this
% induced rotation.
%%
  if x == 1
        for i = 1:9
  Filename = sprintf('Image_0_0%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class0 - control', Filename);
  control(:,:,i) = imadjust(imread(fullFileName));
          end   
        for i = 10:58
  Filename = sprintf('Image_0_%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class0 - control', Filename);
  control(:,:,i) = imadjust(imread(fullFileName));
        end
% Getting all class 1 data - osteo
        for i = 1:9
  Filename = sprintf('Image_1_0%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class1 - osteo', Filename);
  osteo(:,:,i) = imadjust(imread(fullFileName));
        end   
        for i = 10:58
  Filename = sprintf('Image_1_%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class1 - osteo', Filename);
  osteo(:,:,i) = imadjust(imread(fullFileName));
        end
  end
%%
  if x == 2
        for i = 1:9
  Filename = sprintf('Image_0_0%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class0 - control', Filename);
  control(:,:,i) = (imadjust(imread(fullFileName)));
        end   
        for i = 10:40
  Filename = sprintf('Image_0_%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class0 - control', Filename);
  control(:,:,i) = (imadjust(imread(fullFileName)));
        end
        for i = 41:58
  Filename = sprintf('Image_0_%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class0 - control', Filename);
  control(:,:,i) = imnoise(imadjust(imread(fullFileName)),'gaussian');
  control(:,:,i) = wiener2(control(:,:,i),[5 5]);
        end
% Getting all class 1 data - osteo
        for i = 1:9
  Filename = sprintf('Image_1_0%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class1 - osteo', Filename);
  osteo(:,:,i) = (imadjust(imread(fullFileName)));
        end   
        for i = 10:40
  Filename = sprintf('Image_1_%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class1 - osteo', Filename);
  osteo(:,:,i) = (imadjust(imread(fullFileName)));
        end     
        for i = 41:58
  Filename = sprintf('Image_1_%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class1 - osteo', Filename);
  osteo(:,:,i) = imnoise(imadjust(imread(fullFileName)),'gaussian');
  osteo(:,:,i) = wiener2(osteo(:,:,i),[5 5]);
        end  
  end
%%
  if x == 3
        for i = 1:9
  Filename = sprintf('Image_0_0%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class0 - control', Filename);
  control(:,:,i) = imadjust(imread(fullFileName));
        end   
        for i = 10:58
  Filename = sprintf('Image_0_%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class0 - control', Filename);
  control(:,:,i) = imadjust(imread(fullFileName));
        end
% Getting all class 1 data - osteo
        for i = 1:9
  Filename = sprintf('Image_1_0%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class1 - osteo', Filename);
  osteo(:,:,i) = imadjust(imread(fullFileName));
        end   
        for i = 10:58
  Filename = sprintf('Image_1_%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class1 - osteo', Filename);
  osteo(:,:,i) = imadjust(imread(fullFileName));
        end
        for i = 41:58
  control(:,:,i) = angular(control(:,:,i));
  osteo(:,:,i) = angular(osteo(:,:,i));
        end
  end
%%
  if x == 4
        for i = 1:9
  Filename = sprintf('Image_0_0%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class0 - control', Filename);
  control(:,:,i) = (imadjust(imread(fullFileName)));
        end   
        for i = 10:40
  Filename = sprintf('Image_0_%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class0 - control', Filename);
  control(:,:,i) = (imadjust(imread(fullFileName)));
        end
        for i = 41:58
  Filename = sprintf('Image_0_%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class0 - control', Filename);
  control(:,:,i) = imnoise(imadjust(imread(fullFileName)),'gaussian');
  control(:,:,i) = wiener2(control(:,:,i),[5 5]);
        end
% Getting all class 1 data - osteo
        for i = 1:9
  Filename = sprintf('Image_1_0%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class1 - osteo', Filename);
  osteo(:,:,i) = (imadjust(imread(fullFileName)));
        end   
        for i = 10:40
  Filename = sprintf('Image_1_%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class1 - osteo', Filename);
  osteo(:,:,i) = (imadjust(imread(fullFileName)));
        end
         for i = 41:58
  Filename = sprintf('Image_1_%d.tif', i);
  fullFileName = fullfile('F:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class1 - osteo', Filename);
  osteo(:,:,i) = imnoise(imadjust(imread(fullFileName)),'gaussian');
  osteo(:,:,i) = wiener2(osteo(:,:,i),[5 5]);
        end
        for i = 41:58
  control(:,:,i) = angular(control(:,:,i));
  osteo(:,:,i) = angular(osteo(:,:,i));
        end        
  end
    
end 

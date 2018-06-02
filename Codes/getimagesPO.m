function [control, controltest, osteo, osteotest] = getimagesPO(x) 
%%
%The file names must be changed according to the folder where the train and test images are to be stored..
%imadjust is used to enhance the input images (contrast enhancement)
%imnoise is used to induce Gaussian noise with mean 0 and variance 0.01 (default values in MATLAB)
%weiner2 is used to denoise the image 
%angular() is used to induce rotation and to rectify distortion due to this induced rotation.

%%
  if x == 1
        for i = 1:640
  Filename = sprintf('Image_%d.tif', i);
  fullFileName = fullfile('D:\DIP Project\Partitioned Images\Training\Class 0 Control', Filename);
  control(:,:,i) = imadjust(imread(fullFileName));
        end   
        for i = 641:658
  Filename = sprintf('Image_0_%d.tif', i-600);
  fullFileName = fullfile('D:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class0 - control', Filename);
  controltest(:,:,i) = imadjust(imread(fullFileName));
        end
% Getting all class 1 data - osteo
        for i = 1:640
  Filename = sprintf('Image_%d.tif', i);
  fullFileName = fullfile('D:\DIP Project\Partitioned Images\Training\Class 1 Osteo', Filename);
  osteo(:,:,i) = imadjust(imread(fullFileName));
        end   
        for i = 641:658
  Filename = sprintf('Image_1_%d.tif', i-600);
  fullFileName = fullfile('D:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class1 - osteo', Filename);
  osteotest(:,:,i) = imadjust(imread(fullFileName));
        end
  end
  %%
  if x == 2
        for i = 1:640
  Filename = sprintf('Image_%d.tif', i);
  fullFileName = fullfile('D:\DIP Project\Partitioned Images\Training\Class 0 Control', Filename);
  control(:,:,i) = (imadjust(imread(fullFileName)));
        end   
        for i = 641:658
  Filename = sprintf('Image_0_%d.tif', i-600);
  fullFileName = fullfile('D:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class0 - control', Filename);
  controltest(:,:,i) = imnoise(imadjust(imread(fullFileName)),'gaussian');
  controltest(:,:,i) = wiener2(controltest(:,:,i),[5 5]);
        end
% Getting all class 1 data - osteo
        for i = 1:640
  Filename = sprintf('Image_%d.tif', i);
  fullFileName = fullfile('D:\DIP Project\Partitioned Images\Training\Class 1 Osteo', Filename);
  osteo(:,:,i) = (imadjust(imread(fullFileName)));
        end   
        for i = 641:658
  Filename = sprintf('Image_1_%d.tif', i-600);
  fullFileName = fullfile('D:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class1 - osteo', Filename);
  osteotest(:,:,i) = imnoise(imadjust(imread(fullFileName)),'gaussian');
  osteotest(:,:,i) = wiener2(osteotest(:,:,i),[5 5]);
        end
  end
  %%
  if x == 3
       for i = 1:640
  Filename = sprintf('Image_%d.tif', i);
  fullFileName = fullfile('D:\DIP Project\Partitioned Images\Training\Class 0 Control', Filename);
  control(:,:,i) = imadjust(imread(fullFileName));
        end   
        for i = 641:658
  Filename = sprintf('Image_0_%d.tif', i-600);
  fullFileName = fullfile('D:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class0 - control', Filename);
  controltest(:,:,i) = imadjust(imread(fullFileName));
        end
% Getting all class 1 data - osteo
        for i = 1:640
  Filename = sprintf('Image_%d.tif', i);
  fullFileName = fullfile('D:\DIP Project\Partitioned Images\Training\Class 1 Osteo', Filename);
  osteo(:,:,i) = imadjust(imread(fullFileName));
        end   
        for i = 641:658
  Filename = sprintf('Image_1_%d.tif', i-600);
  fullFileName = fullfile('D:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class1 - osteo', Filename);
  osteotest(:,:,i) = imadjust(imread(fullFileName));
        end   
        for i = 641:658
  controltest(:,:,i) = angular(controltest(:,:,i));
  osteotest(:,:,i) = angular(osteotest(:,:,i));
        end
  end
  %%
  if x == 4
        for i = 1:640
  Filename = sprintf('Image_%d.tif', i);
  fullFileName = fullfile('D:\DIP Project\Partitioned Images\Training\Class 0 Control', Filename);
  control(:,:,i) = (imadjust(imread(fullFileName)));
        end   
        for i = 641:658
  Filename = sprintf('Image_0_%d.tif', i-600);
  fullFileName = fullfile('D:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class0 - control', Filename);
  controltest(:,:,i) = imnoise(imadjust(imread(fullFileName)),'gaussian');
  controltest(:,:,i) = wiener2(controltest(:,:,i),[5 5]);
        end
% Getting all class 1 data - osteo
        for i = 1:640
  Filename = sprintf('Image_%d.tif', i);
  fullFileName = fullfile('D:\DIP Project\Partitioned Images\Training\Class 1 Osteo', Filename);
  osteo(:,:,i) = (imadjust(imread(fullFileName)));
        end   
        for i = 641:658
  Filename = sprintf('Image_1_%d.tif', i-600);
  fullFileName = fullfile('D:\DIP Project\TCB_Challenge_Data\TRAIN_TEST_Data\Class1 - osteo', Filename);
  osteotest(:,:,i) = imnoise(imadjust(imread(fullFileName)),'gaussian');
  osteotest(:,:,i) = wiener2(osteotest(:,:,i),[5 5]);
        end       
        for i = 641:658
    controltest(:,:,i) = angular(controltest(:,:,i));
    osteotest(:,:,i) = angular(osteotest(:,:,i));
        end
  end
    
end 



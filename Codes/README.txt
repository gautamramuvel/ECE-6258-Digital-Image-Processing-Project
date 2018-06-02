*******************APPLICATION OF IMAGE PROCESSING TECHNIQUES TO BONE RADIOGRAPH IMAGES FOR OSTEOPOROSIS DIAGNOSIS*******************
Authors:
1. Krishna Kumar Balakrishnan: krishna.kumar:@gatech.edu (GT ID: 903321701)
2. Anirudha Sundaresan: asundaresan7@gatech.edu (GT ID: 903301401)
(The references for the codes have been included as comments in the .m files)

List of all MATLAB files to be added to the path:
1. finalPO.m - The main code. Just run it after changing address in getimagesPO.m
2. getimagesPO.m - Change the address parameters and then run the finalPO.m. This loads the dataset for training and testing.
3. statistical.m - To extract statisitcal features.
4. covar.m - To extract features from analyzing the covariance of the image matrix.
5. curvelet.m - To extract the curvelet features.
6. haarf.m - To extract the haar wavelet features.
7. cnn.m - To use CNN for our image dataset. (pre-trained model in 'modelseriesnetwork.m')
8. overlap.m - The code used to partition our original dataset.
9. validationada.m/validationsvm.m/validationrf.m - To validate the AdaBoost/SVM/Random Forest classifier used.
10. gabor.m - To extract gabor features from an image.
11. adaboost.m - To run the adaboost training & testing code.
12. angular.m - To introduce and rectify the distortion caused by rotation. 
13. alceufc-sfta-cc83d92 - This folder has the code for fractal dimension feature extraction.
14. curvelet-hasan_spic2017_code1 - This folder has the code required to run curvelet.m
15. Gabor Filter code - This folder contains code for extracting gabor features.
Please make sure that all the folders and sub-folders are all added to the path.
(The codes final.m and getimages.m in 'Extra' folder need not be run. This is for the case without partitioning.)

STEPS:
1. In getimagesPO.m: Change the address in 'fullFileName' accordingly. (There are 16 instances of the 'fullFileName' variable.)
2. Run finalPO.m. You will get all results for all cases, as explained in the paper*. This will take approximately 15 mins.(tested on an i7 processor)
3. CNN has been already trained using cnn.m (please use MATLAB2017b) and the model has been saved as 'modelseriesnetwork.mat'**. You can just comment out the codes not required for training (as specified in cnn.m). 

NOTE: 
*The accuracy results for Random Forest and for classification using Gaussian noise may differ because of the randomness involved. 
**Training the cnn from scratch might give different accuracy, due to randomised weight initialization.

*******************THE END*******************
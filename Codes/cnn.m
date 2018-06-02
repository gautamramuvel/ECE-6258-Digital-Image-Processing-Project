%%Ref: https://www.mathworks.com/help/vision/examples/object-detection-using-deep-learning.html
clear all;
clc;
%Get partitioned images for training and 18 images for testing..
[control, controltest, osteo, osteotest] =  getimagesPO(1);
for i = 1:640
trainingI(:,:,i) = control(:,:,i);
trainingLabel(i) = -1;
end
for i = 1:640
trainingI(:,:,i+640) = osteo(:,:,i);
trainingLabel(i+640) = 1;
end
A = randperm(1280);
trainingImages = trainingI(:,:,A);
trainingLabels = trainingLabel(A);

numImageCategories = 2; %osteoporotic and control subject..
trainingLabels = categorical(trainingLabels);
testImages = controltest(:,:,641:658);
for i = 1:18
testImages(:,:,i+18) = osteotest(:,:,640+i);
end
testImages = reshape(testImages,400,400,1,[]);
testImages = imresize(testImages,[100 100]);
testLabels = [-ones(18,1); ones(18,1)];
testLabels = categorical(testLabels);
trainingImages = reshape(trainingImages,100,100,1,[]);
[height, width, numChannels, ~] = size(trainingImages);
imageSize = [height width numChannels];
inputLayer = imageInputLayer(imageSize)

% Convolutional layer parameters
filterSize = [10 10];
numFilters = 15;
middleLayers = [

% The first convolutional layer has a bank of 32 5x5x3 filters. A
% symmetric padding of 2 pixels is added to ensure that image borders
% are included in the processing. This is important to avoid
% information at the borders being washed away too early in the
% network.
convolution2dLayer(filterSize, numFilters, 'Padding', 2)

% Note that the third dimension of the filter can be omitted because it
% is automatically deduced based on the connectivity of the network. In
% this case because this layer follows the image layer, the third
% dimension must be 3 to match the number of channels in the input
% image.

% Next add the ReLU layer:
reluLayer()

% Follow it with a max pooling layer that has a 3x3 spatial pooling area
% and a stride of 2 pixels. This down-samples the data dimensions from
% 32x32 to 15x15.
maxPooling2dLayer(3, 'Stride', 2)

% Repeat the 3 core layers to complete the middle of the network.
convolution2dLayer(filterSize, numFilters, 'Padding', 2)
reluLayer()
maxPooling2dLayer(3, 'Stride',2)

convolution2dLayer(filterSize, 2 * numFilters, 'Padding', 2)
reluLayer()
maxPooling2dLayer(3, 'Stride',2)

]

finalLayers = [

% Add a fully connected layer with 64 output neurons. The output size of
% this layer will be an array with a length of 64.
fullyConnectedLayer(64)

% Add an ReLU non-linearity.
reluLayer

% Add the last fully connected layer. At this point, the network must
% produce 10 signals that can be used to measure whether the input image
% belongs to one category or another. This measurement is made using the
% subsequent loss layers.
fullyConnectedLayer(numImageCategories)

% Add the softmax loss layer and classification layer. The final layers use
% the output of the fully connected layer to compute the categorical
% probability distribution over the image classes. During the training
% process, all the network weights are tuned to minimize the loss over this
% categorical distribution.
softmaxLayer
classificationLayer
]

layers = [
    inputLayer
    middleLayers
    finalLayers
    ]
%%

%CNN has been already trained and the model has been stored in
%'modelseriesnetwork'. You can use it to get the same accuracy as has been 
%reported in the paper..
%Else, CNN can be trained again by uncommenting the following lines..

layers(2).Weights = 0.0001 * randn([filterSize numChannels numFilters]);

opts = trainingOptions('sgdm', ...
    'Momentum', 0.9, ...
    'InitialLearnRate', 0.000001, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropFactor', 0.1, ...
    'LearnRateDropPeriod', 8, ...
    'L2Regularization', 0.004, ...
    'MaxEpochs', 6, ...
    'MiniBatchSize', 70, ...
    'Verbose', true);
% 
%     % Train a network.
% 
     %model = trainNetwork(trainingImages, trainingLabels, layers, opts);
%%

load modelseriesnetwork %must be comented if CNN is to be trained separately..
model


%% Run the network on the test set.
YTest = classify(model, testImages);

% Calculate the accuracy.
accuracy = sum(YTest == testLabels)/numel(testLabels)

close all; clear all; clc
tic
% Training and testing...640 control and 640 osteoporotic images for
% training, and 18 control and 18 osteoporotic images for testing..
% Training iterations for adaboost: 50
% No. of trees for Random Forest ensemble method: 30
% The SVM kernel used is linear - fitcsvm.
% Validation functions return the ROC curves and accuracy of the classifier.

for x=1:4
    %%
    if x==1        
fprintf('NOISELESS CASE.\n\n')        
tic
clearvars -except x;
warning off; % warning off for SVM classifier.
[control, controltest, osteo, osteotest] =  getimagesPO(x);
gaborArray = gaborFilterBank(2,2,39,39); % for general filterbank

% Getting features..
 for i = 1:640 
  [ccon(i),ccor(i),cene(i),chom(i),cmean(i),cstd(i)] = statistical(control(:,:,i));
  [ocon(i),ocor(i),oene(i),ohom(i),omean(i),ostd(i)] = statistical(osteo(:,:,i));
  featureVectorc(:,:,i) = gabor(control(:,:,i),gaborArray);
  featureVectoro(:,:,i) = gabor(osteo(:,:,i),gaborArray);
  [avggc(:,:,i)] = curvelet(control(:,:,i));
  [avggo(:,:,i)] = curvelet(osteo(:,:,i));
  featureC(:,:,i) = covar(control(:,:,i));
  featureO(:,:,i) = covar(osteo(:,:,i));
  chaar(i,:) = haarf(control(:,:,i));  
  ohaar(i,:) = haarf(osteo(:,:,i));
 end
 
 for i = 641:658 
  [ccon(i),ccor(i),cene(i),chom(i),cmean(i),cstd(i)] = statistical(controltest(:,:,i));
  [ocon(i),ocor(i),oene(i),ohom(i),omean(i),ostd(i)] = statistical(osteotest(:,:,i));
  featureVectorc(:,:,i) = gabor(controltest(:,:,i),gaborArray);
  featureVectoro(:,:,i) = gabor(osteotest(:,:,i),gaborArray);
  [avggc(:,:,i)] = curvelet(controltest(:,:,i));
  [avggo(:,:,i)] = curvelet(osteotest(:,:,i));
  featureC(:,:,i) = covar(controltest(:,:,i));
  featureO(:,:,i) = covar(osteotest(:,:,i));
  chaar(i,:) = haarf(controltest(:,:,i));  
  ohaar(i,:) = haarf(osteotest(:,:,i));
 end
 
fprintf('Features have been extracted.\n')
% Training and testing...
Y = [zeros(1,640), ones(1,640)]; %Y - labels for training data for SVM classifier. Yada = [-ones(1,640), ones(1,640)];
Yada = [-ones(1,640), ones(1,640)]; %Yada - labels for training data for Adaboost/Random Forest classifier.
Y1blind = [ones(1,18) zeros(1,18)]; %labels for the blind data for SVM classifer
Y1blindada = [ones(1,18) -ones(1,18)]; %labels for the blind data for Adaboost and Random Forest classifier..
%% ALL FEATURES
for i = 1:640
    TRAINING(i,:) = [extractLBPFeatures(control(:,:,i)) featureVectorc(:,:,i) sfta(control(:,:,i),1) chaar(i,:) avggc(:,:,i) min(min(dct2(double(control(:,:,i))))) ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(control(:,:,i)) featureC(:,:,i)]; %control
    TRAINING(i+640,:) = [extractLBPFeatures(osteo(:,:,i)) featureVectoro(:,:,i) sfta(osteo(:,:,i),1) ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteo(:,:,i))))) ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteo(:,:,i)) featureO(:,:,i)];
end
[estimateclasstotal,model]=adaboost('train',TRAINING,Yada,50);
TreeObject=TreeBagger(30,TRAINING,Yada,'method','classification','NVarToSample','all');
SVMStruct = fitcsvm(TRAINING, Y);%,'kernel_function','polynomial','polyorder',2); 
disp('Training completed for all features');
for i = 641:658 
    testingblind(i - 640,:) = [extractLBPFeatures(osteotest(:,:,i)) featureVectoro(:,:,i) sfta(osteotest(:,:,i),1) ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteotest(:,:,i))))) ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteotest(:,:,i)) featureO(:,:,i)]; 
    testingblind(18 + i - 640,:) = [extractLBPFeatures(controltest(:,:,i)) featureVectorc(:,:,i) sfta(controltest(:,:,i),1) chaar(i,:) avggc(:,:,i) min(min(dct2(double(controltest(:,:,i))))) ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(controltest(:,:,i)) featureC(:,:,i)]; %control
end
for i = 1:36 
 resultblind(i) = predict(SVMStruct, testingblind(i,:));
 resultadaboostblind(i)=adaboost('apply',testingblind(i,:),model);
 resultblindtree(i) = predict(TreeObject,testingblind(i,:));
end
resultblindtree = str2double(resultblindtree);
disp('Blind SVM results for all features');
validationsvm(Y1blind,resultblind,x);
disp('Blind Adaboost results for all features');
validationada(Y1blindada,resultadaboostblind,x);
disp('Blind RF results for all features');
validationrf(Y1blindada,resultblindtree,x);

%% TRANSFORM BASED FEATURES
Y = [zeros(1,640), ones(1,640)]; %Y - labels for training data for SVM classifier. Yada = [-ones(1,640), ones(1,640)];
Yada = [-ones(1,640), ones(1,640)]; %Yada - labels for training data for Adaboost/Random Forest classifier.
Y1blind = [ones(1,18) zeros(1,18)]; %labels for the blind data for SVM classifer
Y1blindada = [ones(1,18) -ones(1,18)]; %labels for the blind data for Adaboost and Random Forest classifier..
clear TRAINING, clear testingblind, clear resultblindtree;

for i = 1:640
    TRAINING(i,:) = [chaar(i,:) avggc(:,:,i) min(min(dct2(double(control(:,:,i)))))]; %control
    TRAINING(i+640,:) = [ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteo(:,:,i)))))];
end

[estimateclasstotal,model]=adaboost('train',TRAINING,Yada,50);
TreeObject=TreeBagger(30,TRAINING,Yada,'method','classification','NVarToSample','all');
SVMStruct = fitcsvm(TRAINING, Y);%,'kernel_function','polynomial','polyorder',2); 
disp('Training completed - Transform based features only');

for i = 641:658 
    testingblind(i - 640,:) = [ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteotest(:,:,i)))))];  
    testingblind(18 + i - 640,:) = [chaar(i,:) avggc(:,:,i) min(min(dct2(double(controltest(:,:,i)))))];
end

for i = 1:36 
 resultblind(i) = predict(SVMStruct, testingblind(i,:));
 resultadaboostblind(i)=adaboost('apply',testingblind(i,:),model);
 resultblindtree(i) = predict(TreeObject,testingblind(i,:));
end

resultblindtree = str2double(resultblindtree);
disp('Blind SVM results for transform based features');
validationsvm(Y1blind,resultblind,x);
disp('Blind Adaboost results for transform based features');
validationada(Y1blindada,resultadaboostblind,x);
disp('Blind RF results for transform based features');
validationrf(Y1blindada,resultblindtree,x);

%% SPATIAL FEATURES

Y = [zeros(1,640), ones(1,640)]; %Y - labels for training data for SVM classifier. Yada = [-ones(1,640), ones(1,640)];
Yada = [-ones(1,640), ones(1,640)]; %Yada - labels for training data for Adaboost/Random Forest classifier.
Y1blind = [ones(1,18) zeros(1,18)]; %labels for the blind data for SVM classifer
Y1blindada = [ones(1,18) -ones(1,18)]; %labels for the blind data for Adaboost and Random Forest classifier..
clear TRAINING, clear testingblind, clear resultblindtree;
for i = 1:640
    TRAINING(i,:) = [ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(control(:,:,i)) sfta(control(:,:,i),1) featureC(:,:,i) featureVectorc(:,:,i) extractLBPFeatures(control(:,:,i))]; %control
    TRAINING(i+640,:) = [ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteo(:,:,i)) sfta(osteo(:,:,i),1) featureO(:,:,i) featureVectoro(:,:,i) extractLBPFeatures(osteo(:,:,i))];
end

[estimateclasstotal,model]=adaboost('train',TRAINING,Yada,50);
TreeObject=TreeBagger(30,TRAINING,Yada,'method','classification','NVarToSample','all');
SVMStruct = fitcsvm(TRAINING, Y);%,'kernel_function','polynomial','polyorder',2); 
disp('Training completed - Spatial features only');

for i = 641:658 
    testingblind(i - 640,:) = [ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteotest(:,:,i)) sfta(osteotest(:,:,i),1) featureO(:,:,i) featureVectoro(:,:,i) extractLBPFeatures(osteotest(:,:,i))];  
    testingblind(18 + i - 640,:) = [ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(controltest(:,:,i)) sfta(controltest(:,:,i),1) featureC(:,:,i) featureVectorc(:,:,i) extractLBPFeatures(controltest(:,:,i))];
end

for i = 1:36 
 resultblind(i) = predict(SVMStruct, testingblind(i,:));
 resultadaboostblind(i)=adaboost('apply',testingblind(i,:),model);
 resultblindtree(i) = predict(TreeObject,testingblind(i,:));
end

resultblindtree = str2double(resultblindtree);
disp('Blind SVM results for spatial features');
validationsvm(Y1blind,resultblind,x);
disp('Blind Adaboost results for spatial features');
validationada(Y1blindada,resultadaboostblind,x);
disp('Blind RF results for spatial features');
validationrf(Y1blindada,resultblindtree,x);

%% DIRECTIONAL TRANSFORM BASED FEATURES
Y = [zeros(1,640), ones(1,640)]; %Y - labels for training data for SVM classifier. Yada = [-ones(1,640), ones(1,640)];
Yada = [-ones(1,640), ones(1,640)]; %Yada - labels for training data for Adaboost/Random Forest classifier.
Y1blind = [ones(1,18) zeros(1,18)]; %labels for the blind data for SVM classifer
Y1blindada = [ones(1,18) -ones(1,18)]; %labels for the blind data for Adaboost and Random Forest classifier..
clear TRAINING, clear testingblind, clear resultblindtree;

for i = 1:640
    TRAINING(i,:) = [avggc(:,:,i)]; %control
    TRAINING(i+640,:) = [avggo(:,:,i)];
end

[estimateclasstotal,model]=adaboost('train',TRAINING,Yada,50);
TreeObject=TreeBagger(30,TRAINING,Yada,'method','classification','NVarToSample','all');
SVMStruct = fitcsvm(TRAINING, Y);%,'kernel_function','polynomial','polyorder',2); 
disp('Training completed - Directional transform based features only');

for i = 641:658 
    testingblind(i - 640,:) = [avggo(:,:,i)];  
    testingblind(18 + i - 640,:) = [avggc(:,:,i)];
end

for i = 1:36 
 resultblind(i) = predict(SVMStruct, testingblind(i,:));
 resultadaboostblind(i)=adaboost('apply',testingblind(i,:),model);
 resultblindtree(i) = predict(TreeObject,testingblind(i,:));
end

resultblindtree = str2double(resultblindtree);
disp('Blind SVM results for directional (frequency domain) features');
validationsvm(Y1blind,resultblind,x);
disp('Blind Adaboost results for directional (frequency domain) features');
validationada(Y1blindada,resultadaboostblind,x);
disp('Blind RF results for directional (frequency domain) features');
validationrf(Y1blindada,resultblindtree,x);


%% NON DIRECTIONAL TRANSFORM BASED FEATURES
Y = [zeros(1,640), ones(1,640)]; %Y - labels for training data for SVM classifier. Yada = [-ones(1,640), ones(1,640)];
Yada = [-ones(1,640), ones(1,640)]; %Yada - labels for training data for Adaboost/Random Forest classifier.
Y1blind = [ones(1,18) zeros(1,18)]; %labels for the blind data for SVM classifer
Y1blindada = [ones(1,18) -ones(1,18)]; %labels for the blind data for Adaboost and Random Forest classifier..
clear TRAINING, clear testingblind, clear resultblindtree;

for i = 1:640
    TRAINING(i,:) = [chaar(i,:) min(min(dct2(double(control(:,:,i)))))]; %control
    TRAINING(i+640,:) = [ohaar(i,:) min(min(dct2(double(osteo(:,:,i)))))];
end

[estimateclasstotal,model]=adaboost('train',TRAINING,Yada,50);
TreeObject=TreeBagger(30,TRAINING,Yada,'method','classification','NVarToSample','all');
SVMStruct = fitcsvm(TRAINING, Y);%,'kernel_function','polynomial','polyorder',2); 
disp('Training completed - Non-directional transform based features only');

for i = 641:658 
    testingblind(i - 640,:) = [ohaar(i,:) min(min(dct2(double(osteotest(:,:,i)))))];  
    testingblind(18 + i - 640,:) = [chaar(i,:) min(min(dct2(double(controltest(:,:,i)))))];
end

for i = 1:36 
 resultblind(i) = predict(SVMStruct, testingblind(i,:));
 resultadaboostblind(i)=adaboost('apply',testingblind(i,:),model);
 resultblindtree(i) = predict(TreeObject,testingblind(i,:));
end

resultblindtree = str2double(resultblindtree);
disp('Blind SVM results for non-directional (frequency domain) features');
validationsvm(Y1blind,resultblind,x);
disp('Blind Adaboost results for non-directional (frequency domain) features');
validationada(Y1blindada,resultadaboostblind,x);
disp('Blind RF results for non-directional (frequency domain) features');
validationrf(Y1blindada,resultblindtree,x);

beep on; beep;
toc      
    end
    
    %%
    if x==2
        fprintf('\n\n\n') 
fprintf('NOISY - GAUSSIAN.\n\n') 
tic
clearvars -except x;;
warning off;  % warning off for SVM classifier.
[control, controltest, osteo, osteotest] =  getimagesPO(x);
gaborArray = gaborFilterBank(2,2,39,39); % for general filterbank

% Getting features..
 for i = 1:640 %6658
  [ccon(i),ccor(i),cene(i),chom(i),cmean(i),cstd(i)] = statistical(control(:,:,i));
  [ocon(i),ocor(i),oene(i),ohom(i),omean(i),ostd(i)] = statistical(osteo(:,:,i));
  featureVectorc(:,:,i) = gabor(control(:,:,i),gaborArray);
  featureVectoro(:,:,i) = gabor(osteo(:,:,i),gaborArray);
  [avggc(:,:,i)] = curvelet(control(:,:,i));
  [avggo(:,:,i)] = curvelet(osteo(:,:,i));
  featureC(:,:,i) = covar(control(:,:,i));
  featureO(:,:,i) = covar(osteo(:,:,i));
  chaar(i,:) = haarf(control(:,:,i));  
  ohaar(i,:) = haarf(osteo(:,:,i));
 end

  for i = 641:658 
  [ccon(i),ccor(i),cene(i),chom(i),cmean(i),cstd(i)] = statistical(controltest(:,:,i));
  [ocon(i),ocor(i),oene(i),ohom(i),omean(i),ostd(i)] = statistical(osteotest(:,:,i));
  featureVectorc(:,:,i) = gabor(controltest(:,:,i),gaborArray);
  featureVectoro(:,:,i) = gabor(osteotest(:,:,i),gaborArray);
  [avggc(:,:,i)] = curvelet(controltest(:,:,i));
  [avggo(:,:,i)] = curvelet(osteotest(:,:,i));
  featureC(:,:,i) = covar(controltest(:,:,i));
  featureO(:,:,i) = covar(osteotest(:,:,i));
  chaar(i,:) = haarf(controltest(:,:,i));  
  ohaar(i,:) = haarf(osteotest(:,:,i));
  end
 
fprintf('Features have been extracted.\n')
% Training and testing...
Y = [zeros(1,640), ones(1,640)]; %Y - labels for training data for SVM classifier. Yada = [-ones(1,640), ones(1,640)];
Yada = [-ones(1,640), ones(1,640)]; %Yada - labels for training data for Adaboost/Random Forest classifier.
Y1blind = [ones(1,18) zeros(1,18)]; %labels for the blind data for SVM classifer
Y1blindada = [ones(1,18) -ones(1,18)]; %labels for the blind data for Adaboost and Random Forest classifier..
%% ALL FEATURES
for i = 1:640
    TRAINING(i,:) = [extractLBPFeatures(control(:,:,i)) featureVectorc(:,:,i) sfta(control(:,:,i),1) chaar(i,:) avggc(:,:,i) min(min(dct2(double(control(:,:,i))))) ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(control(:,:,i)) featureC(:,:,i)]; %control
    TRAINING(i+640,:) = [extractLBPFeatures(osteo(:,:,i)) featureVectoro(:,:,i) sfta(osteo(:,:,i),1) ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteo(:,:,i))))) ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteo(:,:,i)) featureO(:,:,i)];
end

[estimateclasstotal,model]=adaboost('train',TRAINING,Yada,50);
TreeObject=TreeBagger(30,TRAINING,Yada,'method','classification','NVarToSample','all');
SVMStruct = fitcsvm(TRAINING, Y);%'kernel_function','polynomial','polyorder',2); 
disp('Training completed for all features');

for i = 641:658 
    testingblind(i - 640,:) = [extractLBPFeatures(osteotest(:,:,i)) featureVectoro(:,:,i) sfta(osteotest(:,:,i),1) ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteotest(:,:,i))))) ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteotest(:,:,i)) featureO(:,:,i)]; 
    testingblind(18 + i - 640,:) = [extractLBPFeatures(controltest(:,:,i)) featureVectorc(:,:,i) sfta(controltest(:,:,i),1) chaar(i,:) avggc(:,:,i) min(min(dct2(double(controltest(:,:,i))))) ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(controltest(:,:,i)) featureC(:,:,i)]; %control
end

for i = 1:36 
 resultblind(i) = predict(SVMStruct, testingblind(i,:));
 resultadaboostblind(i)=adaboost('apply',testingblind(i,:),model);
 resultblindtree(i) = predict(TreeObject,testingblind(i,:));
end

resultblindtree = str2double(resultblindtree);
disp('Blind SVM results for all features');
validationsvm(Y1blind,resultblind,x);
disp('Blind Adaboost results for all features');
validationada(Y1blindada,resultadaboostblind,x);
disp('Blind RF results for all features');
validationrf(Y1blindada,resultblindtree,x);

beep on; beep;
toc
        
    end
    
    if x==3
%%
fprintf('\n\n\n') 
fprintf('WITH ANGULAR DISTORTION.\n\n')   
tic
clearvars -except x;;
warning off; ; % warning off for SVM classifier.
[control, controltest, osteo, osteotest] =  getimagesPO(x);
gaborArray = gaborFilterBank(2,2,39,39); % for general filterbank

% Getting features..
 for i = 1:640 %6658
  [ccon(i),ccor(i),cene(i),chom(i),cmean(i),cstd(i)] = statistical(control(:,:,i));
  [ocon(i),ocor(i),oene(i),ohom(i),omean(i),ostd(i)] = statistical(osteo(:,:,i));
  featureVectorc(:,:,i) = gabor(control(:,:,i),gaborArray);
  featureVectoro(:,:,i) = gabor(osteo(:,:,i),gaborArray);
  [avggc(:,:,i)] = curvelet(control(:,:,i));
  [avggo(:,:,i)] = curvelet(osteo(:,:,i));
  featureC(:,:,i) = covar(control(:,:,i));
  featureO(:,:,i) = covar(osteo(:,:,i));
  chaar(i,:) = haarf(control(:,:,i));  
  ohaar(i,:) = haarf(osteo(:,:,i));
 end

  for i = 641:658 
  [ccon(i),ccor(i),cene(i),chom(i),cmean(i),cstd(i)] = statistical(controltest(:,:,i));
  [ocon(i),ocor(i),oene(i),ohom(i),omean(i),ostd(i)] = statistical(osteotest(:,:,i));
  featureVectorc(:,:,i) = gabor(controltest(:,:,i),gaborArray);
  featureVectoro(:,:,i) = gabor(osteotest(:,:,i),gaborArray);
  [avggc(:,:,i)] = curvelet(controltest(:,:,i));
  [avggo(:,:,i)] = curvelet(osteotest(:,:,i));
  featureC(:,:,i) = covar(controltest(:,:,i));
  featureO(:,:,i) = covar(osteotest(:,:,i));
  chaar(i,:) = haarf(controltest(:,:,i));  
  ohaar(i,:) = haarf(osteotest(:,:,i));
  end
 
fprintf('Features have been extracted.\n')
% Training and testing...
Y = [zeros(1,640), ones(1,640)]; %Y - labels for training data for SVM classifier. Yada = [-ones(1,640), ones(1,640)];
Yada = [-ones(1,640), ones(1,640)]; %Yada - labels for training data for Adaboost/Random Forest classifier.
Y1blind = [ones(1,18) zeros(1,18)]; %labels for the blind data for SVM classifer
Y1blindada = [ones(1,18) -ones(1,18)]; %labels for the blind data for Adaboost and Random Forest classifier..
%% ALL FEATURES
for i = 1:640
    TRAINING(i,:) = [extractLBPFeatures(control(:,:,i)) featureVectorc(:,:,i) sfta(control(:,:,i),1) chaar(i,:) avggc(:,:,i) min(min(dct2(double(control(:,:,i))))) ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(control(:,:,i)) featureC(:,:,i)]; %control
    TRAINING(i+640,:) = [extractLBPFeatures(osteo(:,:,i)) featureVectoro(:,:,i) sfta(osteo(:,:,i),1) ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteo(:,:,i))))) ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteo(:,:,i)) featureO(:,:,i)];
end

[estimateclasstotal,model]=adaboost('train',TRAINING,Yada,50);
TreeObject=TreeBagger(30,TRAINING,Yada,'method','classification','NVarToSample','all');
SVMStruct = fitcsvm(TRAINING, Y);%'kernel_function','polynomial','polyorder',2); 
disp('Training completed for all features');

for i = 641:658 
    testingblind(i - 640,:) = [extractLBPFeatures(osteotest(:,:,i)) featureVectoro(:,:,i) sfta(osteotest(:,:,i),1) ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteotest(:,:,i))))) ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteotest(:,:,i)) featureO(:,:,i)]; 
    testingblind(18 + i - 640,:) = [extractLBPFeatures(controltest(:,:,i)) featureVectorc(:,:,i) sfta(controltest(:,:,i),1) chaar(i,:) avggc(:,:,i) min(min(dct2(double(controltest(:,:,i))))) ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(controltest(:,:,i)) featureC(:,:,i)]; %control
end

for i = 1:36 
 resultblind(i) = predict(SVMStruct, testingblind(i,:));
 resultadaboostblind(i)=adaboost('apply',testingblind(i,:),model);
 resultblindtree(i) = predict(TreeObject,testingblind(i,:));
end

resultblindtree = str2double(resultblindtree);
disp('Blind SVM results for all features');
validationsvm(Y1blind,resultblind,x);
disp('Blind Adaboost results for all features');
validationada(Y1blindada,resultadaboostblind,x);
disp('Blind RF results for all features');
validationrf(Y1blindada,resultblindtree,x);
beep on; beep;
toc

    end
    
    if x==4
%%
fprintf('\n\n\n') 
fprintf('WITH NOISE AND ANGULAR DISTORTION.\n\n')   
tic
clearvars -except x;;
warning off; % warning off for SVM classifier.
[control, controltest, osteo, osteotest] =  getimagesPO(x);
gaborArray = gaborFilterBank(2,2,39,39); % for general filterbank

% Getting features..
 for i = 1:640 %6658
  [ccon(i),ccor(i),cene(i),chom(i),cmean(i),cstd(i)] = statistical(control(:,:,i));
  [ocon(i),ocor(i),oene(i),ohom(i),omean(i),ostd(i)] = statistical(osteo(:,:,i));
  featureVectorc(:,:,i) = gabor(control(:,:,i),gaborArray);
  featureVectoro(:,:,i) = gabor(osteo(:,:,i),gaborArray);
  [avggc(:,:,i)] = curvelet(control(:,:,i));
  [avggo(:,:,i)] = curvelet(osteo(:,:,i));
  featureC(:,:,i) = covar(control(:,:,i));
  featureO(:,:,i) = covar(osteo(:,:,i));
  chaar(i,:) = haarf(control(:,:,i));  
  ohaar(i,:) = haarf(osteo(:,:,i));
 end
 for i = 641:658 
  [ccon(i),ccor(i),cene(i),chom(i),cmean(i),cstd(i)] = statistical(controltest(:,:,i));
  [ocon(i),ocor(i),oene(i),ohom(i),omean(i),ostd(i)] = statistical(osteotest(:,:,i));
  featureVectorc(:,:,i) = gabor(controltest(:,:,i),gaborArray);
  featureVectoro(:,:,i) = gabor(osteotest(:,:,i),gaborArray);
  [avggc(:,:,i)] = curvelet(controltest(:,:,i));
  [avggo(:,:,i)] = curvelet(osteotest(:,:,i));
  featureC(:,:,i) = covar(controltest(:,:,i));
  featureO(:,:,i) = covar(osteotest(:,:,i));
  chaar(i,:) = haarf(controltest(:,:,i));  
  ohaar(i,:) = haarf(osteotest(:,:,i));
 end
 
fprintf('Features have been extracted.\n')
% Training and testing...
Y = [zeros(1,640), ones(1,640)]; %Y - labels for training data for SVM classifier. Yada = [-ones(1,640), ones(1,640)];
Yada = [-ones(1,640), ones(1,640)]; %Yada - labels for training data for Adaboost/Random Forest classifier.
Y1blind = [ones(1,18) zeros(1,18)]; %labels for the blind data for SVM classifer
Y1blindada = [ones(1,18) -ones(1,18)]; %labels for the blind data for Adaboost and Random Forest classifier..
%% ALL FEATURES
for i = 1:640
    TRAINING(i,:) = [extractLBPFeatures(control(:,:,i)) featureVectorc(:,:,i) sfta(control(:,:,i),1) chaar(i,:) avggc(:,:,i) min(min(dct2(double(control(:,:,i))))) ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(control(:,:,i)) featureC(:,:,i)]; %control
    TRAINING(i+640,:) = [extractLBPFeatures(osteo(:,:,i)) featureVectoro(:,:,i) sfta(osteo(:,:,i),1) ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteo(:,:,i))))) ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteo(:,:,i)) featureO(:,:,i)];
end

[estimateclasstotal,model]=adaboost('train',TRAINING,Yada,50);
TreeObject=TreeBagger(30,TRAINING,Yada,'method','classification','NVarToSample','all');
SVMStruct = fitcsvm(TRAINING, Y);%,'kernel_function','polynomial','polyorder',2); 
disp('Training completed for all features');

for i = 641:658 
    testingblind(i - 640,:) = [extractLBPFeatures(osteotest(:,:,i)) featureVectoro(:,:,i) sfta(osteotest(:,:,i),1) ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteotest(:,:,i))))) ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteotest(:,:,i)) featureO(:,:,i)]; 
    testingblind(18 + i - 640,:) = [extractLBPFeatures(controltest(:,:,i)) featureVectorc(:,:,i) sfta(controltest(:,:,i),1) chaar(i,:) avggc(:,:,i) min(min(dct2(double(controltest(:,:,i))))) ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(controltest(:,:,i)) featureC(:,:,i)]; %control
end

for i = 1:36 
 resultblind(i) = predict(SVMStruct, testingblind(i,:));
 resultadaboostblind(i)=adaboost('apply',testingblind(i,:),model);
 resultblindtree(i) = predict(TreeObject,testingblind(i,:));
end

resultblindtree = str2double(resultblindtree);
disp('Blind SVM results for all features');
validationsvm(Y1blind,resultblind,x);
disp('Blind Adaboost results for all features');
validationada(Y1blindada,resultadaboostblind,x);
disp('Blind RF results for all features');
validationrf(Y1blindada,resultblindtree,x);
beep on; beep;
toc
    end
end
toc

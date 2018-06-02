close all; clear all; clc
tic
% Training and testing...40 control and 40 osteoporotic images for
% training, and 18 control and 18 osteoporotic images for testing..
% Training iterations for adaboost: 50
% No. of trees for Random Forest ensemble method: 30
% A linear SVM is used. - fitcsvm.
% Validation functions return the ROC curves and accuracy of the classifier.

for x = 1:4
    
  if x == 1      
  fprintf('NOISELESS CASE.\n\n')    
  tic
  clearvars -except x;
  [control, osteo] =  getimages(x); %Refer readme file to set path for getimages()
  gaborArray = gaborFilterBank(2,2,39,39); % for general filterbank
% Getting features..
 for i = 1:58 
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
Y = [zeros(1,40), ones(1,40)]; %Y - labels for training data for SVM classifier.
Yada = [-ones(1,40), ones(1,40)]; %Yada - labels for training data for Adaboost/Random Forest classifier.
Y1blind = [ones(1,18) zeros(1,18)]; %labels for the blind data for SVM classifer
Y1blindada = [ones(1,18) -ones(1,18)]; %labels for the blind data for Adaboost and Random Forest classifier..
fprintf('Features have been extracted.\n')
%% ALL FEATURES
for i = 1:40
    TRAINING(i,:) = [extractLBPFeatures(control(:,:,i)) featureVectorc(:,:,i) sfta(control(:,:,i),1) chaar(i,:) avggc(:,:,i) min(min(dct2(double(control(:,:,i))))) ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(control(:,:,i)) featureC(:,:,i)]; %control
    TRAINING(i+40,:) = [extractLBPFeatures(osteo(:,:,i)) featureVectoro(:,:,i) sfta(osteo(:,:,i),1) ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteo(:,:,i))))) ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteo(:,:,i)) featureO(:,:,i)];
end

[estimateclasstotal,model]=adaboost('train',TRAINING,Yada,50); 
TreeObject=TreeBagger(30,TRAINING,Yada,'method','classification','NVarToSample','all'); 
SVMStruct = fitcsvm(TRAINING, Y);
disp('Training completed for all features.');

for i = 41:58 
    testingblind(i - 40,:) = [extractLBPFeatures(osteo(:,:,i)) featureVectoro(:,:,i) sfta(osteo(:,:,i),1) ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteo(:,:,i))))) ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteo(:,:,i)) featureO(:,:,i)]; 
    testingblind(18 + i - 40,:) = [extractLBPFeatures(control(:,:,i)) featureVectorc(:,:,i) sfta(control(:,:,i),1) chaar(i,:) avggc(:,:,i) min(min(dct2(double(control(:,:,i))))) ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(control(:,:,i)) featureC(:,:,i)]; %control
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
Y = [zeros(1,40), ones(1,40)]; %Y - labels for training data for SVM classifier.
Yada = [-ones(1,40), ones(1,40)]; %Yada - labels for training data for Adaboost/Random Forest classifier.
Y1blind = [ones(1,18) zeros(1,18)]; %labels for the blind data for SVM classifer
Y1blindada = [ones(1,18) -ones(1,18)]; %labels for the blind data for Adaboost and Random Forest classifier..
clear TRAINING, clear testingblind, clear resultblindtree;
for i = 1:40
    TRAINING(i,:) = [chaar(i,:) avggc(:,:,i) min(min(dct2(double(control(:,:,i)))))]; %control
    TRAINING(i+40,:) = [ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteo(:,:,i)))))];
end

[estimateclasstotal,model]=adaboost('train',TRAINING,Yada,50);
TreeObject=TreeBagger(30,TRAINING,Yada,'method','classification','NVarToSample','all');
SVMStruct = fitcsvm(TRAINING, Y);
disp('Training completed - Transform based features only');

for i = 41:58 
    testingblind(i - 40,:) = [ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteo(:,:,i)))))];  
    testingblind(18 + i - 40,:) = [chaar(i,:) avggc(:,:,i) min(min(dct2(double(control(:,:,i)))))];
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
Y = [zeros(1,40), ones(1,40)]; %Y - labels for training data for SVM classifier.
Yada = [-ones(1,40), ones(1,40)]; %Yada - labels for training data for Adaboost/Random Forest classifier.
Y1blind = [ones(1,18) zeros(1,18)]; %labels for the blind data for SVM classifer
Y1blindada = [ones(1,18) -ones(1,18)]; %labels for the blind data for Adaboost and Random Forest classifier..
clear TRAINING, clear testingblind, clear resultblindtree;
for i = 1:40
    TRAINING(i,:) = [ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(control(:,:,i)) sfta(control(:,:,i),1) featureC(:,:,i) featureVectorc(:,:,i) extractLBPFeatures(control(:,:,i))]; %control
    TRAINING(i+40,:) = [ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteo(:,:,i)) sfta(osteo(:,:,i),1) featureO(:,:,i) featureVectoro(:,:,i) extractLBPFeatures(osteo(:,:,i))];
end

[estimateclasstotal,model]=adaboost('train',TRAINING,Yada,50);
TreeObject=TreeBagger(30,TRAINING,Yada,'method','classification','NVarToSample','all');
SVMStruct = fitcsvm(TRAINING, Y);
disp('Training completed - Spatial features only');

for i = 41:58 
    testingblind(i - 40,:) = [ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteo(:,:,i)) sfta(osteo(:,:,i),1) featureO(:,:,i) featureVectoro(:,:,i) extractLBPFeatures(osteo(:,:,i))];  
    testingblind(18 + i - 40,:) = [ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(control(:,:,i)) sfta(control(:,:,i),1) featureC(:,:,i) featureVectorc(:,:,i) extractLBPFeatures(control(:,:,i))];
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
Y = [zeros(1,40), ones(1,40)]; %Y - labels for training data for SVM classifier.
Yada = [-ones(1,40), ones(1,40)]; %Yada - labels for training data for Adaboost/Random Forest classifier.
Y1blind = [ones(1,18) zeros(1,18)]; %labels for the blind data for SVM classifer
Y1blindada = [ones(1,18) -ones(1,18)]; %labels for the blind data for Adaboost and Random Forest classifier..
clear TRAINING, clear testingblind, clear resultblindtree;
for i = 1:40
    TRAINING(i,:) = [avggc(:,:,i)]; %control
    TRAINING(i+40,:) = [avggo(:,:,i)];
end

[estimateclasstotal,model]=adaboost('train',TRAINING,Yada,50);
TreeObject=TreeBagger(30,TRAINING,Yada,'method','classification','NVarToSample','all');
SVMStruct = fitcsvm(TRAINING, Y);
disp('Training completed - Directional transform based features only');

for i = 41:58 
    testingblind(i - 40,:) = [avggo(:,:,i)];  
    testingblind(18 + i - 40,:) = [avggc(:,:,i)];
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
Y = [zeros(1,40), ones(1,40)]; %Y - labels for training data for SVM classifier.
Yada = [-ones(1,40), ones(1,40)]; %Yada - labels for training data for Adaboost/Random Forest classifier.
Y1blind = [ones(1,18) zeros(1,18)]; %labels for the blind data for SVM classifer
Y1blindada = [ones(1,18) -ones(1,18)]; %labels for the blind data for Adaboost and Random Forest classifier..
clear TRAINING, clear testingblind, clear resultblindtree;
for i = 1:40
    TRAINING(i,:) = [chaar(i,:) min(min(dct2(double(control(:,:,i)))))]; %control
    TRAINING(i+40,:) = [ohaar(i,:) min(min(dct2(double(osteo(:,:,i)))))];
end

[estimateclasstotal,model]=adaboost('train',TRAINING,Yada,50);
TreeObject=TreeBagger(30,TRAINING,Yada,'method','classification','NVarToSample','all');
SVMStruct = fitcsvm(TRAINING, Y);
disp('Training completed - Non-directional transform based features only');

for i = 41:58 
    testingblind(i - 40,:) = [ohaar(i,:) min(min(dct2(double(osteo(:,:,i)))))];  
    testingblind(18 + i - 40,:) = [chaar(i,:) min(min(dct2(double(control(:,:,i)))))];
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
       
  if x == 2
      fprintf('\n\n\n') 
  fprintf('NOISY - GAUSSIAN.\n\n')   
  tic
  clearvars -except x;
  [control, osteo] =  getimages(x); %Refer readme file to set path for getimages()
  gaborArray = gaborFilterBank(2,2,39,39); % for general filterbank
% Getting features..
 for i = 1:58 %658
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
Y = [zeros(1,40), ones(1,40)]; %Y - labels for training data for SVM classifier.
Yada = [-ones(1,40), ones(1,40)]; %Yada - labels for training data for Adaboost/Random Forest classifier.
Y1blind = [ones(1,18) zeros(1,18)]; %labels for the blind data for SVM classifer
Y1blindada = [ones(1,18) -ones(1,18)]; %labels for the blind data for Adaboost and Random Forest classifier..
  fprintf('Features have been extracted.\n')
% Training and testing...

%% ALL FEATURES
for i = 1:40
    TRAINING(i,:) = [extractLBPFeatures(control(:,:,i)) featureVectorc(:,:,i) sfta(control(:,:,i),1) chaar(i,:) avggc(:,:,i) min(min(dct2(double(control(:,:,i))))) ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(control(:,:,i)) featureC(:,:,i)]; %control
    TRAINING(i+40,:) = [extractLBPFeatures(osteo(:,:,i)) featureVectoro(:,:,i) sfta(osteo(:,:,i),1) ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteo(:,:,i))))) ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteo(:,:,i)) featureO(:,:,i)];
end

[estimateclasstotal,model]=adaboost('train',TRAINING,Yada,50);
TreeObject=TreeBagger(30,TRAINING,Yada,'method','classification','NVarToSample','all');
SVMStruct = fitcsvm(TRAINING, Y);
disp('Training completed for all features.');

for i = 41:58 
    testingblind(i - 40,:) = [extractLBPFeatures(osteo(:,:,i)) featureVectoro(:,:,i) sfta(osteo(:,:,i),1) ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteo(:,:,i))))) ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteo(:,:,i)) featureO(:,:,i)]; 
    testingblind(18 + i - 40,:) = [extractLBPFeatures(control(:,:,i)) featureVectorc(:,:,i) sfta(control(:,:,i),1) chaar(i,:) avggc(:,:,i) min(min(dct2(double(control(:,:,i))))) ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(control(:,:,i)) featureC(:,:,i)]; %control
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
  
  if x == 3
      fprintf('\n\n\n') 
  fprintf('WITH ANGULAR DISTORTION.\n\n')   
  tic
  clearvars -except x;
  [control, osteo] =  getimages(x); %Refer readme file to set path for getimages()
  gaborArray = gaborFilterBank(2,2,39,39); % for general filterbank
% Getting features..
 for i = 1:58 %658
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
Y = [zeros(1,40), ones(1,40)]; %Y - labels for training data for SVM classifier.
Yada = [-ones(1,40), ones(1,40)]; %Yada - labels for training data for Adaboost/Random Forest classifier.
Y1blind = [ones(1,18) zeros(1,18)]; %labels for the blind data for SVM classifer
Y1blindada = [ones(1,18) -ones(1,18)]; %labels for the blind data for Adaboost and Random Forest classifier..
  fprintf('Features have been extracted.\n')
% Training and testing...

%% ALL FEATURES
for i = 1:40
    TRAINING(i,:) = [extractLBPFeatures(control(:,:,i)) featureVectorc(:,:,i) sfta(control(:,:,i),1) chaar(i,:) avggc(:,:,i) min(min(dct2(double(control(:,:,i))))) ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(control(:,:,i)) featureC(:,:,i)]; %control
    TRAINING(i+40,:) = [extractLBPFeatures(osteo(:,:,i)) featureVectoro(:,:,i) sfta(osteo(:,:,i),1) ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteo(:,:,i))))) ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteo(:,:,i)) featureO(:,:,i)];
end

[estimateclasstotal,model]=adaboost('train',TRAINING,Yada,50);
TreeObject=TreeBagger(30,TRAINING,Yada,'method','classification','NVarToSample','all');
SVMStruct = fitcsvm(TRAINING, Y);
disp('Training completed for all features.');

for i = 41:58 
    testingblind(i - 40,:) = [extractLBPFeatures(osteo(:,:,i)) featureVectoro(:,:,i) sfta(osteo(:,:,i),1) ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteo(:,:,i))))) ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteo(:,:,i)) featureO(:,:,i)]; 
    testingblind(18 + i - 40,:) = [extractLBPFeatures(control(:,:,i)) featureVectorc(:,:,i) sfta(control(:,:,i),1) chaar(i,:) avggc(:,:,i) min(min(dct2(double(control(:,:,i))))) ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(control(:,:,i)) featureC(:,:,i)]; %control
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
  
  if x == 4
  fprintf('\n\n\n') 
  fprintf('WITH NOISE AND ANGULAR DISTORTION.\n\n')   
  tic
  clearvars -except x;
  [control, osteo] =  getimages(x); %Refer readme file to set path for getimages()
  gaborArray = gaborFilterBank(2,2,39,39); % for general filterbank
% Getting features..
 for i = 1:58 %658
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
Y = [zeros(1,40), ones(1,40)]; %Y - labels for training data for SVM classifier.
Yada = [-ones(1,40), ones(1,40)]; %Yada - labels for training data for Adaboost/Random Forest classifier.
Y1blind = [ones(1,18) zeros(1,18)]; %labels for the blind data for SVM classifer
Y1blindada = [ones(1,18) -ones(1,18)]; %labels for the blind data for Adaboost and Random Forest classifier..
  fprintf('Features have been extracted.\n')
% Training and testing...
%% ALL FEATURES
for i = 1:40
    TRAINING(i,:) = [extractLBPFeatures(control(:,:,i)) featureVectorc(:,:,i) sfta(control(:,:,i),1) chaar(i,:) avggc(:,:,i) min(min(dct2(double(control(:,:,i))))) ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(control(:,:,i)) featureC(:,:,i)]; %control
    TRAINING(i+40,:) = [extractLBPFeatures(osteo(:,:,i)) featureVectoro(:,:,i) sfta(osteo(:,:,i),1) ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteo(:,:,i))))) ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteo(:,:,i)) featureO(:,:,i)];
end

[estimateclasstotal,model]=adaboost('train',TRAINING,Yada,50);
TreeObject=TreeBagger(30,TRAINING,Yada,'method','classification','NVarToSample','all');
SVMStruct = fitcsvm(TRAINING, Y);
disp('Training completed for all features.');

for i = 41:58 
    testingblind(i - 40,:) = [extractLBPFeatures(osteo(:,:,i)) featureVectoro(:,:,i) sfta(osteo(:,:,i),1) ohaar(i,:) avggo(:,:,i) min(min(dct2(double(osteo(:,:,i))))) ocon(i)+ocor(i) oene(i) ohom(i) omean(i) ostd(i) entropy(osteo(:,:,i)) featureO(:,:,i)]; 
    testingblind(18 + i - 40,:) = [extractLBPFeatures(control(:,:,i)) featureVectorc(:,:,i) sfta(control(:,:,i),1) chaar(i,:) avggc(:,:,i) min(min(dct2(double(control(:,:,i))))) ccon(i)+ccor(i) cene(i) chom(i) cmean(i) cstd(i) entropy(control(:,:,i)) featureC(:,:,i)]; %control
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
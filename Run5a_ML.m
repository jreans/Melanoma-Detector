% %% add in all the attached functions
% 
% % Amanda: make sure you take out my comments when you correct something
% % they tell you to work on.  Delete when done.
% 
% % Amanda: write a word document using endnote to supply references for all your algorithms.  Look at the previous publication as a model.  Please describe the cost function in relation to overfitting.  
% 
clear all
close all

load CurrentWorkingDir

% folderDIR = 'S0_Orig_Marghoob';
cd data
cd(folderDIR)

load Matrix_mIBs_3
% imID is the row number that it appears in the matrix

% cd Track_mIBs

% d = dir('*.mat');
% [x,y]=sort(datenum(char({d.date})));
% most_recent_file=char(d(y(end)).name);
% matrixname = strcat(most_recent_file(1:12), 'Matrix_mIBs');
% gsname = strcat(most_recent_file(1:12),'GSdiag');
% 
% load(matrixname) % Matrix_Out must be a 2-d matrix containing statistically significant biomarkers as columns and images as rows

X = Matrix_Out3;

% old = [1;2;3;4;6;7;8;9;13;14;19;20;24;27;29;30;79;80;81;82;83;84];
% oldind = [1,2,3,5,6,7,11,12,14,19,20,21,27,29,31,33,34,35,36,37,38,40];
% newind = setxor(1:41,oldind);
%             
% X = X(:,newind);


% load(gsname)
load GSdiag

% imorder = zeros(length(imID),1);
%Marghoob
% for i = 1:length(imID)
%     for k = 1:size(Diagnoses,1)
%         if imID(i) == Diagnoses(k,1)
%             imorder(i,1) = k; % order # of pics according to Diagnoses
%         end
%     end
% end
% y = Diagnoses(imID,2);

%MSKCC
y = Diagnoses(imID,2);
y2 = y(y==0 | y==1); % shorter/faster method to extract diagnoses
clear y
y = y2;
clear y2

% problem with Puig imID is dif from # of
if strcmp(folderDIR,'S2_Puig') 
% pics in Matrix_Out3
    y = Diagnoses(imID',15);
    y2 = y(y==0 | y==1); % shorter/faster method to extract diagnoses
    
    counter = 0;
    for i = 1:length(Diagnoses)
        if Diagnoses(i,15) == 1 || Diagnoses(i,15) == 0
            counter = counter + 1;
            imID_2(counter) = i;
        end
    end
    imID3 = intersect(imID_2,imID);
    y = Diagnoses(imID3, 15);
    
end


cd ../..
%%
[M2,ptable] = pvalue3(X,y); % not actually using the outputs, just seeing getting some info on the biomarkers
name = folderDIR; % name of dataset
numiter = 1000; % number of iterations to select a different training and test set to train machine learning model
numalg = 8; % number of algorithms
LoRlambda = 1; % cost function for logistic regression -> higher cost function will penalize more against overfitting
NNlambda = 1; % cost function for neural networks -> higher cost function will penalize more against overfitting
imbalance = 1; % if dataset contains an even ratio of melanoma to nevi, set imbalance to 0; otherwise, set to 1 to select a subset
% containing an even ratio of melanoma/nevi to run machine learning
%% define X, y
idset = []; % filler for balanced datasets
XProcessed = zscore(X);
m = length(y);
v = [1:m]'; % NEED bracket to transpose array


mel = find(y == 1); %only used in imbalance = 1
lim = length(mel);
nev = find(y == 0);

%% Run first time just to debug

% XProcessed is normalized X for already balanced sets (imbalance == 0)
% X is not normalized, is later normalized in each subset (iteration) for imbalanced datasets, ASK
% JOEL IF THIS IS CORRECT APPROACH
ttscore = zeros(size(X,1),numalg,numiter); % Q scores of test set
[tscore,LoRcosttest] = ...
    F29_MLMasterCost(XProcessed,y,m,v,imbalance,X,idset,mel,nev,lim,numalg, ...
    NNlambda, LoRlambda); 

%% Run 1000 iteration using parfor
matdir = dir('*.mat');
mdir = dir('*.m');
p = gcp;
addAttachedFiles(p,{matdir.name,mdir.name});
% addAttachedFiles(p,{'Master','accuracy','checkNNGradients','costFunctionReg', ...
%     'costFunctionROCarea','LinearRegression','LinearRegressionTest', 'nnCostFunction', ...
%     'predictnn','randInitializeWeights'});


parfor iter = 1:numiter  % Dan changed from parfor to for to debug
    [tscore,LoRcosttest] = ...
        F29_MLMasterCost(XProcessed,y,m,v,imbalance,X,idset,mel,nev,lim,numalg, ...
    NNlambda, LoRlambda); 
    ttscore(:,:,iter) = tscore;
    tLoRcosttest(:,:,iter) = LoRcosttest;
end
ttscore = nanmean(ttscore,3);
tLoRcosttest = nanmean(tLoRcosttest,3);
%% Plot ROC curves
figure
stats = zeros(numalg+1,4);

% Dan added 
[Xalg,Yalg,T,AUC] = perfcurve(y,ttscore(:,1),'1');
ROC_X_methods = zeros(length(Xalg)  , size(ttscore,2));
ROC_Y_methods = ROC_X_methods;
for i_alg = 1:size(ttscore,2)
    [Xalg,Yalg,T,AUC] = perfcurve(y,ttscore(:,i_alg),'1');
    [sp,T] = F28_MLgetspwse98(Xalg,Yalg,T);
    stats(i_alg,:) = [i_alg,sp,T,AUC];
    plot(Xalg,Yalg)
    hold on
    
    ROC_X_methods(1:length(Xalg),i_alg) = Xalg; % check it out, Amanda
    ROC_Y_methods(1:length(Yalg),i_alg) = Yalg;
    
end
mttscore = mean(ttscore,2); % Q score
[Xalg,Yalg, T, AUC] = perfcurve(y,mttscore,'1');
plot(Xalg,Yalg,'k-','LineWidth',5)
legend('LoR','NN','SVM','DT','RF','LDA','KNN','NB','All');
xlabel('1-Specificity')
ylabel('Sensitivity')
title(name)


%%
cd Data
cd(folderDIR)

if ~isdir('Track_ML')
    mkdir('Track_ML')
end
cd Track_ML
filename = ['ROC_' folderDIR '_' datestr(now, 'dd-mmm-yyyy')];
saveas(gcf,filename, 'jpg')

StudyName = folderDIR;
save OutPut_ROCs X XProcessed y ttscore ROC_X_methods ROC_Y_methods Xalg Yalg StudyName

cd ..


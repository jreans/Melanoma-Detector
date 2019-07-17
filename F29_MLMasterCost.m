function [tscore,LoRcosttest] = ...
    F29_MLMasterCost(XProcessed,y,m,v,imbalance,X,idset,mel,nev,lim,numalg, ...
    NNlambda, LoRlambda); 
crossval = 0;
%%
% Define train, crossvalidation, and test sets
if crossval
    temp1 = randsample(m, round(m.*0.6));
    imtrain = XProcessed(temp1,:);
    truetrain = y(temp1);
    temp2 = setxor(temp1,v);
    temp3 = randsample(length(temp2), round(m.*0.2));
    temp4 = temp2(temp3);
    imval = XProcessed(temp4,:);
    trueval = y(temp4,:);
    temp5 = union(temp1,temp4);
    temp6 = setxor(temp5,v);
    imtest = XProcessed(temp6,:);
    truetest = y(temp6);
else
    % Define train and test sets
    if imbalance
        temp = randsample(length(nev),lim);
        nevid = nev(temp);
        idset = cat(1,mel,nevid);
        XProcessed = X(idset,:); % defines subset of images
        XProcessed = zscore(XProcessed); % standardize data
        y = y(idset,:); % ground truth for subset of images
        m = length(y);
        v = [1:m]'; % NEED bracket to transpose array

%         id = ones(size(XProcessed,1),1);
%         id(idtrain2) = 0;
%         iddtest2 = find(id == 1);
%         dtest2 = XProcessed(iddtest2,:);
%         % all the images in the MSKCC dataset outside of training
%         y2 = yProcessed(iddtest2);
%         dtscore = NaN(size(XProcessed,1),numalg);
    end
    idtrain = randsample(m, round(m.*0.8));
    imtrain = XProcessed(idtrain,:);
    truetrain = y(idtrain);
    idtest = setxor(idtrain,v);
    imtest = XProcessed(idtest,:);
    truetest = y(idtest);
    
    if imbalance
        idtest = idset(idtest); %reset to get correct indexing for tscore
        idtrain = idset(idtrain);
    else
        sscore = []; % filler for balanced datasets
        dtscore = []; % filler for balanced datasets
    end
end

tscore = NaN(size(X,1),numalg);
%%
% Logistic Regression
initial_theta = zeros(size(imtrain, 2)+1, 1);
% Vary regularization parameter
% Set Options
options = optimset('GradObj', 'on', 'MaxIter', 400);
% Optimize
[theta, costtrain] = ...
    fminunc(@(t)(F15_MLcostFunctionReg(t, imtrain, truetrain, LoRlambda)), initial_theta, options);
% Predict, Accuracy
temp = [ones(size(imtrain,1),1),imtrain];
predtrain = F24_MLsigmoid(temp*theta)>=0.5;
%[LoRactrain, LoRsetrain, LoRsptrain, LoRtrain] = accuracy(predtrain, truetrain);
%temp = [ones(size(imval,1),1),imval];
%predval = sigmoid(temp*theta)>=0.5;
%[LoRacval, LoRseval, LoRspval] = accuracy(predval, trueval);
[LoRcosttest, gradtest] = F15_MLcostFunctionReg(theta, imtest, truetest, LoRlambda);
imtest2 = [ones(size(imtest,1),1),imtest];
LoRpredtest = F24_MLsigmoid(imtest2*theta)>=0.5;
%[LoRactest, LoRsetest, LoRsptest, LoRtest] = accuracy(LoRpredtest, truetest);
tscore(idtest,1) = F24_MLsigmoid(imtest2*theta);
% %%
% % Neural Networks
% input_layer_size = size(imtrain,2);
% hidden_layer_size = ceil(size(imtrain,2)./2); %40;
% num_labels = 2;
% initial_Theta1 = F22_MLrandInitializeWeights(input_layer_size, hidden_layer_size);
% initial_Theta2 = F22_MLrandInitializeWeights(hidden_layer_size, num_labels);
% initial_nn_params = [initial_Theta1(:) ; initial_Theta2(:)];
% % checkNNGradients(3);
% 
% options = optimset('MaxIter', 300);
% costFunction = @(p) F19_MLnnCostFunction(p, ...
%     input_layer_size, ...
%     hidden_layer_size, ...
%     num_labels, imtrain, truetrain, NNlambda);
% [nn_params, cost] = F27_MLfmincg(costFunction, initial_nn_params, options);
% Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
%     hidden_layer_size, (input_layer_size + 1));
% Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end), ...
%     num_labels, (hidden_layer_size + 1));
% [predtrain, NNcosttrain] = F21_MLpredictnn(Theta1, Theta2, imtrain,truetrain,NNlambda);
% %[NNactrain, NNsetrain, NNsptrain, NNtrain] = accuracy(predtrain, truetrain);
% %predval = predict(Theta1, Theta2, imval);
% %[NNacval, NNseval, NNspval] = accuracy(predval, trueval);
% [NNpredtest,NNcosttest] = F21_MLpredictnn(Theta1, Theta2, imtest,truetest,NNlambda);
% %[NNactest, NNsetest, NNsptest, NNtest] = accuracy(NNpredtest, truetest);
% tscore(idtest,2) = F23_MLscorenn(Theta1, Theta2, imtest);
% dscore(:,2) = F23_MLscorenn(Theta1, Theta2, XProcessed);
% if imbalance
%     sscore(:,2) = F23_MLscorenn(Theta1, Theta2, X);
%     dtscore(iddtest2,2) = F23_MLscorenn(Theta1, Theta2, dtest2);
% end
%%
% Neural Networks Matlab Package

NNimtrain = imtrain';
NNtruetrain = [truetrain, 1-truetrain]';
NNimtest = imtest';

net = patternnet(5);
net = train(net,NNimtrain,NNtruetrain);

testY = net(NNimtest);
tscore(idtest,2) = testY(1,:)';
%%
% SVM
mdlSVM = fitcsvm(imtrain,truetrain);
% mdlSVM = fitcsvm(imtrain,truetrain);
% mdlSVM = fitcsvm(imtrain,truetrain,'Standardize',true,'KernelFunction','RBF');
mdlSVM = fitPosterior(mdlSVM);
predtrain = predict(mdlSVM, imtrain);
% SVMStruct = svmtrain(imtrain, truetrain, 'kernel_function', 'rbf');
% predtrain = svmclassify(SVMStruct, imtrain);
%[SVMactrain, SVMsetrain, SVMsptrain, SVMtrain] = accuracy(predtrain, truetrain);
%Group = svmclassify(SVMStruct, imval);
%[SVMacval, SVMseval, SVMspval] = accuracy(Group, trueval);
[SVMpredtest,score] = predict(mdlSVM, imtest);
tscore(idtest,3) = score(:,2);
%%
% Decision Tree
tree = fitctree(imtrain, truetrain,'Cost',[0,1;2,0]);
predtrain = predict(tree, imtrain);
%[DTactrain, DTsetrain, DTsptrain, DTtrain] = accuracy(predtrain, truetrain);
%label = predict(tree, imval);
%[DTacval, DTseval, DTspval] = accuracy(label, trueval);
[DTpredtest,score] = predict(tree, imtest);
%[DTactest, DTsetest, DTsptest, DTtest] = accuracy(DTpredtest, truetest);
tscore(idtest,4) = score(:,2);
%%
% Random Forest
B = TreeBagger(10,imtrain, truetrain,'Cost',[0,1;2,0]);
[label,score,cost] = predict(B, imtest);
[M, I] = max(score, [], 2);
I(I == 1) = 0;
I(I == 2) = 1;
RFpredtest = I;
%[~, RFsetest, RFsptest, RFtest] = accuracy(I, truetest);
tscore(idtest,5) = score(:,2);
%%
% LDA

%MdlLinear = fitcdiscr(imtrain,truetrain);
MdlLinear = fitcdiscr(imtrain,truetrain,'discrimType', 'pseudoLinear');  % Dan Changed to get it to run
[~,score] = predict(MdlLinear,imtest);
tscore(idtest,6) = score(:,2);
%%
% K nearest neighbors

mdlknn = fitcknn(imtrain,truetrain,'NumNeighbors',5);
[~,score] = predict(mdlknn,imtest);
tscore(idtest,7) = score(:,2);
%%
% % Naive Bayes
% Naive Bayes can also be programmed to consider prior probabilities
% instead of using the imbalance simulation
% idtrain = randsample(m, round(m.*0.8));
% imtrain = XProcessed(idtrain,:);
% truetrain = y(idtrain);
% idtest = setxor(idtrain,v);
% imtest = XProcessed(idtest,:);
% truetest = y(idtest);
try  % Dan added Try to get it to run
Mdlnb = fitcnb(imtrain,truetrain);
[~,score] = predict(Mdlnb,imtest);
tscore(idtest,8) = score(:,2);
catch 
end
%%
% % SVM with radial kernel
% mdlSVM = fitcsvm(imtrain,truetrain,'KernelFunction','RBF','Cost',[0,1;2,0]);
% mdlSVM = fitPosterior(mdlSVM);
% predtrain = predict(mdlSVM, imtrain);
% [SVMpredtest,score] = predict(mdlSVM, imtest);
% tscore(idtest,7) = score(:,2);
% [~,score] = predict(mdlSVM, XProcessed);
% dscore(:,7) = score(:,2);
% if imbalance
%     [~,score] = predict(mdlSVM, X);
%     sscore(:,7) = score(:,2);
%     [~,score] = predict(mdlSVM, dtest2);
%     dtscore(iddtest2,7) = score(:,2);
% end
%%
% % Logistic Regression with cost function as maximize ROC area
% initial_theta = zeros(size(imtrain, 2)+1, 1);
% % Vary regularization parameter
% % Set Options
% options = optimset('GradObj', 'on', 'MaxIter', 400);
% % Optimize
% [theta, costtrain] = ...
% 	fminunc(@(t)(costFunctionROCarea(t, imtrain, truetrain, LoR2lambda)), initial_theta, options);
% % Predict, Accuracy
% temp = [ones(size(imtrain,1),1),imtrain];
% predtrain = sigmoid(temp*theta)>=0.5;
% [LoR2actrain, LoR2setrain, LoR2sptrain, LoR2train] = accuracy(predtrain, truetrain);
% %temp = [ones(size(imval,1),1),imval];
% %predval = sigmoid(temp*theta)>=0.5;
% %[LoRacval, LoRseval, LoRspval] = accuracy(predval, trueval);
% costtest = costFunctionReg(theta, imtest, truetest, LoR2lambda);
% imtest2 = [ones(size(imtest,1),1),imtest];
% LoRpredtest = sigmoid(imtest2*theta)>=0.5;
% [LoR2actest, LoR2setest, LoR2sptest, LoR2test] = accuracy(LoRpredtest, truetest);
% tscore(idtest,6) = sigmoid(imtest2*theta);
% XProcessed2 = [ones(size(XProcessed,1),1),XProcessed];
% dscore(:,6) = sigmoid(XProcessed2*theta);
% if imbalance
%     X2 = [ones(size(X,1),1),X];
%     sscore(:,6) = sigmoid(X2*theta);
%     X3 = [ones(size(dtest2,1),1),dtest2];
%     dtscore(iddtest2,6) = sigmoid(X3*theta);
% end
%%
% %%
% % Formats stats about each algorithm
% train = [LiRtrain, LoRtrain, NNtrain, SVMtrain, DTtrain, RFtrain];
% test = [LiRtest, LoRtest, NNtest, SVMtest, DTtest, RFtest];
% total = [train; test];
end
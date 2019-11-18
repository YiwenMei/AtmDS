% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 9/6/2019

%% Functionality
% This function evokes the Random Forest algorithm (RF, Breiman, 2001) for the
% downscaling of precipitation using the algorithm developed by Mei et al. 2019.
% The function has three different modes - predictor variable selection, model
% evaluation, and prediction.
%  1)Predictor Variable Selection Mode is used to determine the predition value
%    of variables on precipitation rate. Measure of prediction value for potential
%    predictors are implemented by either the permuted error or the decrease
%    in impurity due to split of predictors. Based on the results, users can
%    reduce the model complexity by removing variables with relative low prediction
%    value.
%  2)Model Evaluation Mode evokes RF to train classification or regression tree
%    model for precipitation mask or precipitation rate with the selected predictors
%    from Selection mode.
%  3)Prediction Mode is used to estimate the binary precipitation mask or precipitation
%    rate with the trained model from Evaluation mode.

%% Input
% md: mode of the function (6 modes) specified as a char array.
%      - Select PermuteError: use the permuted error as measure of predictor
%      importance metrics in predictor selection;
%      - Select DecreaseImpure: save as the previous but uses the decrease in
%      impurity instead;
%      - Evaluate mask: train the RF model for precipitation mask prediction
%       using Random Forest classification;
%       - Evaluate rate: evaluate the model for precipitation rate prediction
%       using Random Forest regreesion;
%       - Predict mask: predict the precipitation mask using the Random Forest
%       classification model attained from "Evaluate mask";
%       - Predict rate: predict the precipitation rate using the Random Forest
%       regression model attained from "Evaluate rate";

% X: for "Select", this is the potential predictors for precipitation rate;
%    for "Evaluate" and "Predict", this is the selected predictors for precipitation
%     mask or rate;

% Y: for "Select", this is the precipitation rate;
%    for "Evaluate", this is the precipitation mask or precipitation rate;
%    for "Predict" mode, this is the Treebagger object stores the RF classification
%     or regression model attained from the "Evaluate" mode;

% pr1: for "Select" and "Evaluate", this is the index of the witholding sample,
%       the day of year, and the location;
%      for "Predict", this is not requried (for "Predict rate", this can be supplied
%       as the residual of precipitation rate for a time step if the residual
%       is added back to the prediction);

% pr2: for "Select" and "Evaluate", this is the number of tree grown for the
%       RF model;
%      for "Predict mask", this is a binary mask indexing the no-data value of
%       a time step;
%      for "Predict rate", this is a binary mask indexing the "wet pixel" attained
%       from "Predict mask" for a time step;

% pr3: for "Select" and "Evaluate", this is the minimum leaf size;
%      for "Predict", this is not required;

% pr4: for "Select" and "Evaluate", this is a flag to use parallel computing
%       (set to "UseParallel" for parallel computing and "[]" for not using);
%      for "Predict", this is not required;

% pth: for "Select", this is not requried;
%      for "Evaluate", this is the full name of file to store the RF model object;
%      for "Predict", this is the working directory for the code.

%% Output
% O1: for "Select", O1 has two fields O1.PII and O1.I.
%      - O1.PII stores the predictor importance indices measured by the Predictor
%      Permuted Delta Error & the Delta Criterion Decision Split;
%      - O1.I stores the names of predictor sorted by the order of their removal.

%     for "Evaluate Mask", O1 has two fields O1.EM and O1.ss.
%      - O1.EM is the miss classification rate as a function of number of tree
%      for the in-bag, out-of-bag, and validation samples;
%      - O1.ss is the three sample sizes.

%     for "Evaluate Mask", O1 has three fields O1.EM, O1.ysts and O1.ss.
%      - O1.EM is the mean squared error (MSE) for the three samples;
%      - O1.ysts is the mean and variance of the three sample of the response;
%      - O1.ss is the three sample sizes.

%     for "Predict", this is the precipitation mask or precipitation rate for
%      a time step.

% O2: for "Select", this is an NT-by-NP-by-3 matrix stores the MSE for the RF
%      regression models (NT and NP stand for the number of tree and predictor;
%      the 3 records are for the in-bag, out-of-bag, and validation sample);

%     for "Evaluate", this is the four situations of contingency for classification
%      and the model residual for regression;

%     for "Predict", this is not required.

function [O1,O2]=Prec_DS(md,X,Y,pr1,pr2,pr3,pr4,pth)
if strncmp(md,'Select',6)
%% Variable selection mode
  I=cell(1,size(X,2));
  EM=nan(size(X,2),pr2,3);

% Model fitting
  CT=cell(size(I));
  PII=cell(size(I));
  for nt=1:size(X,2)
    CT{nt}=zeros(size(X,2),1);
    for t=1:2*size(X,2)+1 % Reduce the randomness by fitting multiple times
      if strcmp(pr4,'UseParallel')
        paropt=statset('UseParallel',true);
        Mdl=TreeBagger(pr2,X(pr1,:),Y(pr1,:),'Method','regression','OOBPredictorImportance','on','MinLeafSize',pr3,'Options',paropt);
      else
        Mdl=TreeBagger(pr2,X(pr1,:),Y(pr1,:),'Method','regression','OOBPredictorImportance','on','MinLeafSize',pr3);
      end

      if contains(md,'PermuteError')
        PI=Mdl.OOBPermutedPredictorDeltaError; % MDA
      elseif contains(md,'DecreaseImpure')
        PI=Mdl.DeltaCriterionDecisionSplit; % MDI
      end
      PII{nt}=[PII{nt} PI'];
      [~,oi]=min(PI);
      CT{nt}(oi)=CT{nt}(oi)+1;
      if max(CT{nt})>=3 || size(X,2)==1
        [~,oim]=max(CT{nt}); % ID of the least important variable
        break;
      end
    end

% MSE of IB, OoB, and Val sample
    EM(nt,:,1)=error(Mdl.compact,Mdl.X,Mdl.Y,'weights',Mdl.W,'useifort',~Mdl.OOBIndices);
    EM(nt,:,2)=oobError(Mdl);
    EM(nt,:,3)=error(Mdl.compact,X(~pr1,:),Y(~pr1,1),'weights',ones(length(find(~pr1)),1)/length(find(~pr1)));
    if nt~=length(I)
      clear Mdl
    end

% Removal of variable with the lowest PI
    I(nt)=X.Properties.VariableNames(oim);
    X(:,oim)=[];

% Outputs
  end

  y_ib=[];
  y_oob=[];
  ysts=nan(pr2,4);
  for t=1:pr2
    y_ib=[y_ib;Mdl.Y(~Mdl.OOBIndices(:,t))];
    y_oob=[y_oob;Mdl.Y(Mdl.OOBIndices(:,t))];
    ysts(t,:)=[mean(y_ib) var(y_ib) mean(y_oob) var(y_oob)];
  end
  ysts=[ysts repelem([mean(Y{~pr1,1}) var(Y{~pr1,1})],pr2,1)];
  ysts=array2table(ysts,'VariableNames',{'mean_yib','var_yib','mean_yoob','var_yoob','mean_yval','var_yval'});

  O1.PII=PII; % Predictor importance index
  O1.I=I; % Predictor removed per round
  O2.EM=EM; % MSE of model
  O2.ysts=ysts; % statistics of y

elseif strncmp(md,'Evaluate',8)
%% Function evaluation mode
  if contains(md,'mask')
% Model fitting for precipitation mask
    if strcmp(pr4,'UseParallel')
      paropt=statset('UseParallel',true);
      Mdl=TreeBagger(pr2,X(pr1,:),Y(pr1,:),'Method','classification','OOBPredictorImportance','on','MinLeafSize',pr3,'Options',paropt);
    else
      Mdl=TreeBagger(pr2,X(pr1,:),Y(pr1,:),'Method','classification','OOBPredictorImportance','on','MinLeafSize',pr3);
    end

% Miss Classification Rate
    O1.EM=nan(pr2,3);
    O1.EM(:,1)=error(Mdl.compact,Mdl.X,Mdl.Y,'weights',Mdl.W,'useifort',~Mdl.OOBIndices);
    Nib=length(find(~Mdl.OOBIndices))/pr2;
    O1.EM(:,2)=oobError(Mdl);
    Noob=length(find(Mdl.OOBIndices))/pr2;
    O1.EM(:,3)=error(Mdl.compact,X(~pr1,:),Y(~pr1,1),'weights',ones(length(find(~pr1)),1)/length(find(~pr1)));
    Nval=length(find(~pr1));
    O1.EM=array2table(O1.EM,'VariableNames',{'MCR_ib','MCR_oob','MCR_val'});
    O1.ss=[Nib Noob Nval];

    Mdl=compact(Mdl);

% Situations of contingency
    Yh=categorical(Mdl.predict(X));
    O2=nan(size(Y));
    Y=table2array(Y);
    O2(Yh=='0' & Y=='1')=1; % Missing (M)
    O2(Yh=='0' & Y=='0')=2; % Correct Negative (N)
    O2(Yh=='1' & Y=='0')=3; % False Alarm (F)
    O2(Yh=='1' & Y=='1')=4; % Hit (H)

  elseif contains(md,'rate')
% Model fitting for precipitation rate
    if strcmp(pr4,'UseParallel')
      paropt=statset('UseParallel',true);
      Mdl=TreeBagger(pr2,X(pr1,:),Y(pr1,:),'Method','regression','OOBPredictorImportance','on','MinLeafSize',pr3,'Options',paropt);
    else
      Mdl=TreeBagger(pr2,X(pr1,:),Y(pr1,:),'Method','regression','OOBPredictorImportance','on','MinLeafSize',pr3);
    end

% Mean Squared Error
    O1.EM=nan(pr2,3);
    O1.EM(:,1)=error(Mdl.compact,Mdl.X,Mdl.Y,'weights',Mdl.W,'useifort',~Mdl.OOBIndices);
    Nib=length(find(~Mdl.OOBIndices))/pr2;
    O1.EM(:,2)=oobError(Mdl);
    Noob=length(find(Mdl.OOBIndices))/pr2;
    O1.EM(:,3)=error(Mdl.compact,X(~pr1,:),Y(~pr1,1),'weights',ones(length(find(~pr1)),1)/length(find(~pr1)));
    Nval=length(find(~pr1));
    O1.EM=array2table(O1.EM,'VariableNames',{'MSE_ib','MSE_oob','MSE_val'});
    O1.ss=[Nib Noob Nval];

    y_ib=[];
    y_oob=[];
    ysts=nan(pr2,4);
    for t=1:pr2
      y_ib=[y_ib;Mdl.Y(~Mdl.OOBIndices(:,t))];
      y_oob=[y_oob;Mdl.Y(Mdl.OOBIndices(:,t))];
      ysts(t,:)=[mean(y_ib) var(y_ib) mean(y_oob) var(y_oob)];
    end
    ysts=[ysts repelem([mean(Y{~pr1,1}) var(Y{~pr1,1})],pr2,1)];
    ysts=array2table(ysts,'VariableNames',{'mean_yib','var_yib','mean_yoob','var_yoob','mean_yval','var_yval'});
    O1.ysts=ysts;

    Mdl=compact(Mdl);

% Prediction of residual
    Yh=Mdl.predict(X);
    O2=Y{:,1}-Yh;
  end

% Output the RF object
  if ~isempty(pth)
    [pth,nm,~]=fileparts(pth);
    save(fullfile(pth,nm),'Mdl','-v7.3');
  end

elseif strncmp(md,'Predict',7)
%% Prediction mode
  if contains(md,'mask')
% Prediction of precipitation mask
    Pmk=Y.predict(X);
    Pmk=cellfun(@str2double,Pmk);
    O1=nan(size(pr2));
    O1(pr2)=Pmk; % Precipitation mask at high resolution

  elseif contains(md,'rate')
% Prediction of precipitation rate
    Pr=Y.predict(X);
    O1=nan(size(pr2));
    O1(pr2==1)=Pr; % no adjustment of residual

    if ~isempty(pr1) % adjust for residual
      id=find(isnan(pr1)); % Pixels with NaN
% Fill the NaN pixels
      if ~isempty(id)
        xq=fix((id-1)/size(pr1,1))+1;
        yq=id-size(pr1,1)*(xq-1);
        id=find(~isnan(pr1)); % Other pixels
        x=fix((id-1)/size(pr1,1))+1;
        y=id-size(pr1,1)*(x-1);
        [id,d]=knnsearch([x y],[xq yq],'K',4);
        d=d.^2;
        d=d./repmat(sum(d,2),1,size(d,2)); % Weighting factor
        res=pr1(~isnan(pr1));
        res=sum(res(id).*d,2);
        pr1(isnan(pr1))=res;
      end

      pr1=imresize(pr1,size(O1),'bilinear');
      O1=O1+pr1;
    end

% Apply the mask
    O1(pr2==0)=NaN;
  end
  O2=[];
end
end

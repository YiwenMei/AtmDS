% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 11/06/2018

%% Functionality
% This function evokes the Random Forest algorithm (RF, Breiman, 2001) for the
% downscaling of precipitation rate using the algorithm developed by Mei et al. 2018.
% The function has three different modes - predictor variable selection, model
% evaluation, and prediction.
%  1)Predictor Variable Selection Mode is used to determine the predition value
%    of variables on precipitation rate. Measure of prediction value for potential
%    predictors are implemented by either the permuted error or the decrease
%    in impurity due to split of predictors. Based on the results, users can
%    reduce the model complexity by removing variables with relative low prediction
%    value.
%  2)Model Evaluation Mode evokes RF to create classification and regression
%    tree model for precipitation mask and precipitation rate over "wet pixel".
%  3)Prediction Mode is used to estimate the binary precipitation mask and values
%    of precipitation rate.

%% Input
%  md : mode of the function (can be "Select PermuteError", "Select DecreaseImpure",
%       "Evaluate mask", "Evaluate rate", "Predict mask" or "Predict rate". Their
%       meanings are as follow.
%        - Select PermuteError: use the permuted error as measure of predictor
%          importance in predictor selection;
%        - Select DecreaseImpure: use the decrease in impurity due to slit as
%          predictor importance measure in predictor selection;
%        - Evaluate mask: evaluate the model for precipitation mask prediction
%          using Random Forest classification;
%        - Evaluate rate: evaluate the model for precipitation rate prediction
%          using Random Forest regreesion;
%        - Predict mask: predict the precipitation mask using the Random Forest
%          classification model attained from "Evaluate mask";
%        - Predict rate: predict the precipitation rate using the Random Forest
%          regression model attained from "Evaluate rate";
%  X : (potential) predictors for precipitation mask or precipitation rate;
%  Y : for the two "Select" modes, this is the precipitation rate; for the two
%      "Evaluate" mode, this is the precipitation mask or precipitation rate;
%      for "Predict" mode, this is the Treebagger object stores the RF classification
%      or regression model attained from the "Evaluate" mode;
% pr1: in "Select" and "Evaluate", this is the index of the witholding sample
%      used in model validation; in "Predict rate", this is the residual of precipitation
%      rate for a time step (set it to "[]" in the case of "Predict mask");
% pr2: in "Select" and "Evaluate", this is the number of tree grown for the
%      RF model; in "Predict mask", this is a binary mask indexing the no-data
%      value of a time step; in "Predict rate", this is a binary mask indexing
%      the "wet pixel" attained from "Predict mask" for a time step;
% pr3: in "Select" and "Evaluate", this is the min leaf size parameter (set to
%      "[]" for "Predict");
% pr4: in "Select" and "Evaluate", this is a flag to use parallel computing (set
%      to "UseParallel" for parallel computing and "[]" for not using); set it
%      to "[]" for "Predict";
% pth: in the two "Evaluate" modes, this is the full name of file to store the
%      RF model; in the two "Predict" mode, this is the working directory for
%      the code (set to "[]" in "Select" modes).

%% Output
% O1: in "Select", O1 has two fields O1.PII and O1.I.
%      - O1.PII stores the predictor importance indices measured by the Predictor
%        Permuted Delta Error (PPDE) & the Delta Criterion Decision Split (DCDS);
%      - O1.I stores the names of predictor sorted by the order of their removal.
%     in "Evaluate", O1 has two fields O1.cal and O1.val.
%      - O1.cal stores the model performance metrics (additional note) for the
%        calibration sample;
%      - O1.val stores the model performance metrics for the validation sample;
%     in "Predict", this is the precipitation mask or precipitation rate for
%     a time step.
% O2: in "Select", this is the mean Square error (MSE) for the regression tree
%     model as a function of number of predictor and number of tree;
%     in "Evaluate", this is the four situations of contingency for classification
%     and the model residual for regression;
%     in "Predict", set it to "~".

%% Additional note
% Model performance metrics:
%  1)for classification tree model, O2.cal and O2.val stores the number of sample,
%    rate of missing, correct negative, false alarm, and hit;
%  2)for regression tree model, O2.cal and O2.val stores the number of sample,
%    mean relative error, centered root mean square error, correlation coefficient,
%    Nash Sutcliffe Efficiency and Kling-Gupta Efficiency.,'PredictorSelection','curvature'

function [O1,O2]=Prec_DS(md,X,Y,pr1,pr2,pr3,pr4,pth)
if strncmp(md,'Select',6)
%% Variable selection mode
  I=cell(1,size(X,2));
  Err=nan(size(X,2),pr2);

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

      if ~isempty(strfind(md,'PermuteError'))
        PI=Mdl.OOBPermutedPredictorDeltaError; % MDA
      elseif ~isempty(strfind(md,'DecreaseImpure'))
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

% MSE of OoB sample
    Err(nt,:)=oobError(Mdl);
    clear Mdl

% Removal of variable with the lowest PI
    I{nt}=X.Properties.VariableNames(oim);
    X(:,oim)=[];

% Outputs
  end
  O1.PII=PII; % Predictor importance index
  O1.I=I; % Predictor removed per round
  O2=Err; % MSE of model

elseif strncmp(md,'Evaluate',8)
%% Function evaluation mode
  if ~isempty(strfind(md,'mask'))
% Model fitting for precipitation mask
    if strcmp(pr4,'UseParallel')
      paropt=statset('UseParallel',true);
      Mdl=TreeBagger(pr2,X(pr1,:),Y(pr1,:),'Method','classification','OOBPredictorImportance','on','MinLeafSize',pr3,'Options',paropt);
    else
      Mdl=TreeBagger(pr2,X(pr1,:),Y(pr1,:),'Method','classification','OOBPredictorImportance','on','MinLeafSize',pr3);
    end

% Miss Classification Rate
    O1.ib=error(Mdl.compact,Mdl.X,Mdl.Y,'weights',Mdl.W,'useifort',~Mdl.OOBIndices);
    Nib=length(find(~Mdl.OOBIndices))/pr2;
    O1.oob=oobError(Mdl);
    Noob=length(find(Mdl.OOBIndices))/pr2;
    O1.val=error(Mdl.compact,X(~pr1,:),Y(~pr1,1),'weights',ones(length(find(~pr1)),1)/length(find(~pr1)));
    Nval=length(find(~pr1));
    O1.ss=[Nib Noob Nval];

    Mdl=compact(Mdl);

% Situations of contingency
    Yh=categorical(Mdl.predict(X));
    O2=nan(size(Y));
    O2(Yh=='0' & Y{:,1}=='1')=1; % Missing (M)
    O2(Yh=='0' & Y{:,1}=='0')=2; % Correct Negative (N)
    O2(Yh=='1' & Y{:,1}=='0')=3; % False Alarm (F)
    O2(Yh=='1' & Y{:,1}=='1')=4; % Hit (H)

  elseif ~isempty(strfind(md,'rate'))
% Model fitting for precipitation rate
    if strcmp(pr4,'UseParallel')
      paropt=statset('UseParallel',true);
      Mdl=TreeBagger(pr2,X(pr1,:),Y(pr1,:),'Method','regression','OOBPredictorImportance','on','MinLeafSize',pr3,'Options',paropt);
    else
      Mdl=TreeBagger(pr2,X(pr1,:),Y(pr1,:),'Method','regression','OOBPredictorImportance','on','MinLeafSize',pr3);
    end

% Mean Squared Error
    O1.ib=error(Mdl.compact,Mdl.X,Mdl.Y,'weights',Mdl.W,'useifort',~Mdl.OOBIndices);
    Nib=length(find(~Mdl.OOBIndices))/pr2;
    O1.oob=oobError(Mdl);
    Noob=length(find(Mdl.OOBIndices))/pr2;
    O1.val=error(Mdl.compact,X(~pr1,:),Y(~pr1,1),'weights',ones(length(find(~pr1)),1)/length(find(~pr1)));
    Nval=length(find(~pr1));
    O1.ss=[Nib Noob Nval];

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
  if ~isempty(strfind(md,'mask'))
% Prediction of precipitation mask
    Pmk=Y.predict(X);
    Pmk=cellfun(@str2double,Pmk);
    O1=nan(size(pr2));
    O1(pr2)=Pmk; % Precipitation mask at high resolution

  elseif ~isempty(strfind(md,'rate'))
% Prediction of precipitation rate
    Pr=Y.predict(X);
    O1=nan(size(pr2));
    O1(pr2==1)=Pr;

    if isempty(pr1) % no adjustment of residual
      O1=exp(O1);

    else % adjust for residual
      id=find(isnan(pr1)); % Pixels with NaN
      if ~isempty(id)
% Fill the NaN pixels
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
      O1=exp(O1+pr1);
    end

% Apply the mask
    O1(pr2==0)=0;
  end
  O2=[];
end
end

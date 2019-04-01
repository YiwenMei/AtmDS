% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 04/06/2018

%% Functionality
% This function has two functionalities:
%  1)calculates the gridded lapse rate for the study domain using the ambient
%    lapse rate method (Rouf et al. 2019);
%  2)calculates temperature at user-specified height using the grided lapse rate. 

%% Input
% Tfn : funll name of file or array for temperature data (K);
% Hmfn: full name of file or array for the measurement heights of the original
%       wind speed (m);
% Hnfn: full name of file or array for the new measurement height of wind (m);
%  Z  : elevation data over the same area (m);
% pflg: parallel flag (1 - use parallel channel, 0 - squential);
% ndv : no-data value in Z and Tfn (please use the same ndv for inputs);
%  sr : searching radius of block to perform linear regression;
% dcp : double-sided cut-off percentage used to cut off extreme positive/negative
%       values of lapse rate (2-by-1 with the 1st representing the low end and
%       the 2nd for the high end;

%% Output
% LR : lapse rate (K/m);
% EM : performance metrics of the regression model;
% Tn : temperature at new measurement height (K).

%% Additional note
% Order of EM: 1)coefficient of determination, r2=1-SSE/SST;
%              2)root mean square error, rmse=(SSE/(n-m))^0.5;
%              3)intersect of regression model.
% Example of dcp: say if dcp=[.05 .95], the lowest 5% and highest 5% LR values
%                 will be set to LR at the 5% and 95%.

function [LR,EM,Tn,Hn]=LapsRate(Tfn,Hmfn,Hnfn,Z,ndv,pflg,sr,dcp)
%% Check the inputs
switch nargin
    case {1,2,3,5}; error('Not enough number of arguments');
    case 6; sr=1; dcp=[.05 .95];
    case 7; dcp=[.05 .95];
    case 8
    otherwise; error('Too many number of arguments');
end

%% Read the inputs
if ischar(Hmfn)
  Hm=double(imread(Hmfn));
else
  Hm=imresize(Hmfn,size(Z));
end
if ischar(Hnfn)
  Hn=double(imread(Hnfn));
else
  Hn=imresize(Hnfn,size(Z));
end
if ischar(Tfn)
  Tm=double(imread(Tfn));
else
  Tm=Tfn;
end
clear Hmfn Hnfn Tfn
k=Z==ndv | Hm==ndv | Hn==ndv | Tm==ndv;
Z(k)=NaN;
Hm(k)=NaN;
Hn(k)=NaN;
Tm(k)=NaN;
Hn=Z+Hn;
Hm=Z+Hm;

%% Calculation of lapse rate
sdi=size(Z,1);
sdj=size(Z,2);
LR=nan(numel(Z),1);
r2=nan(numel(Z),1);
rms=nan(numel(Z),1);
itc=nan(numel(Z),1);

switch pflg
  case 1 % Use parallel channel
    parfor n=1:numel(Z)
      j=fix((n-1)/sdi)+1;
      i=n-sdi*(j-1);
      if i>sr && i<sdi-sr+1 && j>sr && j<sdj-sr+1 && ~isnan(Z(i,j))
        dHm=Hm(i-sr:sr:i+sr,j-sr:sr:j+sr)-Hm(i,j);
        dTm=Tm(i-sr:sr:i+sr,j-sr:sr:j+sr)-Tm(i,j);
        dHm=reshape(dHm,numel(dHm),1);
        dTm=reshape(dTm,numel(dTm),1);

% Linear regression to find lapse rate based on surrounding pixel
        [lr,~,res,~,stats]=regress(dTm,[ones(length(dHm),1) dHm]);
        r2(n)=stats(1); % r2
        rms(n)=sqrt(sum(res.^2)/(length(res)-length(lr))); % rmse
        itc(n)=lr(1); % Intersect
        LR(n)=lr(2); % Lapse rate
      end
    end

  case 0 % Squential
    for n=1:numel(Z)
      j=fix((n-1)/sdi)+1;
      i=n-sdi*(j-1);
      if i>sr && i<sdi-sr+1 && j>sr && j<sdj-sr+1 && ~isnan(Z(i,j))
        dHm=Hm(i-sr:sr:i+sr,j-sr:sr:j+sr)-Hm(i,j);
        dTm=Tm(i-sr:sr:i+sr,j-sr:sr:j+sr)-Tm(i,j);
        dHm=reshape(dHm,numel(dHm),1);
        dTm=reshape(dTm,numel(dTm),1);

% Linear regression to find lapse rate based on surrounding pixel
        [lr,~,res,~,stats]=regress(dTm,[ones(length(dHm),1) dHm]);
        r2(n)=stats(1); % r2
        rms(n)=sqrt(sum(res.^2)/(length(res)-length(lr))); % rmse
        itc(n)=lr(1); % Intersect
        LR(n)=lr(2); % Lapse rate
      end
    end
  otherwise; error('"pflg" is either 1 or 0');
end

% Reshape the results
LR=reshape(LR,sdi,sdj);
r2=reshape(r2,sdi,sdj);
rms=reshape(rms,sdi,sdj);
itc=reshape(itc,sdi,sdj);
LR=LR(1+sr:end-sr,1+sr:end-sr);
r2=r2(1+sr:end-sr,1+sr:end-sr,:);
rms=rms(1+sr:end-sr,1+sr:end-sr,:);
itc=itc(1+sr:end-sr,1+sr:end-sr,:);

% Control the boundary of LR
nthr=quantile(LR(~isnan(LR)),dcp(1));
pthr=quantile(LR(~isnan(LR)),dcp(2));
r2(LR<nthr | LR>pthr)=NaN;
rms(LR<nthr | LR>pthr)=NaN;
itc(LR<nthr | LR>pthr)=NaN;
EM=cat(3,r2,rms,itc);
LR(LR<nthr)=nthr;
LR(LR>pthr)=pthr;

%% Adjust the measurement height
Hn=Hn(1+sr:end-sr,1+sr:end-sr);
Hm=Hm(1+sr:end-sr,1+sr:end-sr);
Tm=Tm(1+sr:end-sr,1+sr:end-sr);
Z=Z(1+sr:end-sr,1+sr:end-sr);
Tn=LR.*(Hn-Hm)+Tm; % Adjust temperature
Hn=Hn-Z; % Change Hn to m a.g.
end

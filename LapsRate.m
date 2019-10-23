% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 10/16/2019

%% Functionality
% This function has two functionalities:
%  1)calculates the gridded lapse rate for the study domain using the ambient
%    lapse rate method (Rouf et al. 2019);
%  2)calculates temperature at user-specified height using the grided lapse rate. 

%% Input
% T : spatial map class (V2DCls.m) object or workspace variable for air or dew
%      point temperature (K);
% Hm: V2DCls.m object or workspace variable for measurement heights of temperature
%      (m above ground);
% Hn: V2DCls.m object or workspace variable for the new reference height of temperature
%      (m a.g.);
% Z : V2DCls.m object or workspace variable for terrain elevation (m above sea
%      level);

% pflg: parallel flag (false - default, squential; true - parallel);
%  sr : searching radius of block to perform linear regression (default is 1);
% dcp : double-sided cut-off percentage used to cut off extreme positive/negative
%        values of lapse rate (default is [.05 .95] where the lowest 5% and the
%        highest 95% will be cut off);

%% Output
% LR: lapse rate (K/m);
% EM: performance metrics of the regression model;
% Hn: new measurement height (m ag);
% Tn: temperature at new measurement height (K).

%% Additional note
% Order of EM: 1)coefficient of determination, r2=1-SSE/SST;
%              2)root mean square error, rmse=(SSE/(n-m))^0.5;
%              3)intersect of regression model.
% Require V2DCls.m.

function [LR,EM,Tn,Hn]=LapsRate(T,Hm,Hn,Z,varargin)
%% Check the inputs
narginchk(4,7);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'T',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'T'));
addRequired(ips,'Hm',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Hm'));
addRequired(ips,'Hn',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Hn'));
addRequired(ips,'Z',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Z'));

addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));
addOptional(ips,'sr',1,@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'sr'));
addOptional(ips,'dcp',[.05 .95],@(x) validateattributes(x,{'double'},{'numel',2,'>=',0,'<=',1},...
    mfilename,'dcp'));

parse(ips,T,Hm,Hn,Z,varargin{:});
pflg=ips.Results.pflg;
sr=ips.Results.sr;
dcp=ips.Results.dcp;
clear ips varargin

%% Read the inputs
Z=readCls(Z);
Hm=readCls(Hm);
Hm=imresize(Hm,size(Z));
Hn=readCls(Hn);
Hn=imresize(Hn,size(Z));
T=readCls(T);
k=isnan(Z) | isnan(Hm) | isnan(Hn) | isnan(T);
Z(k)=NaN;
Hm(k)=NaN;
Hn(k)=NaN;
T(k)=NaN;
Hn=Z+Hn;
Hm=Z+Hm;

%% Calculation of lapse rate
LR=nan(numel(Z),1);
r2=nan(numel(Z),1);
rms=nan(numel(Z),1);
itc=nan(numel(Z),1);

switch pflg
  case true % Use parallel
    parfor n=1:numel(Z)
      [lr,res,stats]=LapsRate_sub(n,size(Z),sr,Z,Hm,T)
      if ~isempty(lr)
        r2(n)=stats(1); % r2
        rms(n)=sqrt(sum(res.^2)/(length(res)-length(lr))); % rmse
        itc(n)=lr(1); % Intersect
        LR(n)=lr(2); % Lapse rate
      end
    end

  case false % Squential
    for n=1:numel(Z)
      [lr,res,stats]=LapsRate_sub(n,size(Z),sr,Z,Hm,T);
      if ~isempty(lr)
        r2(n)=stats(1); % r2
        rms(n)=sqrt(sum(res.^2)/(length(res)-length(lr))); % rmse
        itc(n)=lr(1); % Intersect
        LR(n)=lr(2); % Lapse rate
      end
    end
end

% Reshape the results
LR=reshape(LR,size(Z));
r2=reshape(r2,size(Z));
rms=reshape(rms,size(Z));
itc=reshape(itc,size(Z));
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
T=T(1+sr:end-sr,1+sr:end-sr);
Z=Z(1+sr:end-sr,1+sr:end-sr);
Tn=LR.*(Hn-Hm)+T; % Adjust temperature
Hn=Hn-Z; % Change Hn to m a.g.
end

function v2d=readCls(vb)
if isa(vb,'V2DCls')
  v2d=vb.readCls;
else
  v2d=vb;
end
end

function [lr,res,stats]=LapsRate_sub(n,sd,sr,Z,Hm,T)
% Linear regression to find lapse rate based on surrounding pixel
j=fix((n-1)/sd(1))+1;
i=n-sd(1)*(j-1);
if i>sr && i<sd(1)-sr+1 && j>sr && j<sd(2)-sr+1 && ~isnan(Z(i,j))
  dHm=Hm(i-sr:sr:i+sr,j-sr:sr:j+sr)-Hm(i,j);
  dTm=T(i-sr:sr:i+sr,j-sr:sr:j+sr)-T(i,j);
  dHm=reshape(dHm,numel(dHm),1);
  dTm=reshape(dTm,numel(dTm),1);

  [lr,~,res,~,stats]=regress(dTm,[ones(length(dHm),1) dHm]);

else
  lr=[];
  res=[];
  stats=[];
end
end

function [LR,Tn,Hn]=LR_Temp(mflg,T,Hm,Z,varargin)
%% Check the inputs
narginchk(4,7);
ips=inputParser;
ips.FunctionName=mfilename;

expmflg={'Slp_avg','Slp_min','Slp_reg'};
msg=cell2mat(cellfun(@(x) [x ', '],expmflg,'UniformOutput',false));
msg=sprintf('Expected mflg to be one of the following %s\n',msg);
addRequired(ips,'mflg',@(x) assert(any(strcmp(x,expmflg)),msg));
addRequired(ips,'T',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'T'));
addRequired(ips,'Hm',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Hm'));
addRequired(ips,'Z',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Z'));

addOptional(ips,'Hn',[],@(x) validateattributes(x,{'double','char'},{},mfilename,'Hn'));
addOptional(ips,'dcp',[.05 .95],@(x) validateattributes(x,{'double'},{'numel',2,'>=',0,'<=',1},...
    mfilename,'dcp'));
addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));

parse(ips,mflg,T,Hm,Z,varargin{:});
dcp=ips.Results.dcp;
Hn=ips.Results.Hn;
pflg=ips.Results.pflg;
clear ips varargin
sr=1;

%% Calculate the lapse rate
Z=readCls(Z);
Hm=readCls(Hm);
Hm=imresize(Hm,size(Z));
T=readCls(T);
Hm=Z+Hm;

ID=reshape(1:numel(T),size(T));
ID=ID(1+sr:end-sr,1+sr:end-sr);
ID=reshape(ID,numel(ID),1);
IDn=[ID-1,ID-1+size(Hm,1),ID+size(Hm,1),ID+1+size(Hm,1),ID+1,ID+1-size(Hm,1),ID-size(Hm,1),...
    ID-1-size(Hm,1)]; % N NE E SE S SW W NW
dT=T(IDn)-T(ID);
dH=Hm(IDn)-Hm(ID);
switch mflg
  case 'Slp_avg'
    LR=dT./dH;
    k=abs(LR)==min(abs(LR),[],2,'omitnan') | abs(LR)==max(abs(LR),[],2,'omitnan');
    LR(~k)=NaN;
    LR=nanmean(LR,2);

  case 'Slp_min'
    LR=dT./dH;
    k=abs(LR)==min(abs(LR),[],2,'omitnan') | abs(LR)==max(abs(LR),[],2,'omitnan');
    LR(~k)=NaN;
    k=abs(LR)==min(abs(LR),[],2,'omitnan');
    LR(~k)=NaN;
    LR=nansum(LR,2); % Collapse the 2nd dimension

  case 'Slp_reg'
    LR=nan(numel(Z)-2*sr,1);
    switch pflg
      case true % Use parallel
        parfor n=1:numel(Z)-2*sr
          LR(n)=LR_Temp_sub(dT,dH,n);
        end
      case false
        for n=1:numel(Z)-2*sr
          LR(n)=LR_Temp_sub(dT,dH,n);
        end
    end
end
LR=reshape(LR,size(Hm)-2*sr);

%% Control the boundary of LR
nthr=quantile(LR(~isnan(LR)),dcp(1));
pthr=quantile(LR(~isnan(LR)),dcp(2));
LR(LR<nthr)=nthr;
LR(LR>pthr)=pthr;

%% Adjust the measurement height
if ~isempty(Hn)
  Hn=readCls(Hn);
  Hn=imresize(Hn,size(Z));
  Hn=Z+Hn;
  Hn=Hn(1+sr:end-sr,1+sr:end-sr);
  Hm=Hm(1+sr:end-sr,1+sr:end-sr);
  T=T(1+sr:end-sr,1+sr:end-sr);
  Z=Z(1+sr:end-sr,1+sr:end-sr);
  Tn=Temp_Adj(T,Hm,Hn,LR);
  Hn=Hn-Z; % Change Hn to m a.g.
else
  Tn=T(1+sr:end-sr,1+sr:end-sr);
  Hn=Hm(1+sr:end-sr,1+sr:end-sr)-Z(1+sr:end-sr,1+sr:end-sr);
end
end

function v2d=readCls(vb)
if isa(vb,'char')
  v2d=matfile(vb);
  vb=cell2mat(who(v2d));
  eval(sprintf('v2d=v2d.%s;',vb));
else
  v2d=vb;
end
end

function [LR,itc,r2,rms]=LR_Temp_sub(dT,dH,n)
[lr,~,res,~,stats]=regress(dT(n,:)',[ones(length(dH(n,:)),1) dH(n,:)']);
r2(n)=stats(1); % r2
rms(n)=sqrt(sum(res.^2)/(length(res)-length(lr))); % rmse
itc(n)=lr(1); % Intersect
LR(n)=lr(2); % Lapse rate
end

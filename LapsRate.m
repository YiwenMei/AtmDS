% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 4/6/2018

%% Functionality
% This function has two functionalities:
%  1)calculates the gridded lapse rate for the study domain using the ambient
%    lapse rate method (Rouf et al. 2019);
%  2)calculates temperature at user-specified height using the grided lapse rate. 

%% Input
% T_fd : details of file or workspace variable for original temperature (K);
% Hm_fd: details of file or workspace variable for measurement heights of the original
%        temperature (m);
% Hn_fd: details of file or workspace variable for the new measurement height
%        of outputted temperature (m);
% Z_fd : details of file or workspace variable for coarse elevation (m asl);

% pflg: parallel flag (true - default, use parallel channel; false - squential);
%  sr : searching radius of block to perform linear regression (default is 1);
% dcp : double-sided cut-off percentage used to cut off extreme positive/negative
%       values of lapse rate (default is [.05 .95] where the lowest 5% and the
%       highest 95% will be cut off);

%% Output
% LR: lapse rate (K/m);
% EM: performance metrics of the regression model;
% Hn: new measurement height (m ag);
% Tn: temperature at new measurement height (K).

%% Additional note
% Order of EM: 1)coefficient of determination, r2=1-SSE/SST;
%              2)root mean square error, rmse=(SSE/(n-m))^0.5;
%              3)intersect of regression model.
% Require read2Dvar.m.

function [LR,EM,Tn,Hn]=LapsRate(T_fd,Hm_fd,Hn_fd,Z_fd,varargin)
%% Check the inputs
narginchk(4,7);
ips=inputParser;
ips.FunctionName=mfilename;
fprintf('%s received 4 required and %d optional inputs\n',mfilename,length(varargin));

addRequired(ips,'T_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'T_fd',1));
addRequired(ips,'Hm_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Hm_fd',2));
addRequired(ips,'Hn_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Hn_fd',3));
addRequired(ips,'Z_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Z_fd',4));

addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg',5));
addOptional(ips,'sr',1,@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'sr',6));
addOptional(ips,'dcp',[.05 .95],@(x) validateattributes(x,{'double'},{'row',2,'>=',0,'<=',1},mfilename,'dcp',7));

parse(ips,T_fd,Hm_fd,Hn_fd,Z_fd,varargin{:});
pflg=ips.Results.pflg;
sr=ips.Results.sr;
dcp=ips.Results.dcp;
clear ips

%% Read the inputs
Z=read2Dvar(Z_fd);
Hm=read2Dvar(Hm_fd);
Hm=imresize(Hm,size(Z));
Hn=read2Dvar(Hn_fd);
Hn=imresize(Hn,size(Z));
Tm=read2Dvar(T_fd);
k=isnan(Z) | isnan(Hm) | isnan(Hn) | isnan(Tm);
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
  case true % Use parallel channel
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

  case false % Squential
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

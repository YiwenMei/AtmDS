% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 11/9/2019

%% Functionality
% This function has two functionalities:
%  1)calculates temperature lapse rate for each grid cell with respect to its
%    eight neighbors and selects the minimum to represent lapse rate of that
%    grid cell;
%  2)adjusts temperature to a user-specified height using the lapse rate.
% The inputs must have an exactly one-grid buffer for the study dowmain.

%% Input
% T : spatial map class (V2DCls.m) object or workspace variable for air or dew
%      point temperature (K);
% Hm: V2DCls.m object or workspace variable for measurement heights of temperature
%      (m above ground);
% Z : V2DCls.m object or workspace variable for terrain elevation (m above sea
%      level);

% Hn : V2DCls.m object or workspace variable for the new reference height of
%       temperature (m a.g.);
% dcp: double-sided cut-off percentage used to cut off extreme positive/negative
%       values of lapse rate (default is [.05 .95] where the lowest 5% and the
%       highest 95% will be cut off).

%% Output
% LR: lapse rate (K/m);
% Tn: temperature at new measurement height (K);
% Hn: new measurement height (m ag).

%% Additional note
% Require V2DCls.m.

function [LR,Tn,Hn]=LR_Temp(T,Hm,Z,varargin)
%% Check the inputs
narginchk(3,5);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'T',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'T'));
addRequired(ips,'Hm',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Hm'));
addRequired(ips,'Z',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Z'));

addOptional(ips,'Hn',[],@(x) validateattributes(x,{'double','V2DCls'},{},mfilename,'Hn'));
addOptional(ips,'dcp',[.05 .95],@(x) validateattributes(x,{'double'},{'numel',2,'>=',0,'<=',1},...
    mfilename,'dcp'));

parse(ips,T,Hm,Z,varargin{:});
dcp=ips.Results.dcp;
Hn=ips.Results.Hn;
clear ips varargin

%% Calculate the lapse rate
Z=readCls(Z);
Hm=readCls(Hm);
Hm=imresize(Hm,size(Z));
T=readCls(T);
Hm=Z+Hm;

ID=reshape(1:numel(T),size(T));
ID=ID(2:end-1,2:end-1);
ID=reshape(ID,numel(ID),1);
IDn=[ID-1,ID-1+size(Hm,1),ID+size(Hm,1),ID+1+size(Hm,1),ID+1,ID+1-size(Hm,1)...
  ID-size(Hm,1),ID-1-size(Hm,1)]; % N NE E SE S SW W NW
LR=(T(IDn)-T(ID))./(Hm(IDn)-Hm(ID));
k=abs(LR)==min(abs(LR),[],2,'omitnan');
LR(~k)=NaN;
LR=nansum(LR,2);
LR=reshape(LR,size(Hm)-2*1);

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

  Hn=Hn(2:end-1,2:end-1);
  Hm=Hm(2:end-1,2:end-1);
  T=T(2:end-1,2:end-1);
  Z=Z(2:end-1,2:end-1);
  Tn=LR.*(Hn-Hm)+T; % Adjust temperature
  Hn=Hn-Z; % Change Hn to m a.g.

else
  Tn=T(2:end-1,2:end-1);
  Hn=Hm(2:end-1,2:end-1)-Z(2:end-1,2:end-1);
end
end

function v2d=readCls(vb)
if isa(vb,'V2DCls')
  v2d=vb.readCls;
else
  v2d=vb;
end
end

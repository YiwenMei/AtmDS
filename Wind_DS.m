% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 8/17/2019

%% Functionality
% This function downscales wind speed assuming log-wind profile.

%% Input
% ws: spatial wind class (Wind2DCls.m) object or workspace variable for total
%      or component wind (m/s);
% Hm: spatial variable class (V2DCls.m) object or workspace variable for measurement
%      height of wind (m ag);
% Hn: V2DCls.m object or workspace variable for new measurement height of outputted
%      wind (m ag);

% z0 : V2DCls.m object or workspace variable for original surface roughness (m ag);
% z0d: V2DCls.m object or workspace variable for downscaled surface roughness
%       (m ag);
% d0 : V2DCls.m object or workspace variable for original zero-plane displacement
%       height (m ag);
% d0d: V2DCls.m object or workspace variable for downscaled zero-plane displacement
%       height (m ag).

%% Output
% wsd: downscaled wind speed (m/s);

% wdd: downscaled wind direction (E is 0, counter-clock's wise is +);
% ud : downscaled eastward wind (m/s);
% vd : downscaled northward wind (m/s).

%% Additional note
% Require V2DCls.m.

function [wsd,wdd,ud,vd]=Wind_DS(ws,Hm,Hn,varargin)
%% Check inputs
narginchk(3,7);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'ws',@(x) validateattributes(x,{'double','V2DCls','Wind2DCls'},{'nonempty'},...
    mfilename,'ws'));
addRequired(ips,'Hm',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Hm'));
addRequired(ips,'Hn',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Hn'));

addOptional(ips,'z0',1,@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'z0'));
addOptional(ips,'z0d',1,@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'z0d'));
addOptional(ips,'d0',0,@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'d0'));
addOptional(ips,'d0d',0,@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'d0d'));

parse(ips,ws,Hm,Hn,varargin{:});
z0=ips.Results.z0;
z0d=ips.Results.z0d;
d0=ips.Results.d0;
d0d=ips.Results.d0d;
clear ips varargin

%% Read the wind records
if isa(ws,'Wind2DCls')
  wty='Component';
else
  wty='Total';
end
switch wty
  case 'Component' % U and V wind
    [ws,wd,~,~]=ws.readCls;
    wd(wd<0)=wd(wd<0)+360; % E is 0, counter-clock's wise is +

  case 'Total' % Total wind
    ws=readCls(ws);
end

Hm=readCls(Hm);
Hn=readCls(Hn);

%% Downscaling of wind speed/direction
if ~isa(z0,'scalar')
% Ratio of shear velocity
  z0=readCls(z0);
  z0d=readCls(z0d);

  if isempty(find(z0<=0, 1)) || isempty(find(z0d<=0, 1))
    z0i=imresize(z0,size(z0d),'bilinear');
    kus=(z0d./z0i).^.09; % Tao & Barrow (2018)
    clear z0i

% Ratio of height
    d0=readCls(d0);
    kh=(Hm-d0)./z0;
    if isempty(find(kh<=0 & ws>0, 1))
      kh(kh<=0)=exp(1);

      d0d=readCls(d0d);
      khd=(imresize(Hn,size(z0d),'bilinear')-d0d)./z0d;
      khd(khd<=0)=exp(1);
      kh=log(khd)./imresize(log(kh),size(z0d),'bilinear');
      kh(kh<0)=0;
      kh=kh/mean(kh(~isnan(kh)));
      kh=kh.*kus; % Combine the weighting factors
      clear khi z0 z0d d0 d0d khd kus

    else
      error('When Hm below d0, ws should be 0');
    end

  else
    error('z0 and z0d must be larger than 0');
  end

else
  kh=Hn./Hm; % kus=1
end
wsd=imresize(ws,size(kh),'bilinear').*kh;
clear kh Hn Hm ws

%% Output wind speed/direction
switch wty
  case 'Component' % U and V wind
    vd=wsd.*imresize(sind(wd),size(wsd),'bilinear');
    ud=wsd.*imresize(cosd(wd),size(wsd),'bilinear');
    wdd=atan2d(vd,ud); % E is 0, counter-clock's wise is +
    wdd(wdd<0)=wdd(wdd<0)+360;

  case 'Total'
    vd=[];
    ud=[];
    wdd=[];
end
end

function v2d=readCls(vb)
if isa(vb,'V2DCls')
  v2d=vb.readCls;
else
  v2d=vb;
end
end

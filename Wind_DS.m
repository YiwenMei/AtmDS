% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 8/17/2019

%% Functionality
% This function downscales wind speed assuming log-wind profile.

%% Input
% ws : details of file or workspace variable for original wind speed (m/s);
% wty: type of ws inputting (possible types are 'UV Wind' or 'Total Wind');
% Hm : details of file or workspace variable for measurement height of original
%      wind speed (m ag);
% Hn : details of file or workspace variable for new measurement height of outputted
%      wind speed (m ag);

% z0 : details of file or workspace variable for original surface roughness (m);
% z0d: details of file or workspace variable for downscaled surface roughness (m);
% d0 : details of file or workspace variable for original zero-plane displacement
%      height (m);
% d0d: details of file or workspace variable for downscaled zero-plane displacement
%      height (m).

%% Output
% wsd: downscaled wind speed (m/s);
% wdd: downscaled wind direction (E is 0, counter-clock's wise is +);
% ud : downscaled eastward wind (m/s);
% vd : downscaled northward wind (m/s).

%% Additional note
% Require read2Dvar.m.

function [wsd,wdd,ud,vd]=Wind_DS(ws,wty,Hm,Hn,varargin)
%% Check inputs
narginchk(4,8);
ips=inputParser;
ips.FunctionName=mfilename;
fprintf('%s received 4 required and %d optional inputs\n',mfilename,length(varargin));

addRequired(ips,'ws',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'ws',1));
addRequired(ips,'wty',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wty',2));
addRequired(ips,'Hm',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Hm',3));
addRequired(ips,'Hn',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Hn',4));

addOptional(ips,'z0',1,@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'z0',5));
addOptional(ips,'z0d',1,@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'z0d',6));
addOptional(ips,'d0',0,@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'d0',7));
addOptional(ips,'d0d',0,@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'d0d',8));

parse(ips,ws,wty,Hm,Hn,varargin{:});
z0=ips.Results.z0;
z0d=ips.Results.z0d;
d0=ips.Results.d0;
d0d=ips.Results.d0d;
clear ips varargin

%% Read the wind records
switch wty
  case 'UV Wind' % U and V wind
    u=read2Dvar(ws(1,:));
    v=read2Dvar(ws(2,:));

    ws=hypot(u,v);
    wd=atan2d(v,u); % E is 0, counter-clock's wise is +
    wd(wd<0)=wd(wd<0)+360;

  case 'Total Wind' % Total wind
    ws=read2Dvar(ws);
end

Hm=read2Dvar(Hm);
Hn=read2Dvar(Hn);
clear u v

%% Downscaling of wind speed/direction
if ~isscalar(z0)
% Ratio of shear velocity
  z0=read2Dvar(z0);
  z0d=read2Dvar(z0d);

  if isempty(find(z0<=0, 1)) || isempty(find(z0d<=0, 1))
    z0i=imresize(z0,size(z0d),'bilinear');
    kus=(z0d./z0i).^.09; % Tao & Barrow (2018)
    clear z0i

% Ratio of height
    d0=read2Dvar(d0);
    kh=(Hm-d0)./z0;
    if isempty(find(kh<=0 & ws>0, 1))
      kh(kh<=0)=exp(1);

      d0d=read2Dvar(d0d);
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
  case 'UV Wind' % U and V wind
    vd=wsd.*imresize(sind(wd),size(wsd),'bilinear');
    ud=wsd.*imresize(cosd(wd),size(wsd),'bilinear');
    wdd=atan2d(vd,ud); % E is 0, counter-clock's wise is +
    wdd(wdd<0)=wdd(wdd<0)+360;

  case 'Total Wind'
    vd=[];
    ud=[];
    wdd=[];
end
end

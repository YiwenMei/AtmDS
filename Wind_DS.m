% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 01/31/2019

%% Functionality
% This function downscales wind speed assuming log-wind profile with a donwscaled
% friction velocity introduced by Tao & Barrow (2017).

%% Input
% Wfn : file full name or array for coarse resolution wind speed (m/s) or a 2-by-1 
%       cell array storing in cell 1 the file full name or array for coarse resolution
%       eastward wind (m/s) and cell 2 the northward wind (m/s);
% Hmfn: file full name or array for measurement height of Wfn (m a.g.l);
% Hnfn: file full name or array for measurement height of outputted wind speed
%       (m a.g.l);
%  Zd : high resolution elevation (m a.s.l.);
% ndv : no-data value for the inputs dataset (use only one ndv for all inputs);
% z0fn: file full name or array for coarse resolution surface roughness (m);
% z0d : downscaled surface roughness (m);
% DHfn: file full name or array for coarse resolution zero-plane displacement
%       height (m).
% DHd : downscaled zero-plane displacement height (m);

%% Output
% wsd: downscaled wind speed (m/s);
% wdd: downscaled wind direction (E is 0, counter-clock's wise is +);
% Ud : downscaled eastward wind (m/s);
% Vd : downscaled northward wind (m/s).

function [wsd,wdd,Ud,Vd]=Wind_DS(Wfn,Hmfn,Hnfn,Zd,ndv,z0fn,z0d,DHfn,DHd)
%% Check the inputs
switch nargin
  case {1:4}; error('Not enough number of arguments');
  case 5; z0=1; DH=0; z0d=1; DHd=0; % Adjustment of measurement height
  case 6; error('Downscaled surface roughness missing');
  case 7; DH=0; DHd=0;
  case 8; error('Downscaled zero-plane displacement height missing');
  case 9
    if ischar(z0fn)
      z0=double(imread(z0fn));
    else
      z0=z0fn;
    end
    z0(z0==ndv)=NaN;
    if ischar(DHfn)
      DH=double(imread(DHfn));
    else
      DH=DHfn;
    end
    DH(DH==ndv)=NaN;
  otherwise; error('Too many number of arguments');
end

%% Read the wind records
if iscell(Wfn) % U and V wind
  Ufn=Wfn{1};
  Vfn=Wfn{2};
  if ischar(Ufn)
    u=double(imread(Ufn));
    v=double(imread(Vfn));
  else
    u=Ufn;
    v=Vfn;
  end
  u(u==ndv)=NaN;
  v(v==ndv)=NaN;
  ws=hypot(u,v);
  wd=atan2d(v,u); % E is 0, counter-clock's wise is +
  wd(wd<0)=wd(wd<0)+360;

elseif ischar(Wfn) % Total wind
  ws=double(imread(Wfn));
  ws(ws==ndv)=NaN;
else
  ws=Wfn;
  ws(ws==ndv)=NaN;
end

%% Parameters
if ischar(Hmfn)
  Hm=double(imread(Hmfn));
  Hm(Hm==ndv)=NaN;
else
  Hm=Hmfn;
end
if ischar(Hnfn)
  Hn=double(imread(Hnfn));
else
  Hn=Hnfn;
end
Hn(Hn==ndv)=NaN;
clear Ufn Vfn Hmfn Hnfn z0fn DHfn

%% Downscaling of wind speed/direction
k=Zd==ndv;
kus=(z0d./imresize(z0,size(Zd),'bilinear')).^.09; % Tao & Barrow (2018)
kus(k)=NaN;

kh=(Hm-DH)./z0;
kh(kh<=0)=NaN;
khi=imresize(kh,size(Zd),'bilinear');
khd=(imresize(Hn,size(Zd),'bilinear')-DHd)./z0d;
khd(khd<=0)=khi(khd<=0);
khd=log(khd)./imresize(log(kh),size(Zd),'bilinear');
khd(k)=NaN;
clear khi u v z0 DH

wsd=imresize(ws,size(Zd),'bilinear').*khd.*kus;
wsd(wsd<0)=0;

if iscell(Wfn) % U and V wind
  Vd=wsd.*imresize(sind(wd),size(Zd),'bilinear');
  Ud=wsd.*imresize(cosd(wd),size(Zd),'bilinear');
  wdd=atan2d(Vd,Ud); % E is 0, counter-clock's wise is +
  wdd(wdd<0)=wdd(wdd<0)+360;

else
  Vd=[];
  Ud=[];
  wdd=[];
end
clear ud vd khd kh kus ws wd
end

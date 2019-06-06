% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 01/31/2019

%% Functionality
% This function adjust for wind speed and direction based on terrain topography
% (Liston & Elder 2006).

%% Input
% Wfn: file full name or array for coarse resolution wind speed (m/s) or a 2-by-1 
%      cell array storing in cell 1 the file full name or array for coarse resolution
%      eastward wind (m/s) and cell 2 the northward wind (m/s);
% Asp: high resolution aspect (deg, N is 0 clock-wise is +);
% Slp: high resolution slope (deg);
% Cpl: high resolution plan curvature of terrain;
% Zd : high resolution elevation (m a.s.l.);
% ndv: no-data value for the inputs dataset (use only one ndv for all inputs).

%% Output:
% wsd: downscaled wind speed (m/s);
% wdd: downscaled wind direction (E is 0, counter-clock's wise is +);
% Ud : downscaled eastward wind (m/s);
% Vd : downscaled northward wind (m/s).

function [wsd,wdd,Ud,Vd]=Wind_Terrain(Wfn,Asp,Slp,Cpl,Zd,ndv)
%% Check the input
switch nargin
  case {1:5}; error('Not enough number of arguments');
  case 6
  otherwise; error('Too many number of arguments');
end

%% Parameters
Asp(Asp==ndv)=NaN;
Asp=90-Asp; % Convert from N is 0, closk's wise is + to E is 0, counter-clock's wise is +
Asp(Asp<0)=Asp(Asp<0)+360;
Slp(Slp==ndv)=NaN;
Cpl(Cpl==ndv)=NaN;

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
  u=imresize(u,size(Zd),'bilinear');
  v(v==ndv)=NaN;
  v=imresize(v,size(Zd),'bilinear');
  ws=hypot(u,v);
  wd=atan2d(v,u); % E is 0, counter-clock's wise is +
  wd(wd<0)=wd(wd<0)+360;

elseif ischar(Wfn) % Total wind
  ws=double(imread(Wfn));
  ws(ws==ndv)=NaN;
  ws=imresize(ws,size(Zd),'bilinear');
else
  ws=Wfn;
  ws(ws==ndv)=NaN;
  ws=imresize(ws,size(Zd),'bilinear');
end
clear Wfn Ufn Vfn u v

%% Adjust wind speed/direction to terrain properties
OSd=-sind(Slp).*cosd(wd-Asp);
OSd(Slp==0)=0;
OSd=OSd/max(abs(OSd(~isnan(OSd))));
OCd=Cpl/max(abs(Cpl(~isnan(Cpl))));
clear Cpl Slp
Wd=1+.5*(OCd+OSd); % Wind speed weight
wsd=ws.*Wd;
wsd(isnan(Wd))=ws(isnan(Wd));

if exist('wd','var')==1 % U and V wind
  thetad=rad2deg(-.5*OSd.*sind(2*(Asp-wd)));
  clear OCd OSd Asp Wd
  wdd=wd+thetad;
  wdd(isnan(thetad))=wd(isnan(thetad));
  wdd(wdd>360)=wdd(wdd>360)-360;
  wdd(wdd<0)=wdd(wdd<0)+360;

  Vd=wsd.*sind(wdd);
  Vd(wsd==0)=0; % ws=0 comes from wd is NaN
  Vd(Zd==ndv)=NaN;
  Ud=wsd.*cosd(wdd);
  Ud(wsd==0)=0;
  Ud(Zd==ndv)=NaN;

else
  wdd=[]; Ud=[]; Vd=[];
end
end

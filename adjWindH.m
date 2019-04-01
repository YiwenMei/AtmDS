% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 11/10/2018

%% Functionality
% This function calculates wind speed at specific measurement heights assuming
% either the logarithmic or the power-law wind profile.

%% Input
% Ufn : full name of files or array for the eastward or total wind record (m/s);
% Vfn : full name of files for the northward wind record (m/s, if Ufn is total
%       wind, set this to []);
% Hmfn: full name of files or array for the measurement heights of the original
%       wind speed (m);
% Hnfn: full name of files or array for the new measurement height of wind (m);
% Prof: type of wind profile assumed (can be either "log" or "pow");
% prfn: full name of files for required parameter for the logarithmic or power-law
%       wind profile:
%        - For "log", this is the surface roughness (m);
%        - For "pow", this is a coefficient dependent on the atmosphere stability
%          (1/7 is suggested for neutral stability conditions;
% dhfn: full name of files for zero-plane displacement height record used in
%       "log" (set to [] in "pow");
% ndv : no-data value for the inputs dataset (use only one ndv for all inputs);

%% Output
% v: wind speed at new measurement height (m/s);
% a: wind direction (eastward is 0 degree and counter-clock's wise is possitive);

function [v,a]=adjWindH(Ufn,Vfn,Hmfn,Hnfn,Prof,prfn,dhfn,ndv)
%% Read the inputs
if ischar(Ufn)
  u=double(imread(Ufn));
else
  u=Ufn;
end
u(u==ndv)=NaN;
if ischar(Vfn)
  v=double(imread(Vfn));
else
  v=Vfn;
end
v(v==ndv)=NaN;
if ischar(prfn)
  pr=double(imread(prfn));
else
  pr=prfn;
end
pr(pr==ndv)=NaN;
if ischar(dhfn)
  dh=double(imread(dhfn));
else
  dh=dhfn;
end
dh(dh==ndv)=NaN;
if ischar(Hmfn)
  Hm=double(imread(Hmfn));
else
  Hm=Hmfn;
end
Hm(Hm==ndv)=NaN;
if ischar(Hnfn)
  Hmn=double(imread(Hnfn));
else
  Hmn=Hnfn;
end
Hmn(Hmn==ndv)=NaN;
clear dhfn Hmfn Hnfn prfn Ufn Vfn

%% Adjust measurement height
% Wind speed and direction
if ~isempty(v)
  a=atan2d(v,u); % E is 0, counter-clock's wise is +
  a(a<0)=a(a<0)+360;
  u=hypot(u,v);
else
  a=[];
end

% Adjust
if strcmp(Prof,'log')
  del_h=Hmn-dh;
  k=del_h>0;
  v(k)=u(k).*log((del_h(k))./pr(k))./log((Hm(k)-dh(k))./pr(k));
  v(~k)=0;
elseif strcmp(Prof,'pow')
  v=u.*(Hmn./Hm).^pr;
end
v(v<0)=0;
end

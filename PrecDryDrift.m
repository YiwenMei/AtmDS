% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 05/19/2018

%% Functionality
% This function compute the "dry drift" of precipitation (Schleiss et al. 2014).
% Dry drift is the distance of a wet pixel to its closest dry pixel (He et al. 2015).

%% Input
%  p : a precipitation map or mask;
% X/Y: X/Y coordinates of grid points (m);
% pct: a percentage to define dry grid cell (set to [] if p is a precipitation
%      mask).

%% Output
% df : distance to the closest dry pixel (i.e. dry drift, m);
% Pmk: precipitation mask (0 -> dry pixel, 1 -> wet pixel).

%% Additional note
% Data must be in projected coordinate.

function [df,Pmk]=PrecDryDrift(p,X,Y,pct,oun)
X=reshape(X,size(X,1)*size(X,2),1);
Y=reshape(Y,size(Y,1)*size(Y,2),1);

if ~isempty(pct) % Define precipitation threshold
  if isempty(find(p==0, 1))
    pt=quantile(reshape(p,numel(p),1),pct);
  else
    pt=0;
  end
else
  pt=0;
end

wp=[X(p>pt) Y(p>pt)]; % Wet pixel
dp=[X(p<=pt) Y(p<=pt)]; % Dry pixel

df=zeros(size(p));
[~,ds]=knnsearch(dp,wp);
df(p>pt)=ds;

Pmk=zeros(size(p));
Pmk(p>pt)=1; % Precipitation mask
Pmk(p<=pt)=0;

% Save or not
if exist('oun','var')==1
  save(oun,'df','Pmk');
end
end

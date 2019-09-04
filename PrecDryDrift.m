% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 8/23/2019

%% Functionality
% This function compute the "dry drift" of precipitation (Schleiss et al. 2014).
%  Dry drift is the distance to its closest dry pixel of a wet pixel over the
%  mean distance of all wet pixels that are closest to the same dry pixel.

%% Input
%  p : details of the precipitation file;
% X/Y: details of the X/Y coordinates of the grids (m);
% pct: a percentage to define dry grid cell (set to [] if p is a precipitation
%      mask);
% oun: full name of file to save the dry drift and precipitation mask variable.

%% Output
% df : distance to the closest dry pixel (i.e. dry drift, m);
% Pmk: precipitation mask (0 -> dry pixel, 1 -> wet pixel).

%% Additional note
% Data must be in projected coordinate.
% Require read2Dvar.m.

function [df,Pmk]=PrecDryDrift(p,X,Y,pct,oun)
%% Define precipitation threshold
p=read2Dvar(p);
if ~isempty(pct)
  if isempty(find(p==0, 1))
    pt=quantile(reshape(p,numel(p),1),pct);
  else
    pt=0;
  end
else
  pt=0;
end

%% Relative closest distance
X=read2Dvar(X);
Y=read2Dvar(Y);
X=reshape(X,size(X,1)*size(X,2),1);
Y=reshape(Y,size(Y,1)*size(Y,2),1);
wp=[X(p>pt) Y(p>pt)]; % Wet pixel
dp=[X(p<=pt) Y(p<=pt)]; % Dry pixel

[id,ds]=knnsearch(dp,wp); % Closest distance
mds=accumarray(id,ds,[],@mean,NaN); % Mean closeset distance
N=accumarray(id,1);
mds=repelem(mds,N);
[~,I]=sort(id);
mds(I)=mds;
ds=ds./mds; % Relative closest distance

df=zeros(size(p));
df(p>pt)=ds;

%% Precipitation mask
Pmk=zeros(size(p));
Pmk(p>pt)=1;
Pmk(p<=pt)=0;

% Save or not
if exist('oun','var')==1
  save(oun,'df','Pmk');
end
end

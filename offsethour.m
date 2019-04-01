% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 1/30/2019

%% Functionality
% This function calculates the offset hours between data and solar elevation
% calculation for a day.

%% Input
% SGfn: file full name or array for coarse resolution incident shortwave radiation
%       (W/s^2);
% ndv : no-data value for the inputs dataset (use only one ndv for all inputs);
% Lat : latitude (in deg) of grids of the study domain;
% Lon : longitude (in deg) of grids of the study domain;
%  Z  : elevation of grids of the domain (m a.s.l.);
%  ds : date string of the day;
% rng : a range of offset hour (e.g. -1:.1:1).

%% Output
% dn : date number of the day;
% ofh: offset hours of the day.

function [dn,ofh]=offsethour(SGfn,ndv,Lat,Lon,Z,ds,rng)
%% Check the input
switch nargin
  case {1:6}; error('Not enough number of arguments');
  case 7
  otherwise; error('Too many number of arguments');
end

%% Find pixels with shortwave equal to 0
k=Z==ndv;
if ischar(SGfn)
  SG=double(imread(SGfn));
else
  SG=SGfn;
end
SG(SG==ndv)=1;
SG=imresize(SG,size(Z),'bilinear');
SG(k)=[];
a=length(find(SG==0));


b=length(Z(~k));
if a>0 && a<b % If a time step is partially enlighted
%% Find pixels with solar elevation lower than 0
  df=nan(length(rng),1);
  parfor i=1:length(rng)
    dn=datenum(ds,'yyyymmddHH')+rng(i)/24;
    jdn=julian(datevec(dn)); % Junlian date number
    sun=sun_positionR(jdn,[Lat(~k) Lon(~k) Z(~k)]');
    El=90-sun.zenith;
    b=length(find(El<=0));

    df(i)=abs(a-b)/a;
  end
%% Compare the amount of pixel and find the minimum difference
  I=find(df==min(df));
  [~,i]=min(abs(I-7));
  i=I(i);
  ofh=rng(i)/24;
  dn=datenum(ds,'yyyymmddHH');

else
  dn=datenum(ds,'yyyymmddHH');
  ofh=NaN;
  fprintf('%s skipped;\n',ds);
  return;
end

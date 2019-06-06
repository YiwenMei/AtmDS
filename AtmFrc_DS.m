% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 10/20/2018

%% Functionality
% Downscaling of air temperature, dew point temperature, air pressure, specific
% and relative humidity.

%% Input
% Tafn : full name of file or array for coarse resolution air temperature (K);
% LRfn : full name of file or array for air temperature lapse rate (K/m);
% Tdfn : full name of file or array for coarse resolution dew point temperature
%        (K);
% LRdfn: full name of file or array for dew point lapse rate (K/m);
% Pafn : full name of file or array for coarse resolution air pressure (Pa);
%  qfn : full name of file or array for coarse resolution specific humidity (g/g);
% RHfn : full name of file or array for coarse resolution relative humidity (%);
%   Z  : coarse resolution elevation (m);
%  Zd  : high resolution elevation (m);
%  ndv : no-data value for the inputs dataset (use only one ndv for all inputs).

%% Output
% Tad: downscaled air temperature (K);
% Tdd: downscaled dew point temperature (K);
% Pad: downscaled air pressure (Pa);
% qd : downscaled specific humidity (g/g).
% RHd: downscaled relative humidity (%).

function [Tad,Tdd,Pad,qd,RHd]=AtmFrc_DS(Tafn,LRfn,Tdfn,LRdfn,Pafn,qfn,RHfn,Z,Zd,ndv)
%% Check the inputs
switch nargin
    case {1:9}; error('Not enough number of arguments');
    case 10
    otherwise; error('Too many number of arguments');
end

%% Constant
R=287.0; % Ideal gass constant J/kg*K
g=9.81; % Gravitational acceleration m/s2
epsi=.62198; % Ratio of molecular weight of water and dry air

%% Parameters
Zd(Zd==ndv)=NaN;
Z(Z==ndv)=NaN;

if ischar(LRfn)
  LR=double(imread(LRfn));
else
  LR=LRfn;
end
if ischar(LRdfn)
  LRd=double(imread(LRdfn));
else
  LRd=LRfn;
end

%% Read the inputs
if ischar(Tafn)
  Ta=double(imread(Tafn));
else
  Ta=Tafn;
end
if ischar(Tdfn)
  Td=double(imread(Tdfn));
else
  Td=Tdfn;
end
if ischar(Pafn)
  Pa=double(imread(Pafn));
else
  Pa=Pafn;
end
if ischar(qfn)
  q=double(imread(qfn));
else
  q=qfn;
end
if ischar(RHfn)
  RH=double(imread(RHfn));
else
  RH=RHfn;
end

k=Ta==ndv | Td==ndv | Pa==ndv | q==ndv | RH==ndv;
Ta(k)=NaN;
Td(k)=NaN;
Td(Td>Ta)=Ta(Td>Ta); % Set Td > Ta to Ta
Pa(k)=NaN;
q(k)=NaN;
RH(k)=NaN;

%% Downscaling
% Air temperature (K)
dZ=Zd-imresize(Z,size(Zd),'bilinear');
Tad=imresize(Ta,size(Zd),'bilinear')+imresize(LR,size(Zd),'bilinear').*dZ;

% Air pressure (Pa)
Tm=(imresize(Ta,size(Zd),'bilinear')+Tad)/2;
Pad=imresize(Pa,size(Zd),'bilinear')./exp(g*(Zd-imresize(Z,size(Zd),'bilinear'))./(R*Tm));

% Saturated vapor pressure (Pa)
esd=Magnus_F(Tad);
es=Magnus_F(Ta);

% Dew point temperature (K)
Tdd=imresize(Td,size(Zd))+imresize(LRd,size(Zd),'bilinear').*dZ;
Tdd(Tdd>Tad)=Tad(Tdd>Tad); % Set Td > Ta to Ta

% Vapor pressure (Pa)
ed=Magnus_F(Tdd);
e=Magnus_F(Td);

% Specific Humidity
w=(Pa-(1-epsi)*e)./e; % q=epsi*e/[Pa-(1-epsi)*e];
w=imresize(w,size(Zd),'bilinear').*ed./(Pad-(1-epsi)*ed); % w=qd/q;
qd=imresize(q,size(Zd),'bilinear').*w;

% Relative Humidity
w=(es./e).*(Pa-e)./(Pa-es); % RH=e/es*(Pa-es)/(Pa-e);
w=imresize(w,size(Zd),'bilinear').*(ed./esd).*(Pad-esd)./(Pad-ed); % w=RHd/RH;
RHd=imresize(RH,size(Zd),'bilinear').*w;
RHd(RHd>100)=100;
end

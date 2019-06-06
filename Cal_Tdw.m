% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 10/10/2018

%% Functionality
% This function calculate the dew point temperature and specific/relative humidity
% if relative/specific humidity is available.

%% Input
% Tafn : full name of file store the air temperature;
% Pafn : full name of file store the air pressure;
% Humfn: full name of file store the relative/specific humidity;
% HumT : type of humidity (either "Specific" or "Relative");
%  ndv : no-data value assigned to the output images.

%% Output
% Td : dew point temperature (K);
% Hum: specific/relative humidity (g/g / %).

function [Td,Humo]=Cal_Tdw(Tafn,Pafn,Humfn,HumT,ndv)
%% Check the inputs
switch nargin
    case {1,2,3}; error('Not enough number of arguments');
    case 4; ndv=-999;
    case 5
    otherwise; error('Too many number of arguments');
end

%% Constant
epsi=.62198; % Ratio of molecular weight of water and dry air
abs0=-273.15;
% Coefficient for Magnus formula adpoted from Buck (1981)
Aw=611.21; % If Ta>5, use coeffcient of the ew2 curve in Buck (1981)
Bw=17.368;
Cw=238.88;
Am=611.21; % If -5<=Ta<=5, use coeffcient of the ew1 curve in Buck (1981)
Bm=17.502;
Cm=240.97;
Ai=611.15; % If Ta<-5, use coeffcient of the ei2 curve in Buck (1981)
Bi=22.452;
Ci=272.55;

%% Read the inputs
if ischar(Tafn)
  Ta=double(imread(Tafn));
else
  Ta=Tafn;
end
if ischar(Pafn)
  Pa=double(imread(Pafn));
else
  Pa=Pafn;
end
if ischar(Humfn)
  Hum=double(imread(Humfn));
else
  Hum=Humfn;
end
clear Tafn Pafn Humfn
k=Ta==ndv | Pa==ndv | Hum==ndv;
Ta(k)=NaN;
Pa(k)=NaN;
Hum(k)=NaN;

%% Calculate Tdw
A=Am*ones(size(Ta));
A(Ta>5)=Aw;
A(Ta<-5)=Ai;
B=Bm*ones(size(Ta));
B(Ta>5)=Bw;
B(Ta<-5)=Bi;
C=Cm*ones(size(Ta));
C(Ta>5)=Cw;
C(Ta<-5)=Ci;
es=Magnus_F(Ta);
ms=epsi*es./(Pa-es); % saturated mixing ratio

switch HumT
  case 'Specific' % If specific humidity is known
    e=Hum.*Pa./(epsi+(1-epsi)*Hum); % from q=epsi*e/(Pa-(1-epsi)*e)
    Td=C.*log(e./A)./(B-log(e./A))-abs0; % from Magnus formula
  
    m=Hum./(1-Hum); % mixing ratio
    Humo=m./ms*100; % Relative humidity
    Humo(Humo>100)=100;

  case 'Relative' % If relative humidity is known
    m=Hum.*ms; % mixing ratio
    e=m.*Pa./(m+epsi); % from m=epsi*e/(Pa-e)
    Td=C.*log(e./A)./(B-log(e./A))-abs0;

    Humo=m./(m+1); % Specific humidity
  otherwise; error('Humidity type missing');
end
end

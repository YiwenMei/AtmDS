% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 05/10/2019

%% Functionality
% Downscaling of air temperature, dew point temperature, air pressure. Specific
%  and relative humidity are calculated based on the results.

%% Input
% Ta : spatial map class (V2DCls.m) object or workspace variable for air temperature
%       (K);
% LR : V2DCls.m object or workspace variable for air temperature lapse rate (K/m);
% Td : V2DCls.m object or workspace variable for original dew point temperature (K);
% LRd: V2DCls.m object or workspace variable for dew point lapse rate (K/m);
% Pa : V2DCls.m object or workspace variable for original air pressure (Pa);
%  Z : V2DCls.m object or workspace variable for coarse resolution elevation (m);
% Zd : V2DCls.m object or workspace variable for high resolution elevation (m);

%% Output
% Tad: downscaled air temperature (K);
% Tdd: downscaled dew point temperature (K);
% Pad: downscaled air pressure (Pa);
% qd : downscaled specific humidity (g/g).
% RHd: downscaled relative humidity (%).

%% Additional note
% Require V2DCls.m and Cal_dew.m.

function [Tad,Tdd,Pad,qd,RHd]=AtmFrc_DS(Ta,LR,Td,LRd,Pa,Z,Zd)
%% Check the inputs
narginchk(7,7);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'Ta',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Ta'));
addRequired(ips,'LR',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'LR'));
addRequired(ips,'Td',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Td'));
addRequired(ips,'LRd',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'LRd'));
addRequired(ips,'Pa',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Pa'));
addRequired(ips,'Z',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Z'));
addRequired(ips,'Zd',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Zd'));
parse(ips,Ta,LR,Td,LRd,Pa,Z,Zd);
clear ips

%% Constant
R=287.0; % Ideal gass constant J/kg*K
g=9.81; % Gravitational acceleration m/s2

%% Read the inputs
Z=readCls(Z);
LR=readCls(LR);
LRd=readCls(LRd);
Ta=readCls(Ta);
Td=readCls(Td);
Pa=readCls(Pa);

k=isnan(Ta) | isnan(Td) | isnan(Pa) | isnan(Z) | isnan(LR) | isnan(LRd);
Ta(k)=NaN;
Td(k)=NaN;
Td(Td>Ta)=Ta(Td>Ta); % Set Td > Ta to Ta
Pa(k)=NaN;
Z(k)=NaN;
LR(k)=NaN;
LRd(k)=NaN;

%% Downscaling
% Air temperature (K)
Zd=readCls(Zd);

dZ=Zd-imresize(Z,size(Zd),'bilinear');
clear Zd Z
Tad=imresize(Ta,size(dZ),'bilinear')+imresize(LR,size(dZ),'bilinear').*dZ;

% Air pressure (Pa)
Tm=(imresize(Ta,size(dZ),'bilinear')+Tad)/2;
Pad=imresize(Pa,size(dZ),'bilinear')./exp(g*dZ./(R*Tm));
clear Tm Ta Pa LR

% Dew point temperature (K)
Tdd=imresize(Td,size(dZ))+imresize(LRd,size(dZ),'bilinear').*dZ;
clear Td LRd
Tdd(Tdd>Tad)=Tad(Tdd>Tad); % Set Td > Ta to Ta

% Humidity
[RHd,qd]=Cal_Tdw(Tad,Pad,'Dew Point',Tdd);
RHd(RHd>100)=100;
end

function v2d=readCls(vb)
if isa(vb,'V2DCls')
  v2d=vb.readCls;
else
  v2d=vb;
end
end

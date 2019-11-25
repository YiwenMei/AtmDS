% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 11/1/2019

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
% Require V2DCls.m, Pair_Adj.m and Hum_Cal.m.

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

%% Downscaling
Z=readCls(Z);
Zd=readCls(Zd);
dZ=Zd-imresize(Z,size(Zd),'bilinear');
clear Z Zd

% Air temperature (K)
LR=readCls(LR);
Ta=readCls(Ta);
Tad=imresize(Ta,size(dZ),'bilinear')+imresize(LR,size(dZ),'bilinear').*dZ;

% Air pressure (Pa)
Pa=readCls(Pa);
Pad=Pair_Adj(imresize(Pa,size(dZ),'bilinear'),0,imresize(Ta,size(dZ),'bilinear'),dZ,LR);
clear Pa LR

% Dew point temperature (K)
LRd=readCls(LRd);
Td=readCls(Td);
Td(Td>Ta)=Ta(Td>Ta); % Set Td > Ta to Ta
Tdd=imresize(Td,size(dZ))+imresize(LRd,size(dZ),'bilinear').*dZ;
clear Td LRd Ta
Tdd(Tdd>Tad)=Tad(Tdd>Tad); % Set Td > Ta to Ta

% Humidity
[RHd,qd]=Hum_Cal(Tad,Pad,'Dew Point',Tdd);
end

function v2d=readCls(vb)
if isa(vb,'V2DCls')
  v2d=vb.readCls;
else
  v2d=vb;
end
end

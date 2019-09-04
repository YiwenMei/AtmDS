% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 05/10/2019

%% Functionality
% Downscaling of air temperature, dew point temperature, air pressure. Specific
%  and relative humidity are calculated based on the results.

%% Input
% Ta_fd : details of file or workspace variable for original air temperature (K);
% LR_fd : details of file or workspace variable for air temperature lapse rate (K/m);
% Td_fd : details of file or workspace variable for original dew point temperature (K);
% LRd_fd: details of file or workspace variable for dew point lapse rate (K/m);
% Pa_fd : details of file or workspace variable for original air pressure (Pa);
%  Z_fd : details of file or workspace variable for coarse resolution elevation (m);
% Zd_fd : details of file or workspace variable for high resolution elevation (m);

%% Output
% Tad: downscaled air temperature (K);
% Tdd: downscaled dew point temperature (K);
% Pad: downscaled air pressure (Pa);
% qd : downscaled specific humidity (g/g).
% RHd: downscaled relative humidity (%).

%% Additional note
% Require read2Dvar.m, Magnus_F.m, and Cal_dew.m.

function [Tad,Tdd,Pad,qd,RHd]=AtmFrc_DS(Ta_fd,LR_fd,Td_fd,LRd_fd,Pa_fd,Z_fd,Zd_fd)
%% Check the inputs
narginchk(7,7);
ips=inputParser;
ips.FunctionName=mfilename;
fprintf('%s received 7 required inputs\n',mfilename);

addRequired(ips,'Ta_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Ta_fd',1));
addRequired(ips,'LR_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'LR_fd',2));
addRequired(ips,'Td_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Td_fd',3));
addRequired(ips,'LRd_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'LRd_fd',4));
addRequired(ips,'Pa_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Pa_fd',5));
addRequired(ips,'Z_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Z_fd',6));
addRequired(ips,'Zd_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Zd_fd',7));
parse(ips,Ta_fd,LR_fd,Td_fd,LRd_fd,Pa_fd,Z_fd,Zd_fd);
clear ips

%% Constant
R=287.0; % Ideal gass constant J/kg*K
g=9.81; % Gravitational acceleration m/s2

%% Read the inputs
FunName=dbstack;
Z=read2Dvar(Z_fd,FunName);
LR=read2Dvar(LR_fd,FunName);
LRd=read2Dvar(LRd_fd,FunName);
Ta=read2Dvar(Ta_fd,FunName);
Td=read2Dvar(Td_fd,FunName);
Pa=read2Dvar(Pa_fd,FunName);
clear Ta_fd Td_fd Pa_fd Z_fd LR_fd LRd_fd

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
Zd=read2Dvar(Zd_fd,FunName);
clear Zd_fd

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

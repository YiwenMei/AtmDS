% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 10/31/2019

%% Functionality
% This function adjust for wind speed and direction based on opographical properties
%  of terrain.

%% Input
% ws : spatial wind class (Wind2DCls.m) object or workspace variable for meridional
%       and zonal wind (m/s);
% Asp: spatial map class (V2DCls.m) object or workspace variable for terrain
%       aspect (deg, N is 0 clock-wise is +);
% Slp: V2DCls.m object or workspace variable for terrain slope (deg);
% Cpl: V2DCls.m object or workspace variable for plan curvature of terrain.

%% Output:
% wsd: downscaled wind speed (m/s);
% wdd: downscaled wind direction (E is 0, counter-clock's wise is +);
% Ud : downscaled eastward wind (m/s);
% Vd : downscaled northward wind (m/s).

function [wsd,wd,Ud,Vd]=Wind_Terrain(ws,Asp,Slp,Cpl)
%% Check the input
narginchk(4,4);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'ws',@(x) validateattributes(x,{'double','Wind2DCls'},{'nonempty'},mfilename,'ws'));
addRequired(ips,'Asp',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Asp'));
addRequired(ips,'Slp',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Slp'));
addRequired(ips,'Cpl',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Cpl'));

parse(ips,ws,Asp,Slp,Cpl);
clear ips

%% Read the wind records
if isa(ws,'Wind2DCls')
  [ws,wd,~,~]=ws.readCls;
  wd(wd<0)=wd(wd<0)+360; % E is 0, counter-clock's wise is +

else
  wd=atan2d(ws(:,:,2),ws(:,:,1)); % E is 0, counter-clock's wise is +
  wd(wd<0)=wd(wd<0)+360;
  ws=hypot(ws(:,:,1),ws(:,:,2));
end

%% Adjust wind speed
Slp=readCls(Slp);
Asp=readCls(Asp);
Asp=90-Asp; % Convert from N is 0, closk's wise is + to E is 0, counter-clock's wise is +
Asp(Asp<0)=Asp(Asp<0)+360;
OSd=-sind(Slp).*cosd(wd-Asp);
OSd(Slp==0)=0;
OSd=OSd/max(abs(OSd(~isnan(OSd))));
clear Slp

Cpl=readCls(Cpl);
OCd=Cpl/max(abs(Cpl(~isnan(Cpl))));
clear Cpl

Wd=1+.5*(OCd+OSd); % Wind speed weight
wsd=ws.*Wd;
wsd(isnan(Wd))=ws(isnan(Wd));
clear Wd OCd ws

%% Adjust wind direction
thetad=rad2deg(-.5*OSd.*sind(2*(Asp-wd)));
clear Asp OSd
wd(~isnan(thetad))=wd(~isnan(thetad))+thetad(~isnan(thetad));
clear thetad
wd(wd>360)=wd(wd>360)-360;
wd(wd<0)=wd(wd<0)+360;

Vd=wsd.*sind(wd);
Ud=wsd.*cosd(wd);
end

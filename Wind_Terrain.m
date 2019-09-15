% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 01/31/2019

%% Functionality
% This function adjust for wind speed and direction based on terrain topography
% (Liston & Elder 2006).

%% Input
% Wfn: file full name or array for coarse resolution wind speed (m/s) or a 2-by-1 
%      cell array storing in cell 1 the file full name or array for coarse resolution
%      eastward wind (m/s) and cell 2 the northward wind (m/s);
% Asp: high resolution aspect (deg, N is 0 clock-wise is +);
% Slp: high resolution slope (deg);
% Cpl: high resolution plan curvature of terrain;
% Zd : high resolution elevation (m a.s.l.);
% ndv: no-data value for the inputs dataset (use only one ndv for all inputs).

%% Output:
% wsd: downscaled wind speed (m/s);
% wdd: downscaled wind direction (E is 0, counter-clock's wise is +);
% Ud : downscaled eastward wind (m/s);
% Vd : downscaled northward wind (m/s).

function [wsd,wd,Ud,Vd]=Wind_Terrain(ws,wty,Asp,Slp,Cpl)
%% Check the input
narginchk(5,5);
ips=inputParser;
ips.FunctionName=mfilename;
fprintf('%s received 5 required inputs\n',mfilename);

addRequired(ips,'ws',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'ws',1));
addRequired(ips,'wty',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wty',2));
addRequired(ips,'Asp',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Asp',3));
addRequired(ips,'Slp',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Slp',4));
addRequired(ips,'Cpl',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Cpl',5));

parse(ips,ws,wty,Asp,Slp,Cpl);
clear ips

%% Read the wind records
switch wty
  case 'UV Wind' % U and V wind
    u=read2Dvar(ws{1});
    v=read2Dvar(ws{2});

    ws=hypot(u,v);
    wd=atan2d(v,u); % E is 0, counter-clock's wise is +
    wd(wd<0)=wd(wd<0)+360;
    clear u v

  case 'Total Wind' % Total wind
    ws=read2Dvar(ws);
end

%% Adjust wind speed/direction to terrain properties
Slp=read2Dvar(Slp);
Asp=read2Dvar(Asp);
Asp=90-Asp; % Convert from N is 0, closk's wise is + to E is 0, counter-clock's wise is +
Asp(Asp<0)=Asp(Asp<0)+360;
OSd=-sind(Slp).*cosd(wd-Asp);
OSd(Slp==0)=0;
OSd=OSd/max(abs(OSd(~isnan(OSd))));
clear Slp

Cpl=read2Dvar(Cpl);
OCd=Cpl/max(abs(Cpl(~isnan(Cpl))));
clear Cpl

Wd=1+.5*(OCd+OSd); % Wind speed weight
wsd=ws.*Wd;
wsd(isnan(Wd))=ws(isnan(Wd));
clear Wd OCd ws

switch wty
  case 'UV Wind' % U and V wind
    thetad=rad2deg(-.5*OSd.*sind(2*(Asp-wd)));
    clear Asp OSd
    wd(~isnan(thetad))=wd(~isnan(thetad))+thetad(~isnan(thetad));
    clear thetad
    wd(wd>360)=wd(wd>360)-360;
    wd(wd<0)=wd(wd<0)+360;

    Vd=wsd.*sind(wd);
    Ud=wsd.*cosd(wd);

  case 'Total Wind' % Total wind
    wd=[];
    Ud=[];
    Vd=[];
end
end

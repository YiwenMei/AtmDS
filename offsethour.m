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

function ofh=offsethour(SG_fd,Lat_fd,Lon_fd,Z_fd,dn,rng)
%% Check the input
narginchk(6,6);
ips=inputParser;
ips.FunctionName=mfilename;
% fprintf('%s received 6 required inputs\n',mfilename);

addRequired(ips,'SG_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'SG_fd',1));
addRequired(ips,'Lat_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Lat_fd',2));
addRequired(ips,'Lon_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Lon_fd',3));
addRequired(ips,'Z_fd',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Z_fd',4));
addRequired(ips,'dn',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'dn',5));
addRequired(ips,'rng',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'rng',6));

parse(ips,SG_fd,Lat_fd,Lon_fd,Z_fd,dn,rng);
clear ips

%% Read the inputs
Z=read2Dvar(Z_fd);
Lat=read2Dvar(Lat_fd);
Lon=read2Dvar(Lon_fd);
SG=read2Dvar(SG_fd);
SG=imresize(SG,size(Z),'bilinear');
SG(isnan(Z))=NaN;
clear SG_fd Lat_fd Lon_fd Z_fd

%% Find the inconsistent pixels
df=nan(length(rng),1);
[~,El]=Call_sunR(dn,Lat,Lon,Z);
if ~isempty(find((El>0 & SG==0) | (El<=0 & SG>0), 1))
  parfor i=1:length(rng)
    dni=dn+rng(i)/24;
    [~,El]=Call_sunR(dni,Lat,Lon,Z);

    df(i)=length(find((El>0 & SG==0) | (El<=0 & SG>0)));
  end
  clear Lat Lon Z SG

%% Compare the amount of pixel and find the minimum difference
  I=find(df==min(df));
  [~,i]=min(abs(I-length(I)));
  i=I(i);
  ofh=rng(i)/24;

else
  ofh=NaN;
  if isempty(find(El>0, 1))
    fprintf('%s all night\n',datestr(dn,'yyyymmddHH'));
  elseif isempty(find(El<=0, 1))
    fprintf('%s all day\n',datestr(dn,'yyyymmddHH'));
  else
    fprintf('%s matched\n',datestr(dn,'yyyymmddHH'));
  end
end
end

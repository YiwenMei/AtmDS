% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 04/15/2018

%% Functionality
% Calculate the shadow mask of terrain based on the horizontal angle at the
% direction of solar azimuth and solor altitude (i.e. Hrz > El -> shadowed,
% Ruiz-Arias et al. 2010). This function is designed for catersian coordinate
% only and it determines points that form the loxodrome based on the x and
% y coordinate. Thus, it is efficient computationally. 

%% Input:
% X/Y: X/Y-coordinate (in m) of grids of the study domain (N by M, where N &
%      M refer to the number of grids in lat and lon);
% Lat: Latitude (in deg) of grids of the study domain (N by M);
% rs : resolution of grids (m);
%  Z : elevation of grids of the domain (m a.s.l., N by M);
% Az : Solar azimuth for the study domain (deg, N is 0 clock's wise is +; can
%      be any dimension representing the Az of the study domain);
% El : Solar altitude for the study domain (deg, can be any dimension);
% ndv: no-data value for the inputs dataset (use only one ndv for all inputs);
% rm : maximum search radius (m, unset it to use the location specific rm).

%% Output:
% Hrz: horizontal angle of terrain at the direction of solar azimuth (deg);
% SMk: a shadow mask for terrain (0 -> shadowed, 1 -> non-shadowed).

%% Additional note:
% Please use a Catesian coordinate for the calculation;
% Bilinear interpolation (aggregation) is used to downscale (upscale) Az/El if
%  the dimensions are different than Z;

function [hrz,SMk]=CastShadow(X,Y,Lat,rs,Z,Az,El,ndv,rm)
%% Parameters
miX=min(X(1,:));
maX=max(X(1,:));
miY=min(Y(:,1));
maY=max(Y(:,1));

Az=imresize(Az,size(Z),'bilinear');
El=imresize(El,size(Z),'bilinear');
k=Z==ndv | El<0;
Az=Az(~k);
El=El(~k);
X=X(~k);
Y=Y(~k);
Lat=Lat(~k);

% Earth radius changes with latitude
r1=6378137;
r2=6356752;
Ro=sqrt(((r1^2*cosd(Lat)).^2+(r2^2*sind(Lat)).^2)./((r1*cosd(Lat)).^2+(r2*sind(Lat)).^2));

%% Maximum search radius
k=find(~k); % recycle k
Zf=Z(k);

switch nargin
  case {1:7}; error('Not enough number of arguments');
  case 8
    rm=sqrt(Zf.*(2*Ro+Zf));
    rm(Zf<0)=0;
  case 9
    rm=rm*ones(size(Zf));
  otherwise; error('Too many number of arguments');
end
N=round(rm/rs);
clear Lat Ro

%% Horizontal angle on solar azimuth direction
ins=ceil(rs*sind(45)); % Increment of searching on loxordrome
Hrz=nan(size(Zf));
for n=1:max(N)
  hrz=nan(size(Zf));
  rmi=nan(size(Zf));
  k1=N<n; % Pick out grid points with max horizon (in N grid) less than the
          % current searching horizon n

% Find coordinate of actual points with min distance to loxordrome
  xn=n*ins*sind(Az)+X; % Ideal point coordinates on loxordrome
  yn=n*ins*cosd(Az)+Y;
  xn(k1)=[];
  yn(k1)=[];

  inx=round((xn-X(~k1))/rs);  % Distance (in N grid) between ideal and actual points
  iny=round((yn-Y(~k1))/-rs);

  xn=X(~k1)+inx*rs; % Coordinate of actual points
  xn(xn<miX)=miX;
  xn(xn>maX)=maX;
  yn=Y(~k1)+iny*-rs;
  yn(yn<miY)=miY;
  yn(yn>maY)=maY;

  idxn=round(xn-miX*ones(size(Zf(~k1))))/rs+1; % Grid code of actual points
  idyn=round(yn-maY*ones(size(Zf(~k1))))/-rs+1;
  id=(idxn-1)*size(Z,1)+idyn;
  clear idxn idyn inx iny

% Horizontal angle of target grid cell to the n-th actual point
  D=hypot(X(~k1)-xn,Y(~k1)-yn); % Distance to target grid cell
  tanHrz=(Z(id)-Zf(~k1))./D; % tanget of horizon angle

  hrz1=atand(tanHrz);
  hrz(~k1)=hrz1;
  Hrz=max([Hrz hrz zeros(size(Zf))],[],2,'omitnan');
  clear hrz1 hrz D

% Reduce the maximum search times using tan(aH_0i)<=tan(aH_0j)=(Z_j-Z_0)/D_j
  tanHrz(isinf(tanHrz) | tanHrz<=0)=NaN;
  rm1=(max(Zf)-Zf(~k1))./tanHrz;
  rmi(~k1)=rm1;
  clear rm1 tanHrz
  rm=min([rmi rm],[],2,'omitnan');
  N=round(rm/rs);
  if n>=max(N,[],'omitnan')
    break
  end
end

% Reshaped horizontal angle
hrz=ndv*ones(size(Z));
hrz(Z~=ndv)=0;
hrz(k)=Hrz;

% Shadow mask
delta=heaviside(El-Hrz);
SMk=ndv*ones(size(Z));
SMk(Z~=ndv)=0;
SMk(k)=delta;
end

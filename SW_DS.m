% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 04/14/2018

%% Functionality
% Downscaling of shortwave radiation with 5 steps:
%  1)Partition of global shortwave into direct and diffuse shortwave based on
%    clear sky index (Ruiz-Arias et al. 2010);
%  2)Adjust direct shortwave for optical air depth difference, local illumination
%    and cast-shadowing (Tao & Barros 2018; Fiddes & Gruber 2014);
%  3)Adjust diffuse shortwave for sky view factor (Ruiz-Arias et al. 2012);
%  4)Calculate reflected shortwave (Tao & Barrow 2017);
%  4)Sum the direct, diffuse and reflected component up for the global shortwave.

%% Input:
% SGfn : full name of file or array for coarse resolution incident shortwave
%        at land surface (W/s^2);
% STfn : full name of file or array for coarse resolution incident shortwave
%        at top-of-atmosphere (W/s^2);
% Albfn: full name of file or array for coarse resolution surface albedo;
% Pafn : full name of file or array for coarse resolution air pressure (Pa);
%  Asp : high resolution terrain aspect (deg, N is 0 clock's wise is +);
%  Slp : high resolution terrain slope (deg);
%  SMk : a high resolution binary mask representing the shadowed area (0: shadowed,
%        1: non-shadowed);
%  SVF : high resolution sky view factor;
%  TCF : high resolution terrain configuration factor (if it is set to [], the
%        algorithm evokes TCF=(1+cosd(Slp))/2-SVF);
%  Pad : downscaled air pressure (Pa);
%   Az : Solar azimuth for the study domain (deg, N is 0 clock's wise is +);
%   El : Solar altitude for the study domain (deg);
%  ndv : no-data value for the inputs dataset (use only one ndv for all inputs);

%% Output:
% SGd: downscaled incident shortwave at land surface (W/s^2);

function [SGd,Sbd,Sdd,Srd]=SW_DS(SGfn,Albfn,Pafn,Asp,Slp,SMk,SVF,Pad,Az,El,ndv,STfn,TCF)
%% Check the input
switch nargin
  case {1:11}; error('Not enough number of arguments');
  case 12; TCF=(1+cosd(Slp))/2-SVF; % Terrain configuration factor
  case 13
  otherwise; error('Too many number of arguments');
end

%% Terrain properties
Asp(Asp==ndv)=NaN;
Slp(Slp==ndv)=NaN;
SVF(SVF==ndv)=NaN;
TCF(TCF==ndv)=NaN;
TCF(TCF<0)=0;
SMk(SMk==ndv)=NaN;
Pad(Pad==ndv)=NaN;

Az=imresize(Az,size(Pad),'bilinear');
El=imresize(El,size(Pad),'bilinear');
El(El>90)=90;
El(El<0)=0;

%% Read the forcing
if ischar(SGfn)
  SG=double(imread(SGfn));
else
  SG=SGfn;
end
if ischar(SGfn)
  ST=double(imread(STfn));
else
  ST=STfn;
end
if ischar(SGfn)
  Alb=double(imread(Albfn));
else
  Alb=ALbfn;
end
k=Alb==ndv;
Alb(k)=0; % ndv exist as place holder during the night
if ischar(SGfn)
  Pa=double(imread(Pafn));
else
  Pa=Pafn;
end
k=SG==ndv | ST==ndv | Pa==ndv;
SG(k)=NaN;
ST(k)=NaN;
Pa(k)=NaN;

%% Partition of Ss
kt=SG./ST; % Clear sky index
kd=.952-1.041*exp(-exp(2.3-4.702*kt)); % Diffuse weight
kd(ST==0)=1;

% Partition of direct and diffuse
Sb=(1-kd).*SG;
Sb=imresize(Sb,size(Pad),'bilinear');
Sd=kd.*SG;
Sd=imresize(Sd,size(Pad),'bilinear');
clear kd

%% Direct shortwave radiation
% Optical depth difference of air mass
kt=(log(SG)-log(ST))./Pa;
kt(ST==0)=0;
kt=imresize(kt,size(Pad),'bilinear');
kt=exp(kt.*(Pad-imresize(Pa,size(Pad),'bilinear')));

% Self-shadowing
zen=90-El; % Zenith angle
zen(zen>85)=85;
cosi=cosd(Slp)+sind(Slp).*tand(zen).*cosd(Az-Asp); % Illumination angle
cosi(Slp==0)=1; % Slp=0 means Asp=NaN except the ocean
cosi(cosi<0)=0;

% Cast-shadowing
SMka=imresize(SG,size(Pad),'bilinear');
SMk(SMka>0 & El==0)=1;

% Adjust for optical depth, illumination and shadow
Sbd=Sb.*cosi.*SMk.*kt;

%% Diffuse shortwave radiation
Sdd=Sd.*SVF;

%% Reflect radiation
Srd=imresize(Alb,size(Pad),'bilinear').*TCF.*(Sdd.*(1-SVF)+Sbd.*cosd(zen));

%% Global shortwave
SGd=Sbd+Sdd+Srd;
end

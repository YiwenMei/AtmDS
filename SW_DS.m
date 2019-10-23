% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 8/2/2019

%% Functionality
% Downscaling of shortwave radiation with 5 steps:
%  1)Partition of global shortwave into direct and diffuse shortwave based on
%    clear sky index (Ruiz-Arias et al. 2010);
%  2)Adjust direct shortwave for optical air depth difference, local illumination
%    and cast-shadowing (Tao & Barros 2018; Fiddes & Gruber 2014);
%  3)Adjust diffuse shortwave for sky view factor (Ruiz-Arias et al. 2012);
%  4)Calculate reflected shortwave (Tao & Barrow 2017);
%  5)Sum the direct, diffuse and reflected component up for the global shortwave.

%% Input:
% SG : spatial map class (V2DCls.m) object or workspace variable for original
%       incident shortwave (W/m2);
% InS: characters specifying methods to calculate broad-band atmospheric transmissivity
%       (this can be the the broad-band atmospheric transmissivity or incident
%       shortwave flux at TOA, W/m2. Possible strings are 'Atm Trans' or 'TOA SW');
% InV: V2DCls.m object or workspace variable for high resolution InS;
% Pa : V2DCls.m object or workspace variable for original air pressure (Pa);
% Pad: V2DCls.m object or workspace variable for downscaled air pressure (Pa);
% Asp: V2DCls.m object or workspace variable for high resolution terrain aspect
%       (deg, N is 0 clock's wise is +);
% Slp: V2DCls.m object or workspace variable for high resolution terrain slope (deg);
% Mk : V2DCls.m object or workspace variable for high resolution binary shadow
%       mask representing the shadowed area (0: shadowed, 1: non-shadowed);
% SVF: V2DCls.m object or workspace variable for high resolution sky view factor;
% Az : Solar azimuth for the study domain (deg, N is 0 clock's wise is +);
% El : Solar altitude for the study domain (deg);
% Ab : V2DCls.m object or workspace variable for original surface albedo;

% BSA: V2DCls.m object or workspace variable for the black-sky albedo;
% WSA: V2DCls.m object or workspace variable for the white-sky albedo.

%% Output:
% SGd: downscaled incident shortwave flux (W/m2);
% Sbd: downscaled incident beam shortwave flux (W/m2);
% Sdd: downscaled incident diffuse shortwave flux (W/m2);
% Srd: downscaled incident reflect shortwave flux (W/m2).

%% Additional note:
% Require V2DCls.m.

function [SGd,Sbd,Sdd,Srd]=SW_DS(SG,InS,InV,Pa,Pad,Asp,Slp,MK,SVF,Az,El,Ab,varargin)
%% Check the inputs
narginchk(12,14);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'SG',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'SG',1));
addRequired(ips,'InS',@(x) any(strcmp(x,{'Atm Trans','TOA SW'})));
addRequired(ips,'InV',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'InV',3));
addRequired(ips,'Pa',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Pa',4));
addRequired(ips,'Pad',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Pad',5));
addRequired(ips,'Asp',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Asp',6));
addRequired(ips,'Slp',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Slp',7));
addRequired(ips,'MK',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'MK',8));
addRequired(ips,'SVF',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'SVF',9));
addRequired(ips,'Az',@(x) isempty(find(x<0 | x>360, 1)));
addRequired(ips,'El',@(x) isempty(find(x>90, 1)));
addRequired(ips,'Ab',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'Ab',12));

addOptional(ips,'BSA',[],@(x) validateattributes(x,{'double','V2DCls'},{},mfilename,'BSA',13));
addOptional(ips,'WSA',[],@(x) validateattributes(x,{'double','V2DCls'},{},mfilename,'WSA',14));
parse(ips,SG,InS,InV,Pa,Pad,Asp,Slp,MK,SVF,Az,El,Ab,varargin{:});
BSA=ips.Results.BSA;
WSA=ips.Results.WSA;
clear ips varargin

%% Partition of direct and diffuse flux
SG=readCls(SG);
Pad=readCls(Pad);
if ~isempty(find(SG>0, 1))
  Pa=readCls(Pa);
  switch InS
    case 'Atm Trans'
      tao_m=readCls(InV);
    case 'TOA SW'
      ST=readCls(InV);
      if ~isempty(find(SG>0 & ST==0, 1))
        error('When SG>0, ST must >0');
      end
      tao_m=SG./ST;
    otherwise
      error('InS must be "Atm Trans" or "TOA SW"');
  end
  k=isnan(SG) | isnan(Pa);
  SG(k)=NaN;
  Pa(k)=NaN;
  tao_m(k)=NaN;

  kd=.952-1.041*exp(-exp(2.3-4.702*tao_m)); % Diffuse weight
  kd(isnan(tao_m))=.5;
  Sb=(1-kd).*SG; % Apply the weight
  Sd=kd.*SG;
  clear ST SG InV

%% Direct shortwave radiation
% Optical depth difference of air mass
  kt=(log(tao_m))./Pa; % Atmospheric attenuation coefficient
  kt(isnan(tao_m))=0;
  kt=imresize(kt,size(Pad),'bilinear');
  kt=exp(kt.*(Pad-imresize(Pa,size(Pad),'bilinear'))); % Adjustment factor 1
  clear Pad Pa tao_m

% Self-shadowing
  Asp=readCls(Asp);
  Slp=readCls(Slp);
  Az=imresize(Az,size(Asp),'bilinear');
  El=imresize(El,size(Asp),'bilinear');
  El(El<0)=0;

  cosi=cosd(Slp)+sind(Slp).*tand(min(90-El,85)).*cosd(Az-Asp); % Adjustment factor 2
% Cosine of illumination angle (Zenith angle is bounded to [0 85])
  cosi(Slp==0)=1; % when Slp is 0 (Asp is undefined), no adjustment except the ocean
  cosi(cosi<0)=0; % when cosi is 0, illumination angle is 0, Sb of the pixel is 0
  clear Az Asp

% Cast-shadowing
  MK=readCls(MK);
  a=El>0; % Fill NaN on the coast and egde of image
  MK(isnan(MK))=a(isnan(MK));
  Sbi=imresize(Sb,size(MK),'bilinear');
  MK(Sbi>0 & El==0)=1;
  clear a Sbi

% Adjust for optical depth, illumination and shadow
  Sb=imresize(Sb,size(MK),'bilinear');
  Sbd=Sb.*cosi.*MK.*kt;
  clear kt cosi Sb MK

%% Diffuse shortwave radiation
  SVF=readCls(SVF);

  Sd=imresize(Sd,size(SVF),'bilinear');
  Sdd=Sd.*SVF;
  clear Sd

%% Reflect radiation
% Downscale albedo
  Ab=readCls(Ab);
  Ab(isnan(Ab))=0;
  Ab=imresize(Ab,size(Slp),'bilinear');
  if ~isempty(BSA) && ~isempty(WSA) % Use MODIS BSA and WSA for shortwave
    BSA=readCls(BSA);
    WSA=readCls(WSA);
    kd=imresize(kd,size(BSA),'bilinear');
    Abd=BSA.*(1-kd)+WSA.*kd;
    clear BSA WSA kd

    a=imhistmatch(Abd(~isnan(Ab) & ~isnan(Abd)),Ab(~isnan(Ab) & ~isnan(Abd)));
    a(a>1)=1;
    a(a<0)=0;
    Abd(~isnan(Ab) & ~isnan(Abd))=a; % Abd=inpaint_nans(Abd,2); too time consuming
    Abd(isnan(Abd))=Ab(isnan(Abd));

  else
    Abd=imresize(Ab,size(Slp),'bilinear');
  end
  clear a Ab
  
  Abd(isnan(Abd))=0; % ndv exist as place holder during the night
  TCF=(1+cosd(Slp))/2-SVF;
  TCF=TCF-1.00001*min(TCF(TCF<0)); % TCF(TCF<0)=0;

  Srd=Abd.*TCF.*(Sdd.*(1-SVF)+Sbd.*cosd(min(90-El,85)));
  clear TCF SVF Abd El k Slp

%% Global shortwave
  SGd=Sbd+Sdd+Srd;
else
  SGd=zeros(size(Pad));
  Sbd=zeros(size(Pad));
  Sdd=zeros(size(Pad));
  Srd=zeros(size(Pad));
end
end

function v2d=readCls(vb)
if isa(vb,'V2DCls')
  v2d=vb.readCls;
else
  v2d=vb;
end
end

function [SGd,Sbd,Sdd,Srd]=SW_DS(SG,InS,InV,Pa,Pad,Asp,Slp,MK,SVF,Az,El,Ab,varargin)
%% Check the inputs
narginchk(12,14);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'SG',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'SG'));
addRequired(ips,'InS',@(x) any(strcmp(x,{'User','Built-in'})));
addRequired(ips,'InV',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'InV'));
addRequired(ips,'Pa',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Pa'));
addRequired(ips,'Pad',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Pad'));
addRequired(ips,'Asp',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Asp'));
addRequired(ips,'Slp',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Slp'));
addRequired(ips,'MK',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'MK'));
addRequired(ips,'SVF',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'SVF'));
addRequired(ips,'Az',@(x) isempty(find(x<0 | x>360, 1)));
addRequired(ips,'El',@(x) isempty(find(x>90, 1)));
addRequired(ips,'Ab',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Ab'));

addOptional(ips,'BSA',[],@(x) validateattributes(x,{'double','char'},{},mfilename,'BSA'));
addOptional(ips,'WSA',[],@(x) validateattributes(x,{'double','char'},{},mfilename,'WSA'));
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
    case 'Built-in'
      S0=1362; % W/m2
      rt_R=1-0.01672*cos(0.9856*(InV-4));
      ST=S0*rt_R^2*cosd(90-El);
      ST=imresize(ST,size(SG));
      ST(ST<0)=0;
      ST(ST<SG)=SG(ST<SG);
    case 'User'
      ST=readCls(InV);
      if ~isempty(find(SG>0 & ST==0, 1))
        error('When SG>0, ST must >0');
      end
  end
  tao_m=SG./ST;
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
if isa(vb,'char')
  v2d=matfile(vb);
  vb=cell2mat(who(v2d));
  eval(sprintf('v2d=v2d.%s;',vb));
else
  v2d=vb;
end
end

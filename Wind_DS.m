function [wsd,wdd,ud,vd]=Wind_DS(ws,Hm,Hn,varargin)
%% Check inputs
narginchk(3,7);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'ws',@(x) validateattributes(x,{'double','char',},{'nonempty'},mfilename,'ws'));
addRequired(ips,'Hm',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Hm'));
addRequired(ips,'Hn',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Hn'));

addOptional(ips,'z0',1,@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'z0'));
addOptional(ips,'z0d',1,@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'z0d'));
addOptional(ips,'d0',0,@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'d0'));
addOptional(ips,'d0d',0,@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'d0d'));

parse(ips,ws,Hm,Hn,varargin{:});
z0=ips.Results.z0;
z0d=ips.Results.z0d;
d0=ips.Results.d0;
d0d=ips.Results.d0d;
clear ips varargin

%% Read the wind records
if isa(ws,'double')
  if size(ws,3)==2
    wty='Component';
    wd=atan2d(ws(:,:,2),ws(:,:,1));
    wd(wd<0)=wd(wd<0)+360; % E is 0, counter-clock's wise is +
    ws=hypot(ws(:,:,1),ws(:,:,2));
  elseif size(ws,3)==1
    wty='Total';
  else
    error('WspTyp must be either "Component" or "Total"');
  end
else
  if size(ws,1)==2
    wty='Component';
    wd=atan2d(readCls(ws(2,:)),readCls(ws(1,:)));
    wd(wd<0)=wd(wd<0)+360; % E is 0, counter-clock's wise is +
    ws=hypot(readCls(ws(1,:)),readCls(ws(2,:)));
  elseif size(ws,1)==1
    wty='Total';
    ws=readCls(ws);
  else
    error('WspTyp must be either "Component" or "Total"');
  end
end

Hm=readCls(Hm);
Hn=readCls(Hn);

%% Downscaling of wind speed/direction
if ~isa(z0,'scalar')
% Ratio of shear velocity
  z0=readCls(z0);
  z0d=readCls(z0d);

  if isempty(find(z0<=0, 1)) || isempty(find(z0d<=0, 1))
    z0i=imresize(z0,size(z0d),'bilinear');
    kus=(z0d./z0i).^.09; % Tao & Barrow (2018)
    clear z0i

% Ratio of height
    d0=readCls(d0);
    kh=(Hm-d0)./z0;
    if isempty(find(kh<=0 & ws>0, 1))
      kh(kh<=0)=exp(1); % set to e to prevent Inf in line 95

      d0d=readCls(d0d);
      khd=(imresize(Hn,size(z0d),'bilinear')-d0d)./z0d;
      khd(khd<=0)=1; % When Hn below d0, the scale factor is 0(=ln(1))
      kh=log(khd)./imresize(log(kh),size(z0d),'bilinear');
      kh(kh<0)=0; % When Hn-d0<z0, the scale factor is 0
      kh=kh.*kus; % Combine the weighting factors
      clear khi z0 z0d d0 d0d khd kus

    else
      error('When Hm below d0, ws should be 0');
    end

  else
    error('z0 and z0d must be larger than 0');
  end

else
  kh=log(Hn)./log(Hm); % kus=1
end
wsd=imresize(ws,size(kh),'bilinear').*kh;
clear kh Hn Hm ws

%% Output wind speed/direction
switch wty
  case 'Component' % U and V wind
    vd=wsd.*imresize(sind(wd),size(wsd),'bilinear');
    ud=wsd.*imresize(cosd(wd),size(wsd),'bilinear');
    wdd=atan2d(vd,ud); % E is 0, counter-clock's wise is +
    wdd(wdd<0)=wdd(wdd<0)+360;

  case 'Total'
    vd=[];
    ud=[];
    wdd=[];
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

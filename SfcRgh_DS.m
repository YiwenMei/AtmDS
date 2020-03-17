function [z0d,d0d]=SfcRgh_DS(z0,dn,LM,LC,LuT,varargin)
%% Check inputs
narginchk(5,7);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'z0',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'z0'));
addRequired(ips,'dn',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'dn'));
addRequired(ips,'LM',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'LM'));
addRequired(ips,'LC',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'LC'));
addRequired(ips,'LuT',@(x) validateattributes(x,{'cell'},{'nrows',2,'ncols',14},mfilename,'LuT'));

addOptional(ips,'VI',1,@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'VI'));
addOptional(ips,'d0',[],@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'d0'));

parse(ips,z0,dn,LM,LC,LuT,varargin{:});
VI=ips.Results.VI;
d0=ips.Results.d0;
clear ips varargin

%% Downscaled surface roughness
ms=datestr(dn,'yyyymmm');
ms=ms(5:7);
LuTH=LuT(1,:);
LuT=LuT(2,:);
k=strcmp(ms,LuTH);

% Surface roughness pattern based on land cover
Z0d=[];
for i=1:length(LuT)
  fn=retrFn(LC,'*',num2str(i,'%02i'),'');
  lcf=readCls(fn);
  Z0d=sum(cat(3,Z0d,lcf*LuT{k}(i)),3);
end
clear lcf LuT dn

% Vegetation index-adjusted surface roughness pattern
VI=readCls(VI);
VI(VI<=0)=NaN;
z0d=Z0d.*VI;
z0d(isnan(z0d))=Z0d(isnan(z0d));
clear VI Z0d

% Apply the pattern to downscale surface roughness
z0=readCls(z0);

z0u=imresize(z0d,size(z0),'bilinear');
z0d=z0d+imresize(z0-z0u,size(z0d),'bilinear');
z0u=imresize(z0,size(z0d),'bilinear');
z0d(z0d<=0)=z0u(z0d<=0);

LM=~isnan(readCls(LM));
z0d(~LM)=NaN;
clear z0u LM

%% Find the downscaled zero-plane displacement height
if isempty(d0)
  d0d=0;
else
  d0=readCls(d0);
  if size(d0)==size(z0)
    w=d0./z0;
    d0d=z0d.*imresize(w,size(z0d),'bilinear');
    d0i=imresize(d0,size(z0d),'bilinear');
    d0d(isnan(d0d))=d0i(isnan(d0d));
    d0d(isnan(z0d))=NaN;
    clear d0i

  else
    error('size of z0 and d0 must be the same');
  end
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

function fn=retrFn(fnPn,TstrFm,ds,ys)
fn=strrep(fnPn,TstrFm,ds);
if strcmp(fileparts(fnPn),'/yyyy') % Have year path or not
  fn=strrep(fn,'yyyy',ys);
end
end

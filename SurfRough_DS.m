% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 02/13/2019

%% Functionality
% This function has two functionalities:
%  1)Downscaling of surface roughness based on vegetation index;
%  2)Calculate zero-plane displacement height based on the an empirical relationship
%    with surface roughness (Kimura et al. 1999).

%% Input
% z0 : surface roughness at coarse resolution (m);
% DH : zero-plane displacement height at coarse resolution (m);
% VId: vegetation index at downscaled resolution;
% VI : vegetation index at coarse resolution;
% ndv : no-data value for the inputs dataset (use only one ndv for all inputs);

%% Output
% z0d : surface roughness at downscaled resolution (m);
% DHd : zero-plane displacement height at downscaled resolution (m);
% c_z0: coefficient (intersect and slope) of the linear model for the surface
%       roughness vs. vegetation index relationship;
% s_z0: rmse, r2, F statistics, p value, and error variance of the multi-linear
%       model;
% c_DH: coefficient (intersect and slope) of the Kimura et al. 1999 model for
%       the surface roughness vs. zero-plane displacement height relationship;
% s_DH: rmse, r2, F statistic, p value, and error variance of the linear model
%       (Kimura et al. 1999);

function [z0d,DHd,c_z0,s_z0,c_DH,s_DH]=SurfRough_DS(z0fn,VId,VI,ndv,DHfn)
%% Check the input
switch nargin
  case {1:3}; error('Not enough number of arguments');

  case 4; DH=0;

  case 5
    if ischar(DHfn)
      DH=double(imread(DHfn));
    else
      DH=DHfn;
    end
    DH(DH==ndv)=NaN;

  otherwise; error('Too many number of arguments');
end

%% Multi-linear regression for z0 vs. LC
if ischar(z0fn)
  z0=double(imread(z0fn));
else
  z0=z0fn;
end
z0(z0==ndv)=NaN;
clear z0fn DHfn

k=~isnan(VI) & ~isnan(z0);
y=z0(k);
X=[ones(length(y),1) reshape(VI(k),length(y),1)];
[b,~,res,~,stat]=regress(log(y),X);
s_z0=[sqrt(sum(res.^2)/(length(res)-length(b))) stat]; % rmse, r2 and others
c_z0=b';

% Residual map
Res=nan(size(VI));
Res(~isnan(VI))=res;
id=find(isnan(Res));
if ~isempty(id)
  xq=fix((id-1)/size(Res,1))+1;
  yq=id-size(Res,1)*(xq-1);
  id=find(~isnan(Res)); % Other pixels
  x=fix((id-1)/size(Res,1))+1;
  y=id-size(Res,1)*(x-1);
  [id,d]=knnsearch([x y],[xq yq],'K',4);
  d=d.^2;
  d=d./repmat(sum(d,2),1,size(d,2)); % Weighting factor
  res=sum(res(id).*d,2);
  Res(isnan(Res))=res;
end

%% Find the downscaled surface roughness
X=[ones(numel(VId(:,:,1)),1) reshape(VId,size(VId,1)*size(VId,2),1)];
z0d=X*b;
z0d=reshape(z0d,size(VId,1),size(VId,2));
z0d=exp(z0d+imresize(Res,size(VId),'bilinear'));

%% Find the downscaled zero-plane displacement height
if size(DH)==size(z0)
  z0=z0(~isnan(z0));
  DH=DH(~isnan(DH));
  k=DH<quantile(DH,.1);
  z0(k)=[];
  DH(k)=[];

  [b,~,res,~,stat]=regress(log(z0),[ones(length(DH),1) log(DH)]);
  s_DH=[sqrt(sum(res.^2)/(length(res)-length(b))) stat]; % rmse, r2 and others
  c_DH=b';

  DHd=exp(log(z0d)-b(1)/b(2));

else
  DHd=DH;
  s_DH=[];
  c_DH=[];
end
end

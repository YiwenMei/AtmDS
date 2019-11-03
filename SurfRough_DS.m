% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 8/16/2019

%% Functionality
% This function is used to downscale surface roughness based on land cover and
%  vegetation index and to calculate the zero-plane displacement height.

%% Input
%  z0  : details of file or workspace variable for original surface roughness (m);
%  dn  : date number of the z0 record time step;
%  LM  : details of file or workspace variable for the high resolution land mask
%        (NaN for water and others for land);
% LCFl : details of file list of the land cover classes fraction;
%  LUT : file name of the z0 look-up-table (use ',' as separator and '.dat' as
%        file extension);

% VI : details of file or workspace variable for the high resolution vegetation
%      index of the time step;
% VIm: details of file or workspace variable for the high resolution vegetation
%      climatology for a period (e.g. mean of VI of a year);
% d0 : details of file or workspace variable for original zero-plane displacement
%      height (m).

%% Output
% z0d : downscaled surface roughness (m);
% d0d : downscaled zero-plane displacement height (m).

%% Additional note
% Require read2Dvar.m, and resizeim.m.

function [z0d,d0d]=SurfRough_DS(z0,dn,LM,LC,LUT,varargin)
%% Check inputs
narginchk(5,8);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'z0',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'z0'));
addRequired(ips,'dn',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'dn'));
addRequired(ips,'LM',@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'LM'));
addRequired(ips,'LC',@(x) validateattributes(x,{'V2DTCls'},{'nonempty'},mfilename,'LC'));
addRequired(ips,'LUT',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'LUT'));

addOptional(ips,'VI',1,@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'VI'));
addOptional(ips,'VIm',1,@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,...
    'VIm'));
addOptional(ips,'d0',[],@(x) validateattributes(x,{'double','V2DCls'},{'nonempty'},mfilename,'d0'));

parse(ips,z0,dn,LM,LC,LUT,varargin{:});
VI=ips.Results.VI;
VIm=ips.Results.VIm;
d0=ips.Results.d0;
clear ips varargin

%% Downscaled surface roughness
ms=datestr(dn,'yyyymmm');
ms=ms(5:7);
LUT=readtable(LUT,'Delimiter',',','ReadVariableNames',true);
k=strcmp(ms,LUT.Properties.VariableNames);

% Surface roughness pattern based on land cover
Z0d=[];
for i=1:size(LUT,1)
  lcf=LC.readCls(i);
  Z0d=sum(cat(3,Z0d,lcf*LUT{i,k}),3);
end
clear lcf LUT dn

% Vegetation index-adjusted surface roughness pattern
VI=readCls(VI);
VI(VI<=0)=NaN;
VIm=readCls(VIm);
z0d=Z0d.*VI./VIm;
z0d(isnan(z0d))=Z0d(isnan(z0d));
clear VI VIm Z0d

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
if isa(vb,'V2DCls')
  v2d=vb.readCls;
else
  v2d=vb;
end
end

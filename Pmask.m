% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 8/23/2019

%% Functionality
% This function compute the precipitation mask with a user-provided quantile
%  and calculate the "dry drift" of precipitation (optional).
%  Dry drift is the distance to its closest dry pixel of a wet pixel over the
%  mean distance of all wet pixels that are closest to the same dry pixel.

%% Input
%  p : details of the precipitation file;

% pct : a percentage of precipitation value where grid cell with p value less
%       than this is defined as dry grid cell;
% Dflg: a flag to indicate whether to calculate dry drift or not (calculate if
%       true and not if false);
% X/Y : details of the X/Y coordinates of the grids (the grid index is used if
%       not specified).

%% Output
% Pmk: precipitation mask (0 -> dry pixel, 1 -> wet pixel);

% df : distance to the closest dry pixel (i.e. dry drift, m).

%% Additional note
% Require read2Dvar.m.

function [Pmk,df]=Pmask(p,varargin)
%% Check inputs
narginchk(1,5);
ips=inputParser;
ips.FunctionName=mfilename;
fprintf('%s received 1 required and %d optional inputs\n',mfilename,length(varargin));

addRequired(ips,'p',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'p',1));

addOptional(ips,'pct',0,@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'pct',2));
addOptional(ips,'Dflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'Dflg',3));
addOptional(ips,'X',[],@(x) validateattributes(x,{'double'},{},mfilename,'X',4));
addOptional(ips,'Y',[],@(x) validateattributes(x,{'double'},{},mfilename,'Y',5));

parse(ips,p,varargin{:});
pct=ips.Results.pct;
Dflg=ips.Results.Dflg;
X=ips.Results.X;
Y=ips.Results.Y;
clear ips varargin

%% Precipitation mask
p=read2Dvar(p);
if pct~=0
  pt=quantile(reshape(p,numel(p),1),pct); % Define precipitation threshold
end

Pmk=nan(size(p));
Pmk(p>pt)=1;
Pmk(p<=pt)=0;
clear p

%% Relative closest distance
if Dflg
% X and Y
  if isempty(X) || isempty(Y) % Use number of grid
    X=repmat(1:size(p,2),size(p,1),1);
    Y=repmat(1:size(p,1),1,size(p,2))';
  else
    X=read2Dvar(X);
    Y=read2Dvar(Y);
  end

  X=reshape(X,size(X,1)*size(X,2),1);
  Y=reshape(Y,size(Y,1)*size(Y,2),1);
  wp=[X(Pmk==1) Y(Pmk==1)]; % Wet pixel
  dp=[X(Pmk==0) Y(Pmk==0)]; % Dry pixel

  [id,ds]=knnsearch(dp,wp); % Closest distance
  mds=accumarray(id,ds,[],@mean,NaN); % Mean closeset distance
  N=accumarray(id,1);
  mds=repelem(mds,N);
  [~,I]=sort(id);
  mds(I)=mds;
  ds=ds./mds; % Relative closest distance

  df=nan(size(p));
  df(Pmk==1)=ds;
  df(Pmk==0)=0;
end
end

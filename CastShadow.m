function [Hrz,a_h,Mk]=CastShadow(X,Y,Lat,Z,Az,El,varargin)
%% Check the inputs
narginchk(6,7);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'X',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'X'));
addRequired(ips,'Y',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Y'));
addRequired(ips,'Lat',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Lat'));
addRequired(ips,'Z',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'Z'));
addRequired(ips,'Az',@(x) isempty(find(x<0 | x>360, 1)));
addRequired(ips,'El',@(x) isempty(find(x>90, 1)));

addOptional(ips,'rm',[],@(x) validateattributes(x,{'double'},{'>',0},mfilename,'rm'));
parse(ips,X,Y,Lat,Z,Az,El,varargin{:});
rm=ips.Results.rm;
clear ips varargin

%% Initialize the outputs
Z=readCls(Z);

a_h=nan(size(Z));
a_h(~isnan(Z))=0;
Hrz=nan(size(Z));
Hrz(~isnan(Z))=0;
Mk=nan(size(Z));
Mk(~isnan(Z))=0;

%% Parameters
Az=imresize(Az,size(Z),'bilinear');
El=imresize(El,size(Z),'bilinear');
k=isnan(Z) | El<0;
if ~isempty(find(~k, 1))
  Az=Az(~k);
  El=El(~k);

  X=readCls(X);
  Y=readCls(Y);
  if length(unique(abs(diff(X,[],2))))~=1 || length(unique(abs(diff(Y,[],1))))~=1
    error('Both X and Y must be monotonic with constant increments')
  else
    if unique(abs(diff(X,[],2)))~=unique(abs(diff(Y,[],1)))
      error('Resolution of X and Y must be the same');
    else
      rs=unique(abs(diff(X,[],2)));
    end
  end
  miX=min(X(1,:));
  maX=max(X(1,:));
  miY=min(Y(:,1));
  maY=max(Y(:,1));
  X=X(~k);
  Y=Y(~k);

% Earth radius changes with latitude
  Lat=readCls(Lat);
  Lat=Lat(~k);

  r1=6378137;
  r2=6356752;
  Ro=sqrt(((r1^2*cosd(Lat)).^2+(r2^2*sind(Lat)).^2)./((r1*cosd(Lat)).^2+(r2*sind(Lat)).^2));

%% Maximum search radius
  k=find(~k); % recycle k
  Zf=Z(k);

  if isempty(rm)
    Zx=max(Zf);
    Zi=min(Zf);
    rm=sqrt((Zx-Zi).*(2*Ro+Zx-Zi));
  else
    if rm<rs
      error('Searching radius rm must be greater than the resolution of grids rs');
    else
      rm=rm*ones(size(Zf));
    end
  end
  N=round(rm./rs); % Searching radius in num. of grids
  clear Lat Ro rm

%% Horizontal angle on solar azimuth direction
  ins=ceil(rs*sind(45)); % Increment of searching on loxordrome
  a_Hrz=nan(size(Zf));
  hrz=nan(size(Zf));
  for n=1:max(N)
    a_h=nan(size(Zf));
    D=nan(size(Zf));
    a_hrz_m=nan(size(Zf));
    k1=isnan(N); % Remove grids with horizon angle greater than the maximum possible horizon angle

% Progress indication
%     if rem(n,10)==0 && n~=max(N)
%       fprintf('%.2f%%. Size reduced %.2f%%. \n',100*n/max(N),100*(1-length(find(~k1))/length(N)));
%     elseif rem(n-1,10)==0 && n~=max(N)
%       fprintf('%.2f%%',100*n/max(N));
%     elseif n==max(N)
%       fprintf('%.2f%%. Size reduced %.2f%%. \n',100*n/max(N),100*(1-length(find(~k1))/length(N)));
%     else
%       fprintf('--');
%     end

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
    d=hypot(X(~k1)-xn,Y(~k1)-yn); % Distance to target grid cell
    D(~k1)=d;
    a_hrz1=atand((Z(id)-Zf(~k1))./d); % tanget of horizon angle
    a_h(~k1)=a_hrz1;
    [a_Hrz,i]=max([a_Hrz a_h],[],2,'omitnan');
    i(a_Hrz==a_h)=2;

    hrz(i==1)=hrz(i==1);
    hrz(i==2)=D(i==2);
    clear a_hrz id xn yn i D

% Reduce the search times using the maximum possible horizon angle
    a_hrz_m1=atand((max(Zf)-Zf(~k1))./d); % Maximum possible horizon angle
    a_hrz_m(~k1)=a_hrz_m1;
    N(a_Hrz>=a_hrz_m)=NaN;
  end
  clear a_hrz_m a_hrz_m1 N Az k1 X Y d

%% Outputs
  a_h(k)=a_Hrz;
  a_h(a_h<0)=0;
  Hrz(k)=hrz;
  delta=heaviside(El-a_Hrz); % Shadow mask
  Mk(k)=delta;
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

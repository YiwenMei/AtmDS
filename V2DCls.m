classdef V2DCls
% V2DCls is a class for spatial raster image.

properties
  Fnm % Name of the file(s)
  vtp % Type of variable
  ndv % No data value for the variable
  Ulm % Physical upper limit
  Llm % Physical lower limit
  Gtg % Geographic information of the variable ([xl yt;xr yb;Rx Ry] where x/y/R is
      %  the horizontal/vertical/resolution, l/r/b/t stands for left/right/bottom/top)

  Vnm % Name of the variable
  unt % Unit of the variable
end

methods
%% Object building
  function obj=V2DCls(Fnm,vtp,ndv,Ulm,Llm,Gtg,varargin)
    narginchk(6,8);
    ips=inputParser;
    ips.FunctionName=mfilename;

    addRequired(ips,'Fnm',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'Fnm'));
    addRequired(ips,'vtp',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'vtp'));
    addRequired(ips,'ndv',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ndv'));
    addRequired(ips,'Ulm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Ulm'));
    addRequired(ips,'Llm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Llm'));
    addRequired(ips,'Gtg',@(x) validateattributes(x,{'double'},{'size',[3,2]},mfilename,'Gtg'));

    addOptional(ips,'Vnm','',@(x) validateattributes(x,{'char'},{},mfilename,'Vnm'));
    addOptional(ips,'unt','-',@(x) validateattributes(x,{'char'},{},mfilename,'unt'));

    parse(ips,Fnm,vtp,ndv,Ulm,Llm,Gtg,varargin{:});
    Vnm=ips.Results.Vnm;
    unt=ips.Results.unt;
    clear ips varargin

    obj.Fnm=Fnm;
    obj.vtp=vtp;
    obj.ndv=ndv;
    obj.Ulm=Ulm;
    obj.Llm=Llm;
    obj.Gtg=Gtg;

    obj.Vnm=Vnm;
    obj.unt=unt;
  end

%% Forcing variable reading
  function v2d=readCls(obj)
    [~,nm,fex]=fileparts(obj.Fnm);
    switch fex
      case {'.tif','tiff'} % compatable for .tiff
        v2d=double(imread(obj.Fnm));
        nm=[nm fex];

      case {'.nc4','nc'} % compatable for .nc
        v2d=double(ncread(obj.Fnm,obj.Vnm))';
        nm=[nm fex ':' obj.Vnm];

      case {'.hdf','hdf5'} % compatable for .hdf5
        v2d=double(hdfread(obj.Fnm,obj.Vnm));
        nm=[nm fex ':' obj.Vnm];

      case {'.asc','.txt'}
        v2d=double(readmatrix(obj.Fnm,'Delimiter',obj.Vnm,'NumHeaderLines',5));
        nm=[nm fex];

      case '.mat'
        v2d=matfile(obj.Fnm);
        eval(sprintf('v2d=v2d.%s;',obj.Vnm));
        nm=[nm fex ':' obj.Vnm];
    end
    v2d(v2d==obj.ndv)=NaN;

% Check the boundary
    validateattributes(v2d(~isnan(v2d)),{'double'},{'<=',obj.Ulm,'>=',obj.Llm},'',nm);
  end

%% Grids of variable
  function [X,Y]=GridCls(obj)
    X=obj.Gtg(1,1)+obj.Gtg(3,1)/2:obj.Gtg(3,1):obj.Gtg(2,1)-obj.Gtg(3,1)/2;
    Y=obj.Gtg(1,2)-obj.Gtg(3,2)/2:-obj.Gtg(3,2):obj.Gtg(2,2)+obj.Gtg(3,2)/2;
    [X,Y]=meshgrid(X,Y);
  end
end
end

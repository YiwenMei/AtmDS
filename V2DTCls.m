classdef V2DTCls
% V2DTCls is a class for a stack of spatial raster image with a time line.

properties
  Fnm % Name of the file(s)
  vtp % Type of variable
  ndv % No data value for the variable
  Ulm % Physical upper limit
  Llm % Physical lower limit
  Gtg % Geographic information of the variable ([xl yt;xr yb;Rx Ry] where x/y/R is
      %  the horizontal/vertical/resolution, l/r/b/t stands for left/right/bottom/top)
  ofs % offset to UTC in hour
  TmC % Time-window convention
  TmR % Time resolution in number of hours
  TmF % Format of the file names

  Vnm % Name of the variable
  unt % Unit of the variable
end

methods
%% Object building
  function obj=V2DTCls(Fnm,vtp,ndv,Ulm,Llm,Gtg,ofs,TmC,TmR,TmF,varargin)
    narginchk(10,12);
    ips=inputParser;
    ips.FunctionName=mfilename;

    addRequired(ips,'Fnm',@(x) validateattributes(x,{'cell'},{'nonempty'},mfilename,'Fnm'));
    addRequired(ips,'vtp',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'vtp'));
    addRequired(ips,'ndv',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ndv'));
    addRequired(ips,'Ulm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Ulm'));
    addRequired(ips,'Llm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Llm'));
    addRequired(ips,'Gtg',@(x) validateattributes(x,{'double'},{'size',[3,2]},mfilename,'Gtg'));
    addRequired(ips,'ofs',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ofs'));
    addRequired(ips,'TmC',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'TmC'));
    addRequired(ips,'TmR',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'TmR'));
    addRequired(ips,'TmF',@(x) validateattributes(x,{'cell'},{'numel',2},mfilename,'TmF'));

    addOptional(ips,'Vnm','',@(x) validateattributes(x,{'char'},{},mfilename,'Vnm'));
    addOptional(ips,'unt','-',@(x) validateattributes(x,{'char'},{},mfilename,'unt'));

    parse(ips,Fnm,vtp,ndv,Ulm,Llm,Gtg,ofs,TmC,TmR,TmF,varargin{:});
    Vnm=ips.Results.Vnm;
    unt=ips.Results.unt;
    clear ips varargin

    obj.Fnm=Fnm;
    obj.vtp=vtp;
    obj.ndv=ndv;
    obj.Ulm=Ulm;
    obj.Llm=Llm;
    obj.Gtg=Gtg;
    obj.ofs=ofs;
    obj.TmC=TmC;
    obj.TmR=TmR;
    obj.TmF=TmF;

    obj.Vnm=Vnm;
    obj.unt=unt;
  end

%% Forcing variable reading
  function v2d=readCls(obj,n)
    [~,nm,fex]=fileparts(obj.Fnm{n});
    switch fex
      case {'.tif','tiff'} % compatable for .tiff
        v2d=double(imread(obj.Fnm{n}));
        nm=[nm fex];

      case {'.nc4','nc'} % compatable for .nc
        v2d=double(ncread(obj.Fnm{n},obj.Vnm))';
        nm=[nm fex ':' obj.Vnm];

      case {'.hdf','hdf5'} % compatable for .hdf5
        v2d=double(hdfread(obj.Fnm{n},obj.Vnm));
        nm=[nm fex ':' obj.Vnm];

      case {'.asc','.txt'}
        v2d=double(readmatrix(obj.Fnm{n},'Delimiter',obj.Vnm,'NumHeaderLines',5));
        nm=[nm fex];

      case '.mat'
        v2d=matfile(obj.Fnm{n});
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

%% Extract the time line
  function Ttg=TimeCls(obj)
    [~,ds,~]=cellfun(@(X) fileparts(X),obj.Fnm,'UniformOutput',false);
    ds=cellfun(@(X) X(length(obj.TmF{1})+1:end),ds,'UniformOutput',false);
    Ttg=datenum(cell2mat(ds),obj.TmF{2});
  end
end
end

classdef Wind2DTCls
% Wind2DTCls is a class for a stack of spatial raster image of wind with a time line.

properties
  Typ % Type of wind class ('Total' or 'Component' if the total wind or UV wind
      % is supplied
  Fnm % Cell array stores the names of the input files (Each cell has one or
      %  two rows for, if Typ is 'Total', the name of the total wind file or
      %  both the names of the total wind and wind direction file; if Typ is
      %  'Component', the two rows store file names for the U and V wind)
  vtp % Type of variable
  ndv % No data value for the variables
  Ulm % Physical upper limit of the total wind
  Llm % Physical lower limit of the total wind
  Gtg % Geographic information of the variables ([xl yt;xr yb;Rx Ry] where x/y/R
      %  is the horizontal/vertical/resolution, l/r/b/t stands for left/right/bottom/top)
  ofs % offset to UTC in hour
  TmC % Time-window convention
  TmR % Time resolution in number of hours
  TmF % Format of the file names

  Vnm % Cell array stores the names of the variables (must have the same dimension
      %  as Fnm within each cell for the corresponding variables except the case
      %  of geotiff file)
  unt % Unit of the variable (if it is not specified, m/s is used)
end

methods
%% Object building
  function obj=Wind2DTCls(Typ,Fnm,ndv,Ulm,Llm,Gtg,ofs,TmC,TmR,TmF,varargin)
    narginchk(10,12);
    ips=inputParser;
    ips.FunctionName=mfilename;

    addRequired(ips,'Typ',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'Typ'));
    addRequired(ips,'Fnm',@(x) validateattributes(x,{'cell'},{'nonempty'},mfilename,'Fnm'));
    addRequired(ips,'ndv',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ndv'));
    addRequired(ips,'Ulm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Ulm'));
    addRequired(ips,'Llm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Llm'));
    addRequired(ips,'Gtg',@(x) validateattributes(x,{'double'},{'size',[3,2]},mfilename,'Gtg'));
    addRequired(ips,'ofs',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ofs'));
    addRequired(ips,'TmC',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'TmC'));
    addRequired(ips,'TmR',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'TmR'));
    addRequired(ips,'TmF',@(x) validateattributes(x,{'cell'},{'numel',2},mfilename,'TmF'));

    addOptional(ips,'Vnm','',@(x) validateattributes(x,{'cell'},{},mfilename,'Vnm'));
    addOptional(ips,'unt','m/s',@(x) validateattributes(x,{'char'},{},mfilename,'unt'));

    parse(ips,Typ,Fnm,ndv,Ulm,Llm,Gtg,ofs,TmC,TmR,TmF,varargin{:});
    Vnm=ips.Results.Vnm;
    unt=ips.Results.unt;
    clear ips varargin

    obj.Typ=Typ;
    obj.Fnm=Fnm;
    obj.vtp='Wspd';
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
  function [v1,v2]=readCls(obj,n)
    v1=[];
    for v=1:size(obj.Fnm{n},1)
      [~,~,fex]=fileparts(obj.Fnm{n}(v,:));
      switch fex
        case {'.tif','tiff'} % compatable for .tiff
          v1=cat(3,v1,double(imread(obj.Fnm{n}(v,:))));

        case {'.nc4','nc'} % compatable for .nc
          v1=cat(3,v1,double(ncread(obj.Fnm{n}(v,:),obj.Vnm{v}))');

        case {'.hdf','hdf5'} % compatable for .hdf5
          v1=cat(3,v1,double(hdfread(obj.Fnm{n}(v,:),obj.Vnm{v})));

        case {'.asc','.txt'}
          v1=cat(3,v1,double(readmatrix(obj.Fnm{n}(v,:),'Delimiter',obj.Vnm,'NumHeaderLines',5)));

        case '.mat'
          v2=matfile(obj.Fnm{n}(v,:));
          eval(sprintf('v2=v2.%s;',obj.Vnm{v}));
          v1=cat(3,v1,v2);
      end
    end
    v1(v1==obj.ndv)=NaN;

% Calculat the total wind
    switch obj.Typ
      case 'Total'
        if size(v1,3)==2
          v2=v1(:,:,2);
          v1=v1(:,:,1);
        else
          v2=[];
        end

      case 'Component'
        v2=atan2d(v1(:,:,2),v1(:,:,1)); % E is 0, counter-clock's wise is +
        v1=hypot(v1(:,:,1),v1(:,:,2));
        v2(v2<0)=v2(v2<0)+360;
    end

% Check the output
    validateattributes(v1(~isnan(v1)),{'double'},{'<=',obj.Ulm,'>=',obj.Llm},'','Total wind');
  end

%% Grids of variable
  function [X,Y]=GridCls(obj)
    X=obj.Gtg(1,1)+obj.Gtg(3,1)/2:obj.Gtg(3,1):obj.Gtg(2,1)-obj.Gtg(3,1)/2;
    Y=obj.Gtg(1,2)-obj.Gtg(3,2)/2:-obj.Gtg(3,2):obj.Gtg(2,2)+obj.Gtg(3,2)/2;
    [X,Y]=meshgrid(X,Y);
  end

%% Extract the time line
  function Ttg=TimeCls(obj)
    [~,ds,~]=cellfun(@(X) fileparts(X(1,:)),obj.Fnm,'UniformOutput',false);
    ds=cellfun(@(X) X(length(obj.TmF{1})+1:end),ds,'UniformOutput',false);
    Ttg=datenum(cell2mat(ds),obj.TmF{2});
  end
end
end

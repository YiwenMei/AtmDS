classdef VarCls
%UNTITLED2 Summary of this class goes here
%   Detailed explanation goes here

properties
  Fnm
  ndv
  Ulm
  Llm
  Vnm
end

methods
  function obj = VarCls(v1,v2,v3,v4,varargin)
%
    narginchk(4,5);
    ips=inputParser;
    ips.FunctionName=mfilename;

    addRequired(ips,'v1',@(x) validateattributes(x,{'char','cell'},{'nonempty'},mfilename,'v1'));
    addRequired(ips,'v2',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'v2'));
    addRequired(ips,'v3',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'v3'));
    addRequired(ips,'v4',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'v4'));

    addOptional(ips,'v5','',@(x) validateattributes(x,{'char'},{},mfilename,'v5'));
    parse(ips,v1,v2,v3,v4,varargin{:});
    v5=ips.Results.v5;

    obj.Fnm = v1;
    obj.ndv = v2;
    obj.Ulm = v3;
    obj.Llm = v4;
    obj.Vnm = v5;
  end

%% Forcing variable reading
  function v2d = read2Dcls(obj)
%METHOD1 Summary of this method goes here
%   Detailed explanation goes here
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
end
end


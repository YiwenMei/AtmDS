classdef TSCls
% TSCls is a class for raster image or a stack of raster image. It contains
%  seven properties and five methods.

properties
  vtp % Type of variable
  Ulm % Physical upper limit
  Llm % Physical lower limit
  Gtg % Geographic information of the time series (Gtg.ID Gtg.X Gtg.Y where ID
      %  is the name of the location and X/Y is the horizontal/vertical coordinate)
  gid % Group ID of the file
  os1 % offset in hours to UTC
  TS1 % Values of the variable
  TL1 % Time line of the variable
  TR1 % Time resolution
  TC1 % Time-window convention
  os2 % offset in hours to UTC
  TS2 % Values of the variable
  TL2 % Time line of the variable
  TR2 % Time resolution
  TC2 % Time-window convention

  unt % Unit of the variable
  OFn % Output name of the file
end

methods
%% Object building
  function obj=TSCls(vtp,Ulm,Llm,Gtg,gid,os1,TS1,TL1,TR1,TC1,os2,TS2,TL2,TR2,TC2,varargin)
% Check inputs
    narginchk(15,17);
    ips=inputParser;
    ips.FunctionName=mfilename;

    addRequired(ips,'vtp',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'vtp'));
    addRequired(ips,'Ulm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Ulm'));
    addRequired(ips,'Llm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Llm'));
    addRequired(ips,'Gtg',@(x) validateattributes(x,{'struct'},{'nonempty'},mfilename,'Gtg'));
    addRequired(ips,'gid',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'gid'));
    addRequired(ips,'os1',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'os1'));
    addRequired(ips,'TS1',@(x) validateattributes(x,{'double'},{'vector'},mfilename,'TS1'));
    addRequired(ips,'TL1',@(x) validateattributes(x,{'double'},{'vector'},mfilename,'TL1'));
    addRequired(ips,'TR1',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'TR1'));
    addRequired(ips,'TC1',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'TC1'));
    addRequired(ips,'os2',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'os2'));
    addRequired(ips,'TS2',@(x) validateattributes(x,{'double'},{'vector'},mfilename,'TS2'));
    addRequired(ips,'TL2',@(x) validateattributes(x,{'double'},{'vector'},mfilename,'TL2'));
    addRequired(ips,'TR2',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'TR2'));
    addRequired(ips,'TC2',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'TC2'));

    addOptional(ips,'unt','-',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'unt'));
    addOptional(ips,'OFn','',@(x) validateattributes(x,{'cell'},{},mfilename,'OFn'));

    parse(ips,vtp,Ulm,Llm,Gtg,gid,os1,TS1,TL1,TR1,TC1,os2,TS2,TL2,TR2,TC2,varargin{:});
    unt=ips.Results.unt;
    OFn=ips.Results.OFn;
    clear ips varargin

% Assign values
    obj.vtp=vtp;
    obj.Ulm=Ulm;
    obj.Llm=Llm;
    obj.Gtg=Gtg;
    obj.gid=gid;
    obj.os1=os1;
    obj.TS1=TS1;
    obj.TL1=TL1;
    obj.TR1=TR1;
    obj.TC1=TC1;
    obj.os2=os2;
    obj.TS2=TS2;
    obj.TL2=TL2;
    obj.TR2=TR2;
    obj.TC2=TC2;

    obj.unt=unt;
    obj.OFn=OFn;

% Check the boundary
    ts=[obj.TS1;obj.TS2];
    validateattributes(ts,{'double'},{'<=',obj.Ulm,'>=',obj.Llm,'nonnan'});
  end

%% Unify the time line
  function [TL1,TL2]=UniTL(obj,CTp)
    switch CTp
      case 'UTC'
        TL1=obj.TL1+obj.os1;
        TL2=obj.TL2+obj.os2;
      case 'tg'
        TL1=obj.TL1;
        TL2=obj.TL2+obj.os2-obj.os1;
      case 'rf'
        TL1=obj.TL1+obj.os1-obj.os2;
        TL2=obj.TL2;
    end
  end
end
end

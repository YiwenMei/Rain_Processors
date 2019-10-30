% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 10/7/2019

%% Functionality
% This function is used to process the near-real-time or post-real-time high
%  resolution CMORPH satellite precipitation product (Joyce et al. 2004). Its
%  main functionalities are
%   1)unzip the CMORPH record (e.g. CMORPH_V1.0_RAW_8km-30min_2008070100.bz2);
%   2)read the unzipped record; and
%   3)extract the record based on an given lat/lon box.
%  Its optional functionalities are
%   4)resample the extracted record;
%   5)project the extracted record to another coordinate system; and
%   6)crop the projected record to a retangular box.
%  If any of the optional functionality is evoked, the processed record will
%  be outputted as geotif in the output directory.

%% Input
% fname: full name of the input high resolution CMORPH record (e.g.,
%        G:\CMORPH\HRC\200807\CMORPH_V1.0_RAW_8km-30min_2008070100.bz2;
%        G:\CMORPH\CHRC\200807\CMORPH_V1.0_ADJ_8km-30min_2008070100.bz2;);
% wkpth: working directory of the code;
%  xl  : longitude (in the range of [-180 180]) of the west boundary (xl can
%        have an optional second element to represent the west boundary coordinate
%        in the unit of the output coordinate system);
%  xr  : similar to xl but for the east boundary;
%  yb  : latitude (in the range of [-60 60]) of the south boundary (yb can have
%        an optional second element to represent the south boundary coordinate
%        in the unit of the output coordinate system);
%  yt  : similar to yb but for the north boundary;

% opth: output directory to store the outputed .tif files;
% onm : a user-assigned name for the outputted CMORPH files as character;
% ors : output coordinate system (e.g. EPSG:102009);
%  rs : x and y resolution of the outputted precipitation.

%% Output
%  P : processed precipitation with lat/lon grids (if any of the optional functionalities
%       was evoked, the lat/lon grids are implicitly included in the outputted
%       .tif files and are not included in P; the orientation of P follows the
%       human reading convention);

% ofn: name of the outputted file with the form of onmyyyymmddhh.tif stored under
%       opth (if no file is outputted, ofn is []);

% onmyyyymmddhh.tif - geotiff file stores the processed precipitation in opth
%  (it only appears if any of the optional functionalty is evoked).

%% Additional note
% 1)The input files are unzipped from the original monthly ".tar" file;
% 2)If any optional functionality is required, please make sure to have GDAL
%   installed;
% 3)The no-data value of HRC/CHRC are preseved in the .tif file; and
% 4)Require matV2tif.m.

function [P,ofn]=HRC_process(fname,wkpth,xl,xr,yb,yt,varargin)
%% Check the inputs
narginchk(6,10);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'fname',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'fname'));
addRequired(ips,'wkpth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wkpth'));
msg=sprintf('Size of xl, xr, yb, yt must be 1 or 2');
addRequired(ips,'xl',@(x) assert(~isempty(x) & length(x)<3,msg));
addRequired(ips,'xr',@(x) assert(~isempty(x) & length(x)<3,msg));
addRequired(ips,'yb',@(x) assert(~isempty(x) & length(x)<3,msg));
addRequired(ips,'yt',@(x) assert(~isempty(x) & length(x)<3,msg));

addOptional(ips,'opth','',@(x) validateattributes(x,{'char'},{},mfilename,'opth'));
addOptional(ips,'onm','CMORPH-',@(x) validateattributes(x,{'char'},{'nonempty'},...
    mfilename,'onm'));
addOptional(ips,'ors','wgs84',@(x) validateattributes(x,{'char'},{},mfilename,'ors'));
addOptional(ips,'rs',[],@(x) validateattributes(x,{'double'},{},mfilename,'rs'));

parse(ips,fname,wkpth,xl,xr,yb,yt,varargin{:});
opth=ips.Results.opth;
onm=ips.Results.onm;
ors=ips.Results.ors;
rs=ips.Results.rs;
clear ips msg varargin

%% Lat/lon grids and other info of CMORPH
rs_lon=360/4948;
rs_lat=120/1649;
Lon=0:rs_lon:360;
Lat=60:-rs_lat:-60;

ndv=-999; % no-data value of CMORPH

%% unzip and read the file
[~,nm,~]=fileparts(fname);
ufn=fullfile(wkpth,nm);
system(sprintf('bunzip2 -c "%s" > "%s"',fname,ufn));
% system(sprintf('7z e "%s" -o"%s" * -r',fname,wkpth));

fid=fopen(ufn);
p=fread(fid,'float32','l');
fclose(fid);

p=reshape(p,4948,1649,length(p)/4948/1649);
p(p==ndv)=NaN;
p=nanmean(p,3);
p=rot90(p);
p(isnan(p))=ndv;
delete(ufn);

%% Extracting the file
% Index of interested domain
if xl(1)<0 % Convert longitude to the range of [0 360];
  xl(1)=xl(1)+360;
end
if xr(1)<0
  xr(1)=xr(1)+360;
end

cl=find(xl(1)-Lon>=0,1,'last'); % left column
cr=find(xr(1)-Lon<=0,1,'first')-1; % right column
rt=find(yt(1)-Lat<=0,1,'last'); % top row
rb=find(yb(1)-Lat>=0,1,'first')-1; % bottom row
p=p(rt:rb,cl:cr); % extract
lon=Lon(cl:cr);
lat=Lat(rt:-1:rb)';
P.lat=lat;
P.lon=lon;

%% Resample/project/crop the file
pr2=[];
ofn=[];

xll=(cl-1)*rs_lon; % longitude of lower left corner
yll=60-rb*rs_lat; % latitude of lower left corner

ds=cell2mat(regexp(nm,'_(\d{10})','tokens','once'));
if ~isempty(opth)
  if system('gdalinfo --version')~=0
    error('GDAL is not detected. Please install GDAL to evoke the optional functionalities.\n');
  end

  tfn=fullfile(wkpth,sprintf('%s%s.tif',onm,ds));
  matV2tif(tfn,p,xll,yll,rs_lat,ndv,'wgs84',wkpth);

  fun='gdalwarp -overwrite -of GTiff -r bilinear -q'; % GDAL function
  pr1=sprintf('-t_srs %s',ors); % Project
  if ~isempty(rs) % Resample
    pr2=sprintf('-tr %i %i',rs(1),rs(2));
  end
  if length(xl)==2 && length(xr)==2 && length(yt)==2 && length(yb)==2
    pr3=sprintf('-te %i %i %i %i',xl(2),yb(2),xr(2),yt(2));
  elseif length(xl)==1 && length(xr)==1 && length(yt)==1 && length(yb)==1
    pr3=[]; % Crop projected image
  else
    error('sizes of xl, yb, xr, and yt must be the same and equal to 1 or 2.');
  end

  ofn=fullfile(opth,sprintf('%s%s.tif',onm,ds));
  system(sprintf('%s %s %s %s "%s" "%s"',fun,pr1,pr2,pr3,tfn,ofn));
  delete(tfn);

  P=double(imread(ofn));
end
p(p==ndv)=NaN;
P.p=p;
end

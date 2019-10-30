% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 10/11/2019

%% Functionality
% This function is used to process the CHIRP/CHIRPS precipitation product (Funk
%  et al. 2015). Its main functionalities are
%   1)read one time step (one day) within the yearly record and
%   2)extract the precipitaiton record of the time step based on a given lat/lon
%      box.
%  The extracted records will be outputted as .mat files to the output directory.
%  The optional functionalities of this function are
%   3)resample the extracted record;
%   4)project the extracted record to another coordinate system;
%   5)crop the projected record to a retangular box.
%  The operations are looped for all time steps contain in the yearly record
%   and the processed files are outputted as .tif in the output directory.

%% Input
% fname: full name with path of the input CHIRP or CHIRPS yearly record (e.g.
%         G:\CHP\CHIRP\chirp.2007.days_p05.nc; G:\CHP\CHIRPS\chirps-v2.0.2007.days_p05.nc);
% wkpth: working directory of the code;
% opth : output directory to store the outputed .mat/.tif files;
%  onm : a user-assigned name for the outputted CHIRP/CHIRPS files as character;
%  xl  : longitude (in the range of [-180 180]) of the west boundary (xl can
%        have an optional second element to represent the west boundary coordinate
%        in the unit of the output coordinate system);
%  xr  : similar to xl but for the east boundary;
%  yb  : latitude (in the range of [-50 50]) of the south boundary (yb can have
%        an optional second element to represent the south boundary coordinate
%        in the unit of the output coordinate system);
%  yt  : similar to yb but for the north boundary;

% pflg: parallel flag (false - default, squential; true - parallel);
% ors : output coordinate system (e.g. EPSG:102009);
%  rs : x and y resolution of the outputted precipitation.

%% Output
% Ofn: name list of the output .mat or .tif files stores in opth;

% onmyyyymmddhh.mat / onmyyyymmddhh.tif - .mat or .tif files of the processed
%  precipitation in opth;
% Gridyyyy.mat - the lat/lon grids of the .mat file.

%% Additional note
% 1)If any optional functionality is required, please make sure to have GDAL
%   installed;
% 2)The no-data value of CHIRP/CHIRPS are preseved in the .tif files; and
% 3)Require matV2tif.m and doy2date.m.

function Ofn=CHP_process(fname,wkpth,opth,onm,xl,xr,yb,yt,varargin)
%% Check the inputs
narginchk(8,11);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'fname',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'fname'));
addRequired(ips,'wkpth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wkpth'));
addRequired(ips,'opth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'opth'));
addRequired(ips,'onm',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'onm'));
msg=sprintf('Size of xl, xr, yb, yt must be 1 or 2');
addRequired(ips,'xl',@(x) assert(~isempty(x) & length(x)<3,msg));
addRequired(ips,'xr',@(x) assert(~isempty(x) & length(x)<3,msg));
addRequired(ips,'yb',@(x) assert(~isempty(x) & length(x)<3,msg));
addRequired(ips,'yt',@(x) assert(~isempty(x) & length(x)<3,msg));

addOptional(ips,'ors','wgs84',@(x) validateattributes(x,{'char'},{},mfilename,'ors'));
addOptional(ips,'rs',[],@(x) validateattributes(x,{'double'},{},mfilename,'rs'));
addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},...
    mfilename,'pflg'));

parse(ips,fname,wkpth,opth,onm,xl,xr,yb,yt,varargin{:});
ors=ips.Results.ors;
rs=ips.Results.rs;
pflg=ips.Results.pflg;
clear ips msg varargin

%% Lat/lon grids and other info of CHIRP/CHIRPS
Lat=double(ncread(fname,'latitude'));
Lon=double(ncread(fname,'longitude'));
rs_lon=360/length(Lon);
rs_lat=100/length(Lat);
Lon=0:rs_lon:360;
Lat=50:-rs_lat:-50;

[~,ys,~]=fileparts(fname);
ys=cell2mat(regexp(ys,'.(\d{4}).','tokens','once'));

T=double(ncread(fname,'time'));
ndv=ncinfo(fname,'precip');
ndv=double(ndv.FillValue); % no-data value

%% Index of interested domain
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

Grid.lat=Lat(rt:-1:rb)';
Grid.lon=Lon(cl:cr);
save(fullfile(opth,['Grid' ys '.mat']),'Grid');

xll=(cl-1)*rs_lon; % longitude of lower left corner
yll=50-rb*rs_lat; % latitude of lower left corner

%% Parameters for resample/project/crop the file
pr1=sprintf('-t_srs %s',ors); % Project
pr2=[];
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

%% Read, crop, and resample/project/crop the record
Ofn={};
switch pflg
  case true
    parfor d=1:length(T)
      ofn=CHP_process_sub(fname,[1,1,d],[length(Lon)-1 length(Lat)-1 1],rt,rb,cl,cr,...
          ndv,d,ys,opth,onm,pr1,pr2,pr3,xll,yll,rs_lat,wkpth);
      Ofn=[Ofn;{ofn}];
    end

  case false
    for d=1:length(T)
      ofn=CHP_process_sub(fname,[1,1,d],[length(Lon)-1 length(Lat)-1 1],rt,rb,cl,cr,...
          ndv,d,ys,opth,onm,pr1,pr2,pr3,xll,yll,rs_lat,wkpth);
      Ofn=[Ofn;{ofn}];
    end
end
end

function ofn=CHP_process_sub(fname,sid,cts,rt,rb,cl,cr,ndv,d,ys,opth,onm,pr1,pr2,pr3,...
    xll,yll,rso,wkpth)
p=rot90(ncread(fname,'precip',sid,cts));
p=[p(:,3601:7200) p(:,1:3600)]; % convert from [-180 180] to [0 360]
p=p(rt:rb,cl:cr); % crop

ds=datestr(doy2date(d,str2double(ys)),'yyyymmdd');
nm=sprintf('%s%s',onm,ds);
ofn=fullfile(opth,[nm '.mat']);
save(ofn,'p');

if ~contains(pr1,'wgs84') || ~isempty(pr2) || ~isempty(pr3)
  if system('gdalinfo --version')~=0
    error('GDAL is not detected. Please install GDAL to evoke the optional functionalities.\n');
  end
  delete(fullfile(opth,['Grid' ys '.mat'])); % Delete the Grid<ys>.mat

  tfn=fullfile(wkpth,[nm '.tif']);
  p(isnan(p))=ndv;
  matV2tif(tfn,p,xll,yll,rso,ndv,'wgs84',wkpth);
  delete(ofn); % delete the .mat file

  fun='gdalwarp -overwrite -of GTiff -r bilinear -q'; % GDAL function
  ofn=fullfile(opth,[nm '.tif']);
  system(sprintf('%s %s %s %s "%s" "%s"',fun,pr1,pr2,pr3,tfn,ofn)); % resample/project/crop
  delete(tfn);
end
end

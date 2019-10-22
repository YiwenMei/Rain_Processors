% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 10/11/2019

%% Functionality
% This function is used to process the CHIRP or CHIRPS precipitation product (Funk
%  et al. 2015). Its main functionalities are
%   1)read one time step (one day) within the yearly record (e.g. chirp.2007.days_p05.nc);
%   2)extract the precipitaiton record of the time step based on a given lat/lon
%      box;
%  Its optional functionalities are
%   3)resample the extracted record;
%   4)project the extracted record to another coordinate system;
%   5)crop the projected record to a retangular box.
%  The operations are looped for all time steps contain in the yearly record.

%% Input
% fname: full name with path of the input CHIRP or CHIRPS yearly record (e.g.
%         G:\CHP\CHIRP\chirp.2007.days_p05.nc; G:\CHP\CHIRPS\chirps-v2.0.2007.days_p05.nc);
%  cty : type of record as character (CHIRP or CHIRPS);
% wkpth: working directory of the code;
% opth : path to store the outputed .tif files;
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

% CHIRPyyyymmddhh.mat/CHIRPyyyymmddhh.tif/CHIRPSyyyymmddhh.mat/CHIRPSyyyymmddhh.tif:
%  .mat or .tif files stores the processed precipitation in opth.

%% Additional note
% 1)Please make sure to have GDAL installed and the out_path set if you want
%   to reproject the data into other coordinate system.
% 2)If reprojection is not required but record outputted as .asc is wanted, outpth
%   is required to set.
% 3)The no-data value of HRC/CHRC are preseved in the .asc and .tif record.

function Ofn=CHP_process(fname,cty,wkpth,opth,xl,xr,yb,yt,varargin)
%% Check the inputs
narginchk(8,11);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'fname',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'fname'));
expInS={'CHIRP','CHIRPS'};
msg=cell2mat(cellfun(@(x) [x ', '],expInS,'UniformOutput',false));
msg=sprintf('Expected InS to be one of the following %s\n',msg);
addRequired(ips,'cty',@(x) assert(any(strcmp(x,expInS)),msg));
addRequired(ips,'wkpth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wkpth'));
addRequried(ips,'opth',[],@(x) validateattributes(x,{'char'},{},mfilename,'opth'));
msg=sprintf('Size of xl, xr, yb, yt must be 1 or 2');
addRequired(ips,'xl',@(x) assert(~isempty(x) & length(x)<3,msg));
addRequired(ips,'xr',@(x) assert(~isempty(x) & length(x)<3,msg));
addRequired(ips,'yb',@(x) assert(~isempty(x) & length(x)<3,msg));
addRequired(ips,'yt',@(x) assert(~isempty(x) & length(x)<3,msg));

addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));
addOptional(ips,'ors','wgs84',@(x) validateattributes(x,{'char'},{},mfilename,'ors'));
addOptional(ips,'rs',[],@(x) validateattributes(x,{'double'},{},mfilename,'rs'));

parse(ips,fname,cty,wkpth,xl,xr,yb,yt,varargin{:});
pflg=ips.Results.pflg;
ors=ips.Results.ors;
rs=ips.Results.rs;
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
ndv=ndv.FillValue; % no-data value

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

xll=(cl-1)*rs_lon; % longitude of lower left corner
yll=50-rb*rs_lat; % latitude of lower left corner

%% Parameters for resample/project/crop the file
pr1=sprintf('-t_srs %s ',ors); % Project
pr2=[];
if ~isempty(rs) % Resample
  pr2=sprintf('-tr %i %i ',rs(1),rs(2));
end
pr3=[]; % Crop projected image
if length(xl)==2
  pr3=sprintf('-te %i %i %i %i ',xl(2),yb(2),xr(2),yt(2));
end

%% Read, crop, and resample/project/crop the record
Ofn={};
switch pflg
  case true
    parfor d=1:length(T)
      ofn=CHP_process_sub(fname,[1,1,d],[length(Lon)-1 length(Lat)-1 1],rt,rb,cl,cr,...
          ndv,d,ys,opth,cty,pr1,pr2,pr3,xll,yll,rs_lat,wkpth);
      Ofn=[Ofn;{ofn}];
    end

  case false
    for d=1:length(T)
      ofn=CHP_process_sub(fname,[1,1,d],[length(Lon)-1 length(Lat)-1 1],rt,rb,cl,cr,...
          ndv,d,ys,opth,cty,pr1,pr2,pr3,xll,yll,rs_lat,wkpth);
      Ofn=[Ofn;{ofn}];
    end
end
end

function ofn=CHP_process_sub(fname,sid,cts,rt,rb,cl,cr,ndv,d,ys,opth,cty,pr1,pr2,pr3,...
    xll,yll,rso,wkpth)
p=rot90(ncread(fname,'precip',sid,cts));
p=[p(:,3601:7200) p(:,1:3600)]; % convert from [-180 180] to [0 360]
p=p(rt:rb,cl:cr); % crop
p(isnan(p))=ndv;

ds=datestr(doy2date(d,str2double(ys)),'yyyymmdd');
nm=fullfile(opth,sprintf('%s%s',cty,ds));
ofn=[nm '.mat'];
save(ofn,'p')

if strcmp(pr1,'wgs84') || ~isempty(pr2) || ~isempty(pr3)
  if system('gdalinfo --version')~=0
    error('GDAL is not detected. Please install GDAL to evoke the optional functionalities.\n');
  end

  nm=[nm '.tif'];
  tfn=fullfile(wkpth,nm);
  matV2tif(tfn,p,xll,yll,rso,ndv,'wgs84',wkpth);

  fun='gdalwarp -overwrite -of GTiff -r bilinear '; % GDAL function
  ofn=fullfile(opth,nm);
  par=[pr1 pr2 pr3];
  system([fun par tfn ' ' ofn]); % resample/project/crop
  delete(tfn);
end
end

% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 6/23/2018

%% Functionality
% This function is used to process the GSMaP-MVK and the gauge-adjusted GSMaP
% satellite precipitation product (Ushio et al. 2009; Mega et al. 2014). It
% includes several functionalities:
%  1)unzip the GSMaP record (e.g., gsmap_mvk.20000301.0000.v5.222.1.dat.gz);
%  2)read and crop the record based on an given lat/lon box;
%  3)output the cropped record as .asc file (optional);
%  4)project the record in .asc to another projection and output the record as
%    GTiff (optional).

%% Input
% infname: full name with path of the input high resolution CMORPH record (e.g.,
%          G:\GSMaP\V5MVK\hourly\2000\03\gsmap_mvk.20000301.0000.v5.222.1.dat.gz;
%          G:\GSMaP\V5G\2000\03\gsmap_gauge.20000302.0000.v5.222.1.40.dat.gz;);
%   cty  : type of GSMaP record (0 for MVK, 1 for gauged);
% workpth: path to store the unzipped record;
%  xl/xr : left/right longitude of the boundary (xl/xr can have either 1 or 2 number(s)
%          where the first one represents the longitude and must be in the range of
%          [-180 180]; the second one is boundary in the projected coordinate unit);
%  yb/yt : bottom/top latitude of the boundary (yb/yt can have either 1 or 2 number(s)
%          where the first one represents the latitude and must be in the range of
%          [60 -60]; the second one is boundary in the projected coordinate unit);
% outpth : path to store the .asc and, if have, the .tif files (set it to "[]" if
%          no need to output record in .asc format);
% out_pj : output coordinate system (e.g., EPSG:102009; set it to "[]" if no
%          reprojection is required);
%   rs   : x and y resolution of the projected image (set it to "[]" if no
%          reprojection is required).

%% Output
%      p  /  p1      : cropped precipitation map in original/new projection
%                      (the orientation follows the human reading convention);
%  MGSMyyyymmddhh.asc: precipitaiton map in original projection and resolution
%                      outputted to outpth as .asc file;
% MGSMpyyyymmddhh.tif: precipitaiton map in new projection, resolution and extend
%                      outputted to outpth as .tif file;

%% Additional note
% 1)Please make sure to have GDAL installed and the out_path set if you want to
%   reproject the data into other coordinate system.
% 2)If reprojection is not required (i.e., out_pj and rs are "[]") but record
%   outputted as .asc is wanted, out_path is required to set.
% 3)The no-data values of MGSM/CGSM are preseved in the .asc and .tif record.

function [p,p1]=GSM_process(infname,cty,workpth,xl,xr,yb,yt,outpth,out_pj,rs)
% Lat/lon grids and other info of GSMaP
rs_lon=360/3600;
rs_lat=120/1200;
Lon=0:rs_lon:360;
Lat=60:-rs_lat:-60;

ndv=-99; % no-data value

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

nr=length(rt:rb); % number of row
nc=length(cl:cr); % number of column
xll=(cl-1)*rs_lon; % longitude of lower left corner
yll=60-rb*rs_lat; % latitude of lower left corner

% unzip and read the input
% system(sprintf('7z e "%s" -o"%s" * -r',infname,workpth));
gunzip(infname,workpth); % unzip

[~,nm,~]=fileparts(infname);
uz_fn=[workpth nm];
fid=fopen(uz_fn);
p=fread(fid,'float32','l');
fclose(fid);

p=reshape(p,3600,1200,length(p)/3600/1200)';

p=p(rt:rb,cl:cr); % crop
delete(uz_fn);

p1=[];
if ~isempty(out_pj)
% Creat the asc files
  ds=nm([11:11+7 20:21]); % pay careful attention to this.
  if cty==0
    name=[outpth 'MGSM' ds '.asc'];
  else
    name=[outpth 'CGSM' ds '.asc'];
  end

  fid=fopen(name,'w');
  fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',['ncols ' num2str(nc)],['nrows '...
    num2str(nr)],['xllcorner ' num2str(xll,8)],['yllcorner ' num2str(yll,8)],...
    ['cellsize ' num2str(rs_lat)],['NODATA_value ' num2str(ndv)]);
  fclose(fid);
  dlmwrite(name,p,'delimiter',' ','-append'); % output .asc

% Project to a new coordinate
  fun='gdalwarp -overwrite -of GTiff -r bilinear '; % GDAL function
  pr1='-s_srs wgs84 '; % Projection of original record
  pr2=['-t_srs ' out_pj ' '];
  pr3=[];
  if ~isempty(rs)
    pr3=sprintf('-tr %i %i ',rs(1),rs(2));
  end
  pr4=[];
  if length(xl)==2
    pr4=sprintf('-te %i %i %i %i ',xl(2),yb(2),xr(2),yt(2));
  end

  par=[pr1 pr2 pr3 pr4];
  inv=['"' name '" '];
  if cty==0
    ouv=['"' outpth 'MGSMp' ds '.tif"'];
  else
    ouv=['"' outpth 'CGSMp' ds '.tif"'];
  end
  system([fun par inv ouv]); % project

  delete(name);

  if cty==0
    p1=double(imread([outpth 'MGSMp' ds '.tif']));
  else
    p1=double(imread([outpth 'CGSMp' ds '.tif']));
  end
  p1(p1==ndv)=NaN;

else
  if ~isempty(outpth)
% Creat the asc files
    ds=nm(length(nm)-9:length(nm)); % pay careful attention to this.
    if cty==0
      name=[outpth 'MGSM' ds '.asc'];
    else
      name=[outpth 'CGSM' ds '.asc'];
    end

    fid=fopen(name,'w');
    fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',['ncols ' num2str(nc)],['nrows '...
      num2str(nr)],['xllcorner ' num2str(xll,8)],['yllcorner ' num2str(yll,8)],...
      ['cellsize ' num2str(rs_lat)],['NODATA_value ' num2str(ndv)]);
    fclose(fid);
    dlmwrite(name,p,'delimiter',' ','-append');
  end
end
p(p==ndv)=NaN;
end

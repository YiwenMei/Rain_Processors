% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 11/10/2017

%% Functionality
% This function is used to process the PERSIANN-CCS satellite precipitation
% product (Hong et al. 2007). It includes several functionalities:
%  1)unzip the PERSIANN-CCS record (rgccs1hyydddhh.bin.gz);
%  2)read and crop the record based on an given lat/lon box;
%  3)output the cropped record as .asc file (optional);
%  4)project the record in .asc to another projection and output the record as
%    GTiff (optional).

%% Input
% infname: full name with path of the input PERSIANN-CCS file (e.g.,
%          G:\PERSIANN\CCS\2003\01\rgccs1h0303123.bin.gz);
% workpth: path to store the unzipped record;
%  xl/xr : left/right longitude and coordinate of the boundary (field one must be
%          in the range of [-180 180], field two is in the projected coordinate unit,
%          set field two to NaN if no need to crop the projected image);
%  yb/yt : bottom/top latitude and coordinate of the boundary (field one must be
%          in the range of [60 -60], refer to descriptions for xl/xr for field two);
% outpth : path to store the .asc and, if have, the .tif files (set it to "[]" if
%          no need to output record in .asc format);
% out_pj : output coordinate system (e.g., EPSG:102009; set it to "[]" if no
%          reprojection is required);
%   rs   : x and y resolution of the projected image (set it to "[]" if no
%          reprojection is required).

%% Output
%       p / p1      : cropped precipitation map in original/new projection (the
%                     orientation follows the human reading convention);
%  CCSyyyymmddhh.asc: precipitaiton map in original projection outputted to
%                     out_path as .asc file;
% CCSpyyyymmddhh.tif: precipitaiton map in new projection outputted to out_path
%                     as .tif file;

%% Additional note
% 1)Please make sure to have GDAL installed and the out_path set if you want
%   reproject the data into other coordinate system.
% 2)If reprojection is not required (i.e., out_pj is "[]") but record outputted
%   as .asc is wanted, out_path is required to set.
% 3)The doy2date.m function is required and may be downloaded from
%   https://github.com/sfoerster/matlab/blob/master/supervised_learning/doy2date.m;
% 4)The scale factor and no-data value of PERSIANN-CCS are preseved in the .asc
%   and .tif record.
% 4)No scale factor is preserved in p and p1 and no-data value is replaced by NaN.

function [p,p1]=CCS_process(infname,workpth,xl,xr,yb,yt,outpth,out_pj,rs)
% Lat/lon grids and other info of PERSIANN-CCS
rs_lon=360/9000;
rs_lat=120/3000;
Lon=0:rs_lon:360;
Lat=60:-rs_lat:-60;

ndv=-9999; % no-data value
scf=100; % Scale factor

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
% system(sprintf('7z e "%s" -o"%s" * -r',in_fname,work_path));
gunzip(infname,workpth); % unzip

[~,nm,~]=fileparts(infname);
uz_fn=[workpth nm];
fid=fopen(uz_fn,'r'); % read
p=fread(fid,[9000 3000],'int16','b'); % Format of PERSIANN-CCS
fclose(fid);

p=p';
p=p(rt:rb,cl:cr); % crop

delete(uz_fn);

p1=[];
if ~isempty(out_pj)
% Creat the asc files
  y=str2double(['20' nm(8:9)]);
  jd=str2double(nm(10:12));
  ds=[datestr(doy2date(jd,y),'yyyymmdd') nm(13:14)];
  name=[outpth 'CCS' ds '.asc'];

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

  delete(name);

  par=[pr1 pr2 pr3 pr4];
  inv=['"' name '" '];
  ouv=['"' outpth '\CCSp' ds '.tif"'];
  system([fun par inv ouv]); % project

  p1=double(imread([outpth '\CCSp' ds '.tif']))/scf;
  p1(p1==ndv/scf)=NaN;

else
  if ~isempty(outpth)
% Creat the asc files
    y=str2double(['20' nm(8:9)]);
    jd=str2double(nm(10:12));
    ds=[datestr(doy2date(jd,y),'yyyymmdd') nm(13:14)];
    name=[outpth '\CCS' ds '.asc'];

    fid=fopen(name,'w');
    fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',['ncols ' num2str(nc)],['nrows '...
      num2str(nr)],['xllcorner ' num2str(xll,8)],['yllcorner ' num2str(yll,8)],...
      ['cellsize ' num2str(rs_lat)],['NODATA_value ' num2str(ndv)]);
    fclose(fid);
    dlmwrite(name,p,'delimiter',' ','-append');
  end
end
p(p==ndv)=NaN;
p=p/scf;
end

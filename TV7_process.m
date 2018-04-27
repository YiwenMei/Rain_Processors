% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 3/2/2018

%% Functionality
% This function is used to process the TMPA 3B42-V7 satellite precipitation
% product (Huffuman et al. 2007). It includes several functionalities:
%  1)unzip the TMPA 3B42-V7 record (3B42.2000030100.7R2.bin.gz);
%  2)read and crop the record based on an given lat/lon box;
%  3)output the cropped record as .asc file (optional);
%  4)project the record in .asc to another projection and output the record
%    in GTiff (optional).

%% Input
% infname: full name with path of the input TMPA 3B42-V7 file (e.g.,
%          G:\TMPA\3B42V7\2000\3B42.20000301.00.7A.HDF.Z);
% workpth: path to store the unzipped record;
%  xl/xr : left/right longitude of the boundary (xl/xr can have either 1 or 2 number(s)
%          where the first one represents the longitude and must be in the range of
%          [-180 180]; the second one is boundary in the projected coordinate unit);
%  yb/yt : bottom/top latitude of the boundary (yb/yt can have either 1 or 2 number(s)
%          where the first one represents the latitude and must be in the range of
%          [50 -50]; the second one is boundary in the projected coordinate unit);
% outpth : path to store the .asc and, if have, the .tif files (set it to "[]" if
%          no need to output record in .asc format);
% out_pj : output coordinate system (e.g., EPSG:102009; set it to "[]" if no
%          reprojection is required);
%   rs   : x and y resolution of the projected image (set it to "[]" if no
%          reprojection is required).

%% Output
%         p / p1       : cropped precipitation map in original/new projection (the
%                        orientation follows the human reading convention);
%  3B42V7yyyymmddhh.asc: precipitaiton map in original projection and resolution
%                        outputted to outpth as .asc file;
% 3B42V7pyyyymmddhh.tif: precipitaiton map in new projection, resolution and extend
%                        outputted to outpth as .tif file;

%% Additional note
% 1)Please make sure to have GDAL installed and the out_path set if you want
%   reproject the data into other coordinate system.
% 2)If reprojection is not required (i.e., out_pj and rs are "[]") but record
%   outputted as .asc is wanted, out_path is required to set.
% 3)The no-data value of TMPA 3B42V7 are changed (from -9.9999...e3) to -9999 in
%   the .asc and .tif record.

function [p,p1]=TV7_process(infname,workpth,xl,xr,yb,yt,outpth,out_pj,rs)
% Lat/lon grids and other info of TMPA 3B42V7
res_lon=360/1440;
res_lat=100/400;
Lon=-180:res_lon:180;
Lat=50:-res_lat:-50;

ndv=-9999; % no-data value

% Index of interested domain
cl=find(xl(1)-Lon>=0,1,'last'); % left column
cr=find(xr(1)-Lon<=0,1,'first')-1; % right column
rt=find(yt(1)-Lat<=0,1,'last'); % top row
rb=find(yb(1)-Lat>=0,1,'first')-1; % bottom row

nr=length(rt:rb); % number of row
nc=length(cl:cr); % number of column
xll=-180+(cl-1)*res_lon; % longitude of lower left corner
yll=50-rb*res_lat; % latitude of lower left corner

% unzip and read the input
system(sprintf('7z e "%s" -o"%s" * -r',infname,workpth));
% gunzip(infname,workpth); % unzip

[~,nm,~]=fileparts(infname);
uz_fn=[workpth nm];
p=hdfread(uz_fn,'precipitation');
p=rot90(p); % Original upper-left corner is (S,W). Rotate 90 degree counter
p(p<0)=ndv; %  clock-wise to (W,N).
p=p(rt:rb,cl:cr); % crop
delete(uz_fn);

p1=[];
if ~isempty(out_pj)
% Creat the asc files
  ds=nm([6:13 15:16]); % pay careful attention to this.
  name=[outpth '3B42V7' ds '.asc'];

  fid=fopen(name,'w');
  fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',['ncols ' num2str(nc)],['nrows '...
    num2str(nr)],['xllcorner ' num2str(xll,8)],['yllcorner ' num2str(yll,8)],...
    ['cellsize ' num2str(res_lat)],['NODATA_value ' num2str(ndv)]);
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
  ouv=['"' outpth '3B42V7p' ds '.tif"'];
  system([fun par inv ouv]); % project

  delete(name);

  p1=double(imread([outpth '3B42V7p' ds '.tif']));
  p1(p1==ndv)=NaN;

else
  if ~isempty(outpth)
% Creat the asc files
    ds=nm(8:end-8);
    name=[outpth '3B42V7' ds '.asc'];

    fid=fopen(name,'w');
    fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',['ncols ' num2str(nc)],['nrows '...
      num2str(nr)],['xllcorner ' num2str(xll,8)],['yllcorner ' num2str(yll,8)],...
      ['cellsize ' num2str(res_lat)],['NODATA_value ' num2str(ndv)]);
    fclose(fid);
    dlmwrite(name,p,'delimiter',' ','-append');
  end
end
p(p==ndv)=NaN;
end

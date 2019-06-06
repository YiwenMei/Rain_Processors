% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 9/17/2018

%% Functionality
% This function is used to process the APHRODITE precipitation product (Yatagai
% et al. 2011). It includes several functionalities:
%  1)unzip the APHRODITE record (APHRO_MA_025deg_V1101.2001.nc.gz);
%  2)read and crop the record based on a given lat/lon box;
%  3)output the cropped record as .asc file (optional);
%  4)project the record in .asc to another projection and output the record as
%    GTiff (optional).

%% Input
% infname: full name with path of the input APHRODITE record (e.g.,
%          G:\APHRODITE\APHRO_MA_025deg_V1101.2001.nc.gz;
% workpth: path to store the unzipped record;
%   fty  : field of the record needs to be extracted (0 for the ratio of station,
%          1 for precipitation rate);
%  xl/xr : left/right longitude of the boundary (xl/xr can have either 1 or 2
%          number(s) where the first one represents the longitude in degree;
%          the second one is boundary in the projected coordinate unit;
%  yb/yt : bottom/top latitude of the boundary (yb/yt can have either 1 or 2
%          number(s) where the first one represents the latitude in degree; the
%          second one is boundary in the projected coordinate unit);
% outpth : path to store the .asc and, if have, the .tif files (set it to "[]"
%          if no need to output record in .asc format);
% out_pj : output coordinate system (e.g., EPSG:102009; set it to "[]" if no
%          reprojection is required);
%   rs   : x and y resolution of the projected image (set it to "[]" if no
%          reprojection is required).

%% Output
% APH_pryyyymmdd.asc : precipitaiton map in original projection and resolution
%                      outputted to outpth as .asc file;
% APHp_pryyyymmdd.tif: precipitaiton map in new projection, resolution and extend
%                      outputted to outpth as .tif file;
% APH_rsyyyymmdd.asc : ratio of station map in original projection and resolution
%                      outputted to outpth as .asc file;
% APHp_rsyyyymmdd.tif: ratio of station map in new projection, resolution and
%                      extend outputted to outpth as .tif file;

%% Additional note
% 1)Please make sure to have GDAL installed and the out_path set if you want
%   to reproject the data into other coordinate system.
% 2)If reprojection is not required but record outputted as .asc is wanted, outpth
%   is required to set.
% 3)The no-data value of APHRODITE are preseved in the .asc and .tif record;
% 4)Spatial interpolation for ratio of station is nearest.

function APH_process(infname,workpth,fty,xl,xr,yb,yt,outpth,out_pj,rs)
% Unzip the record
gunzip(infname,workpth);
[~,nm,~]=fileparts(infname);
uz_fn=[workpth nm];
ys=nm(23:26); % pay careful attention to this.

% Lat/lon grids and other info of APHRODITE
Lat=double(ncread(uz_fn,'lat'));
Lon=double(ncread(uz_fn,'lon'));
rso=.25;
Lat=flipud(unique([Lat-rso/2,Lat+rso/2]));
Lon=unique([Lon-rso/2,Lon+rso/2]);

if fty==1
  fn='precip'; % Field of interest
elseif fty==0
  fn='rstn';
end
T=double(ncread(uz_fn,'time'));
ndv=ncinfo(uz_fn,fn);
ndv=double(ndv.Attributes(1).Value); % no-data value

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
xll=min(Lon)+(cl-1)*rso; % longitude of lower left corner
yll=max(Lat)-rb*rso; % latitude of lower left corner

% Read the record
parfor d=1:length(T)
  p=double(rot90(ncread(uz_fn,fn,[1,1,d],[length(Lon)-1 length(Lat)-1 1])));
  p=p(rt:rb,cl:cr); % crop

  if ~isempty(out_pj)
% Creat the asc files
    ds=datestr(doy2date(d,str2double(ys)),'yyyymmdd');
    name=[outpth 'APH_' fn(1:2) ds '.asc'];

    fid=fopen(name,'w');
    fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',['ncols ' num2str(nc)],['nrows '...
      num2str(nr)],['xllcorner ' num2str(xll,8)],['yllcorner ' num2str(yll,8)],...
      ['cellsize ' num2str(rso)],['NODATA_value ' num2str(ndv)]);
    fclose(fid);
    dlmwrite(name,p,'delimiter',' ','-append'); % output .asc

% Project to a new coordinate
    if fty==1
      fun='gdalwarp -overwrite -of GTiff -r bilinear '; % GDAL function
    elseif fty==0
      fun='gdalwarp -overwrite -of GTiff '; % GDAL function
    end
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
    ouv=['"' outpth 'APHp_' fn(1:2) ds '.tif"'];
    system([fun par inv ouv]); % project

    delete(name);

  else
    if ~isempty(outpth)
% Creat the asc files
      ds=datestr(doy2date(d,str2double(ys)),'yyyymmdd');
      name=[outpth 'APH_' fn(1:2) ds '.asc'];

      fid=fopen(name,'w');
      fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',['ncols ' num2str(nc)],['nrows '...
        num2str(nr)],['xllcorner ' num2str(xll,8)],['yllcorner ' num2str(yll,8)],...
        ['cellsize ' num2str(rso)],['NODATA_value ' num2str(ndv)]);
      fclose(fid);
      dlmwrite(name,p,'delimiter',' ','-append');
    end
  end
end
delete(uz_fn);
end

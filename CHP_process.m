% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 8/11/2018

%% Functionality
% This function is used to process the satellite-only or satellite-gauge CHIRPS
% precipitation product (Funk et al. 2015). It includes several functionalities:
%  1)read and crop the record based on a given lat/lon box;
%  2)output the cropped record as .asc file (optional);
%  3)project the record in .asc to another projection and output the record as
%    GTiff (optional).

%% Input
% infname: full name with path of the input CHIRP or CHIRPS record (e.g.,
%          G:\CHP\CHIRP\chirp.2007.days_p05.nc;
%          G:\CHP\CHIRPS\chirps-v2.0.2007.days_p05.nc);
%   cty  : type of CHP record (0 for CHIRP, 1 for CHIRPS);
%  xl/xr : left/right longitude of the boundary (xl/xr can have either 1 or 2
%          number(s) where the first one represents the longitude and must be
%          in the range of [-180 180]; the second one is boundary in the projected
%          coordinate unit);
%  yb/yt : bottom/top latitude of the boundary (yb/yt can have either 1 or 2
%          number(s) where the first one represents the latitude and must be
%          in the range of [50 -50]; the second one is boundary in the projected
%          coordinate unit);
% outpth : path to store the .asc and, if have, the .tif files (set it to "[]"
%          if no need to output record in .asc format);
% out_pj : output coordinate system (e.g., EPSG:102009; set it to "[]" if no
%          reprojection is required);
%   rs   : x and y resolution of the projected image (set it to "[]" if no
%          reprojection is required).

%% Output
%  CHIRPyyyymmddhh.asc : precipitaiton map in original projection and resolution
%                        outputted to outpth as .asc file;
% CHIRPSpyyyymmddhh.tif: precipitaiton map in new projection, resolution and
%                        extend outputted to outpth as .tif file;

%% Additional note
% 1)Please make sure to have GDAL installed and the out_path set if you want
%   to reproject the data into other coordinate system.
% 2)If reprojection is not required but record outputted as .asc is wanted, outpth
%   is required to set.
% 3)The no-data value of CHIRP/CHIRPS are preseved in the .asc and .tif record.

function CHP_process(infname,cty,xl,xr,yb,yt,outpth,out_pj,rs)
% Lat/lon grids and other info of CHIRP/CHIRPS
Lat=double(ncread(infname,'latitude'));
Lon=double(ncread(infname,'longitude'));
rs_lon=360/length(Lon);
rs_lat=100/length(Lat);
Lon=0:rs_lon:360;
Lat=50:-rs_lat:-50;

[~,nm,~]=fileparts(infname);
if cty==0
  ys=nm(7:10); % pay careful attention to this.
else
  ys=nm(13:16);
end

T=double(ncread(infname,'time'));
ndv=ncinfo(infname,'precip');
ndv=ndv.FillValue; % no-data value

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
yll=50-rb*rs_lat; % latitude of lower left corner

% Read the record
for d=1:length(T)
  p=rot90(ncread(infname,'precip',[1,1,d],[length(Lon)-1 length(Lat)-1 1]));
  p=[p(:,3601:7200) p(:,1:3600)]; % convert from [-180 180] to [0 360]
  p=p(rt:rb,cl:cr); % crop
  p(isnan(p))=ndv;

  if ~isempty(out_pj)
% Creat the asc files
    ds=datestr(doy2date(d,str2double(ys)),'yyyymmdd');
    if cty==0
      name=[outpth 'CHIRP' ds '.asc'];
    else
      name=[outpth 'CHIRPS' ds '.asc'];
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
      ouv=['"' outpth 'CHIRPp' ds '.tif"'];
    else
      ouv=['"' outpth 'CHIRPSp' ds '.tif"'];
    end
    system([fun par inv ouv]); % project

    delete(name);

  else
    if ~isempty(outpth)
% Creat the asc files
      ds=datestr(doy2date(d,str2double(ys)),'yyyymmdd');
      if cty==0
        name=[outpth 'CHIRP' ds '.asc'];
      else
        name=[outpth 'CHIRPS' ds '.asc'];
      end

      fid=fopen(name,'w');
      fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',['ncols ' num2str(nc)],['nrows '...
        num2str(nr)],['xllcorner ' num2str(xll,8)],['yllcorner ' num2str(yll,8)],...
        ['cellsize ' num2str(rs_lat)],['NODATA_value ' num2str(ndv)]);
      fclose(fid);
      dlmwrite(name,p,'delimiter',' ','-append');
    end
  end
end
end

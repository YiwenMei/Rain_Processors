function [rc,rc1]=CHP_Gnum_process(cfn,wkpth,xl,xr,yb,yt,opth,out_pj,rs)
% Lat/lon grids and other info of CMORPH
xl0=0;
xr0=360;
yt0=90;
yb0=-90;
rs_lat=0.5;
rs_lon=0.5;
Lon=xl0:rs_lon:xr0;
Lat=yt0:-rs_lat:yb0;
ndv=-999;
ors='wgs84';

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
yll=yt0-rb*rs_lat; % latitude of lower left corner

% Read the csv record
fid=fopen(cfn);
rc=textscan(fid,'%s%s%s%s%s%s%s','Delimiter',{',',',',',',',',',',',',','});
fclose(fid);

% Construct the spatial density map of gauge
lat=cellfun(@str2double,rc{3}(2:end));
lon=cellfun(@str2double,rc{4}(2:end));
lon(lon<0)=lon(lon<0)+360; % Convert longitude to 0 360 deg
N_lat=lat<(yt0:-rs_lat:yb0+rs_lat) & lat>=(yt0-rs_lat:-rs_lat:yb0);
N_lon=lon>=(xl0:rs_lon:xr0-rs_lon) & lon<(xl0+rs_lon:rs_lon:xr0);
rc=double(N_lat')*double(N_lon);

% Crop gauge density record
rc=rc(rt:rb,cl:cr);

[~,nm,~]=fileparts(cfn);
name=[wkpth 'CHPgn' nm([1:4 6 7]) '.asc'];
ouv=[opth 'CHPgnp' nm([1:4 6 7]) '.tif'];

if ~isempty(out_pj)
% Creat the asc files
  fid=fopen(name,'w');
  fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',['ncols ' num2str(nc)],['nrows '...
      num2str(nr)],['xllcorner ' num2str(xll,8)],['yllcorner ' num2str(yll,8)],...
      ['cellsize ' num2str(rs_lat)],['NODATA_value ' num2str(ndv)]);
  fclose(fid);
  dlmwrite(name,rc,'delimiter',' ','-append'); % output .asc

% Project to a new coordinate
  fun='gdalwarp -overwrite -of GTiff -r near '; % GDAL function
  pr1=['-s_srs ' ors ' ']; % Projection of original record
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
  system([fun par '"' name '" "' ouv '"']); % project
  delete(name);

  rc1=double(imread(ouv));
  rc1(rc1==ndv)=NaN;

else
  if ~isempty(opth)
% Creat the asc files
    matV2tif(ouv,rc,xll,yll,rs_lat,ndv,ors,wkpth);
  end
end
end

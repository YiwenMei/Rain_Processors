function [p,p1]=CPCU_process(infname,cty,wkpth,xl,xr,yb,yt,oupth,out_pj,rs)
% Lat/lon grids and other info of CMORPH
rs_lon=360/720;
rs_lat=180/360;
Lon=0:rs_lon:360;
Lat=90:-rs_lat:-90;

ndv=-999; % no-data value

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
yll=90-rb*rs_lat; % latitude of lower left corner

% unzip and read the input
bfn=cell2mat(gunzip(infname));
movefile(bfn,wkpth);
nm=cell2mat(regexp(bfn,'(?<year>\d{4})(?<month>\d{2})(?<date>\d{2})','match'));
bfn=dir([wkpth '*' nm '*']);
bfn=fullfile(wkpth,bfn.name);

fid=fopen(bfn);
p=fread(fid,'float32','l');
fclose(fid);
p=reshape(p,720,360,2);

switch cty
  case 'pr'
    p=rot90(p(:,:,1));
    p=p(rt:rb,cl:cr); % crop
    itm='bilinear';
    scf=.1;
  case 'gn'
    p=rot90(p(:,:,2));
    p=p(rt:rb,cl:cr); % crop
    itm='near';
    scf=1;
end
delete(bfn);

p1=[];
if ~isempty(out_pj)
% Creat the asc files
  name=[wkpth 'CPCU' cty nm '.asc'];

  fid=fopen(name,'w');
  fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',['ncols ' num2str(nc)],['nrows '...
    num2str(nr)],['xllcorner ' num2str(xll,8)],['yllcorner ' num2str(yll,8)],...
    ['cellsize ' num2str(rs_lat)],['NODATA_value ' num2str(ndv)]);
  fclose(fid);
  dlmwrite(name,p,'delimiter',' ','-append'); % output .asc

% Project to a new coordinate
  fun='gdalwarp -overwrite -of GTiff -s_srs wgs84 '; % GDAL function
  pr1=['-r ' itm ' ']; % Projection of original record
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
  ouv=[oupth 'CPCU' cty 'p' nm '.tif'];
  system([fun par '"' name '" "' ouv '"']); % project
  delete(name);

  p1=double(imread(ouv));
  p1(p1==ndv)=NaN;
  p1=p1*scf;

else
  if ~isempty(oupth)
% Creat the asc files
    name=[wkpth 'CPCU' cty nm '.asc'];

    fid=fopen(name,'w');
    fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',['ncols ' num2str(nc)],['nrows '...
      num2str(nr)],['xllcorner ' num2str(xll,8)],['yllcorner ' num2str(yll,8)],...
      ['cellsize ' num2str(rs_lat)],['NODATA_value ' num2str(ndv)]);
    fclose(fid);
    dlmwrite(name,p,'delimiter',' ','-append');
  end
end
p(p==ndv)=NaN;
p=p*scf;
end

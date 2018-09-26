% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 09/25/2018

%% Functionality
% This function includes several functionalities:
%  1)read all the tiles of MOD11A1 Band 31 & 32 emissivity for the study area;
%  2)calculate the emissivity based on the mean between the two bands;
%  3)mosaic all the emissivity tiles;
%  4)reproject and crop the mosaiced emissivity.

%% Input
%   vfl  : strings of the list of MOD11A1 files;
%  oupth : path to store the output emissivity geotiff;
% rxo/ryo: resolution of interest in output coordinate unit;
%  xl/xr : left/right corner coordinate in output coordinate unit;
%  yb/yt : same as xl/xr but for the bottom and top corner coordinate;
%  coor  : coordinate system of interested;
%   ndv  : no-data value assigned to the output image.

%% Output
% Output is a mosaiced emissivity map stored in oupth. The naming convention
%  is EMSyyyyddd.tif the same as the input date.

function Emis_process(efl,oupth,rxo,ryo,xl,xr,yb,yt,coor,ndv)
% Domain info
hi=nan(size(efl,1),1);
vi=nan(size(efl,1),1);
for n=1:size(efl,1)
  vfn=efl(n,:);
  [~,ns,~]=fileparts(vfn);
  hi(n)=str2double(ns(19:20));
  vi(n)=str2double(ns(22:23));
end
ntr=length(unique(vi)); % num. of tile in row
ntc=length(unique(hi)); % num. of tile in column
hi=min(hi)-1; % initial column tile num.
vi=min(vi)-1; % initial row tile num.

% VI tile info
hif=hdfinfo(vfn);
ivn1=hif.Vgroup.Vgroup(1).SDS(9).Name; % Band 31 emissivity
ivn2=hif.Vgroup.Vgroup(1).SDS(10).Name; % Band 32 emissivity
hif=hif.Vgroup.Vgroup(1).SDS(9);
nr=hif.Dims(1).Size; % num. of grid in row of a tile
nc=hif.Dims(2).Size; % num. of grid in column of a tile
scf=double(hif.Attributes(6).Value); % scale factor
ofs=double(hif.Attributes(7).Value); % additive offset
ndv_o=double(hif.Attributes(4).Value); % no-data value of VI
hif=hdfinfo(vfn,'eos');
hif=hif.Grid;
rtx=floor(rxo*hif.Columns/(hif.LowerRight(1)-hif.UpperLeft(1))); % Upscale ratio of x
rty=floor(ryo*hif.Rows/(hif.UpperLeft(2)-hif.LowerRight(2))); % Upscale ratio of y

%% Mosaic the EMS tiles
Em=ndv*ones(ntr*nr/rty,ntc*nc/rtx);
xyll=nan(size(efl,1),4);

for n=1:size(efl,1)
% Read the emissivity
  vfn=efl(n,:);
  em1=double(hdfread(vfn,ivn1));
  em1(em1==ndv_o)=NaN;
  em2=double(hdfread(vfn,ivn2));
  em2(em2==ndv_o)=NaN;

% max emissivity
  em=nanmean(cat(3,em1,em2),3)*scf+ofs;
  em(isnan(em))=ndv;
  clear em1 em2

% Upscale emissivity
  hif=hdfinfo(vfn,'eos');
  hif=hif.Grid;
  xlt=hif.UpperLeft(1);
  xrt=hif.LowerRight(1);
  rx=(xrt-xlt)/hif.Columns; % resolution of a EMS tile
  ybt=hif.LowerRight(2);
  ytt=hif.UpperLeft(2);

% Mosaic the EMS tiles
  [~,ns,~]=fileparts(vfn);
  h=str2double(ns(19:20))-hi;
  v=str2double(ns(22:23))-vi;
  Em((v-1)*size(em,1)+1:v*size(em,2),(h-1)*size(em,1)+1:h*size(em,1))=em;

  xyll(n,:)=[xlt ytt xrt ybt];
end

%% Project the EMS
% Write asc
ncn=[oupth 'p.asc'];

fid=fopen(ncn,'w');
fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',['ncols ' num2str(size(Em,2))],...
    ['nrows ' num2str(size(Em,1))],['xllcorner ' num2str(min(xyll(:,1)),12)],...
    ['yllcorner ' num2str(min(xyll(:,4)),12)],['cellsize ' num2str(rtx*rx,12)],...
    ['NODATA_value ' num2str(ndv)]);
dlmwrite(ncn,Em,'delimiter',' ','-append');
fclose(fid);

clear Em ems

% Project the EMS
fun='gdalwarp -overwrite -of GTiff -r bilinear '; % GDAL function
pr1=['-t_srs ' coor ' '];
pr2=sprintf('%s %i %i %i %i ','-te',xl,yb,xr,yt);
pr3=sprintf('%s %i %i ','-tr',rxo,ryo);
pr4='-s_srs "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" ';
pr5=sprintf('%s %i ','-srcnodata',ndv);

par=[pr1 pr2 pr3 pr4 pr5];
ouv=[oupth 'EMS' ns(10:16) '.tif'];
system([fun par '"' ncn '" "' ouv '"']);

delete(ncn);
end

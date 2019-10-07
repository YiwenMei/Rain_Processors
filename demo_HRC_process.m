%% Untar the file
inpth='/scratch/clewis22/Cropping/Data';
zfl=dir(fullfile(inpth,'*.tar'));
zfn=fullfile(inpth,zfl.name);
untar(zfn,inpth);
fn=cell2mat(regexp(zfn,'_(\d{6}).tar','tokens','once'));
bzl=dir(fullfile(inpth,fn,'*.bz2'));

%% Example 1: Read and extract the file
bzn=fullfile(inpth,fn,bzl(3).name);
wkpth=fullfile(inpth,'wkdir');
xl=61.27;
xr=90.15;
yb=19.87;
yt=41.16;

p1=HRC_process(bzn,'HRC',wkpth,xl,xr,yb,yt);

%% Example 2: Read, extract, and output the file as geotiff
bzn=fullfile(inpth,fn,bzl(6).name);
opth=fullfile(inpth,'outpath');

p2=HRC_process(bzn,'HRC',wkpth,xl,xr,yb,yt,opth);

%% Example 3: Read, extract, project, and output the file as geotiff
bzn=fullfile(inpth,fn,bzl(7).name);
ors='EPSG:102012';

p3=HRC_process(bzn,'HRC',wkpth,xl,xr,yb,yt,opth,ors);

%% Example 4: Read, extract, resample, and output the file as geotiff
bzn=fullfile(inpth,fn,bzl(8).name);
ors='wgs84';
rs=[.1 .2];

p4=HRC_process(bzn,'HRC',wkpth,xl,xr,yb,yt,opth,ors,rs);

%% Example 5: Read, extract, resample, project, and output the file as geotiff
bzn=fullfile(inpth,fn,bzl(10).name);
ors='EPSG:102012';
rs=[10000 10000];

p5=HRC_process(bzn,'HRC',wkpth,xl,xr,yb,yt,opth,ors,rs);

%% Example 6: Read, extract, resample, project, crop, and output the file as geotiff
bzn=fullfile(inpth,fn,bzl(11).name);
xl=[61.27 -3850000];
xr=[90.15 -1500000];
yb=[19.87 3450000];
yt=[41.16 5100000];
ors='EPSG:102012';
rs=[10000 10000];

p6=HRC_process(bzn,'HRC',wkpth,xl,xr,yb,yt,opth,ors,rs);

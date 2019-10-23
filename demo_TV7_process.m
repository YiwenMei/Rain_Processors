inpth='/scratch/ymei2/Example';

%% Example 1 - Extracting
fn=fullfile(inpth,'3B42.20100101.00.7A.HDF.Z');
wkpth=fullfile(inpth,'wkdir');
xl=61.27;
xr=90.15;
yb=19.87;
yt=41.16;

p1=TV7_process(fn,wkpth,xl,xr,yb,yt);

%% Example 2 - Extracting and Outputting as geotif
opth=fullfile(inpth,'outfile');

[p2,ofn]=TV7_process(fn,wkpth,xl,xr,yb,yt,opth);

%% Example 3 - Extracting, Resampling, and Outputting as geotif
rs=[.5 .5];

[p3,~]=TV7_process(fn,wkpth,xl,xr,yb,yt,opth,'wgs84',rs);

%% Example 4 - Extracting, Projecting, and Outputting as geotif
ors='EPSG:102012'; % Asia Lambert Conformal Conic

[p4,~]=TV7_process(fn,wkpth,xl,xr,yb,yt,opth,ors);

%% Example 5 - Extracting, Projecting, Resampling, and Outputting as geotif
rs=[40000 40000];

[p5,~]=TV7_process(fn,wkpth,xl,xr,yb,yt,opth,ors,rs);

%% Example 6 - Extracting, Projecting, Cropping, and Outputting as geotif
xl=[61.27 -3850000];
xr=[90.15 -1500000];
yb=[19.87 3450000];
yt=[41.16 5100000];

p6=TV7_process(fn,wkpth,xl,xr,yb,yt,opth,ors);

%% Example 7 - Extracting, Projecting, Cropping, Resampling, and Outputting as geotif
p7=TV7_process(fn,wkpth,xl,xr,yb,yt,opth,ors,rs);

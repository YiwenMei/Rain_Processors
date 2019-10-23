inpth='/scratch/ymei2/Example';
fn=fullfile(inpth,'chirp.2010.days_p05.nc');
% fn=fullfile(inpth,'chirps-v2.0.2010.days_p05.nc');

%% Example 1 - Extracting and Outputting as .mat
wkpth=fullfile(inpth,'wkdir');
opth=fullfile(inpth,'outfile');
xl=61.27;
xr=90.15;
yb=19.87;
yt=41.16;

Ofn=CHP_process(fn,wkpth,opth,'CHIRP',xl,xr,yb,yt);

%% Example 2 - Extracting, Resampling, and Outputting as .tif
rs=[.1 .2];

Ofn=CHP_process(fn,wkpth,opth,'CHIRP',xl,xr,yb,yt,'wgs84',rs);

%% Example 3 - Extracting, Projecting, and Outputting as geotif
ors='EPSG:102012';

Ofn=CHP_process(fn,wkpth,opth,'CHIRP',xl,xr,yb,yt,ors);

%% Example 4 - Extracting, Projecting, Resampling, and Outputting as geotif
rs=[12000 10000];

Ofn=CHP_process(fn,wkpth,opth,'CHIRP',xl,xr,yb,yt,ors,rs);

%% Example 5 - Extracting, Projecting, Cropping, and Outputting as geotif
xl=[61.27 -3850000];
xr=[90.15 -1500000];
yb=[19.87 3450000];
yt=[41.16 5100000];

Ofn=CHP_process(fn,wkpth,opth,'CHIRP',xl,xr,yb,yt,ors);

%% Example 6 - Extracting, Projecting, Cropping, Resampling, and Outputting as geotif
Ofn=CHP_process(fn,wkpth,opth,'CHIRP',xl,xr,yb,yt,ors,rs);

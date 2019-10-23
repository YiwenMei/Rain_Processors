%% Untar the file
inpth='/scratch/ymei2/Example';
zfl=dir(fullfile(inpth,'*.tar'));
zfn=fullfile(inpth,zfl.name);
fn=untar(zfn,inpth);
fn=fn(2:end);

%% Example 1 - Extracting
wkpth=fullfile(inpth,'wkdir');
xl=61.27;
xr=90.15;
yb=19.87;
yt=41.16;

p1=HRC_process(fn{1},wkpth,xl,xr,yb,yt);

%% Example 2 - Extracting and Outputting as geotif
opth=fullfile(inpth,'outfile');

p2=HRC_process(fn{1},wkpth,xl,xr,yb,yt,opth,'HRC');

%% Example 3 - Extracting, Resampling, and Outputting as geotif
rs=[.1 .2];

p3=HRC_process(fn{1},wkpth,xl,xr,yb,yt,opth,'HRC','wgs84',rs);

%% Example 4 - Extracting, Projecting, and Outputting as geotif
ors='EPSG:102012';

p4=HRC_process(fn{1},wkpth,xl,xr,yb,yt,opth,'HRC',ors);

%% Example 5 - Extracting, Projecting, Resampling, and Outputting as geotif
rs=[12000 10000];

p5=HRC_process(fn{1},wkpth,xl,xr,yb,yt,opth,'HRC',ors,rs);

%% Example 6 - Extracting, Projecting, Cropping, and Outputting as geotif
xl=[61.27 -3850000];
xr=[90.15 -1500000];
yb=[19.87 3450000];
yt=[41.16 5100000];

p6=HRC_process(fn{1},wkpth,xl,xr,yb,yt,opth,'HRC',ors);

%% Example 7 - Extracting, Projecting, Cropping, Resampling, and Outputting as geotif
p7=HRC_process(fn{1},wkpth,xl,xr,yb,yt,opth,'HRC',ors,rs);

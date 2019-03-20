%% Set working directory and program path 
cd 'C:\Brandon\Manuscripts\Overstreet2018_SurfReflect\Analysis'
addpath(genpath('C:\Brandon\mfiles\'))



% bring in the roughness data
filename = 'SurfRfl_CASI2014_RoughnessV4.mat';
rcR      = matfile(filename,'Writable',true);
% check what's inside
whos(rcR)

% bring in the morpho data 
filename = 'SurfRfl_CASI2014_Morphov2.mat';
rcM      = matfile(filename,'Writable',true);
% check what's inside
whos(rcM)
CLProf = rcM.CLProf;
% Load the model outputs
MOut = rcM.MOut;

filename = 'SR2014_RC_ImgData.mat';
rcI      = matfile(filename,'Writable',true);
% check what's inside
whos(rcI)

% The morpho structure contains a number of images which are a culmination
% of the modeling analysis. We won't load those in yet but they are:
% imgMU = morphologic unit map
% imgD  = model derived water depth
% imgNIR = band 34 (846 nm)
% imgSTD = variance of NIR band calculated using stdfilt function
% imgU  = modeled velocity
% imgS  = modeled shear stress
% ImgD = table containing decimated img-derived depths

% Pull in the commonly-used image variables
imXdata = rcM.imXdata;
imYdata = rcM.imYdata;
wvl     = rcM.wvl;

% bring in the standard colormap (1,pool,blue;2,run,yellow;3 riffle,purple;
% 4 shallow,green; 5 backwater, light blue
map = rcM.map;


%% Clip the flow model to the active channel
LidarPoly = shaperead('Lidar\RCLidarChNoZ.shp');
figure
plot(LidarPoly(1).X,LidarPoly(1).Y);

%Remove points outside the channel boundary
iFP = ~inpoly([MOut.X,MOut.Y],[LidarPoly(1).X(~isnan(LidarPoly(1).X))',LidarPoly(1).Y(~isnan(LidarPoly(1).Y))']);

%Pull out the islands
for i = 2:length(LidarPoly)
    iIs = inpoly([MOut.X,MOut.Y],[LidarPoly(i).X(~isnan(LidarPoly(i).X))',LidarPoly(i).Y(~isnan(LidarPoly(i).Y))']);
    iFP(iIs) = 1;
end

MOut(iFP,:) = [];

%% Remove points with zero depth
iZ = MOut.Depth == 0;
MOut(iZ,:) = [];





%% MU classification
% Assign Velocity Codes (1-5)to each grid cell
% 0 - 0.10 = 1;
% 0.11 - 0.25 = 2;
% 0.26 - 0.40 = 3;
% 0.41 - 0.75 = 4;
% >0.75 = 5;
% v1Max = 0.1;
% v2Max = 0.25;
% v3Max = 0.4;
% v4Max = 0.75;

v1Max = 1.1;
v2Max = 1.3;
v3Max = 1.4;
v4Max = 1.5;

%find indices cooresponding to velocity ranges
i1 = MOut.Velocitymagnitude > 0 & MOut.Velocitymagnitude <= v1Max;
i2 = MOut.Velocitymagnitude > v1Max & MOut.Velocitymagnitude <= v2Max;
i3 = MOut.Velocitymagnitude > v2Max & MOut.Velocitymagnitude <= v3Max;
i4 = MOut.Velocitymagnitude > v3Max & MOut.Velocitymagnitude <= v4Max;
i5 = MOut.Velocitymagnitude > v4Max;

% Create a new variable ncV and assign code
ncV = zeros(size(MOut.Velocitymagnitude));

ncV(i1) = 1;
ncV(i2) = 2;
ncV(i3) = 3;
ncV(i4) = 4;
ncV(i5) = 5;

clearvars v1Max v2Max v3Max v4Max
%% Assign Depth Codes (1-5) to each grid cell
% 0.01 - 0.4 = 5;
% 0.41 - 0.8 = 4;
% 0.81 - 1.20 = 3;
% 1.21 - 1.50 = 4;
% >1.50       = 1;

% d1Max = 0.4;
% d2Max = 0.8;
% d3Max = 1.2;
% d4Max = 1.5;

d1Max = 1.1;
d2Max = 1.3;
d3Max = 1.6;
d4Max = 2;

%find indices cooresponding to depth ranges
i5 = MOut.Depth > 0 & MOut.Depth <= d1Max;
i4 = MOut.Depth > d1Max & MOut.Depth <= d2Max;
i3 = MOut.Depth > d2Max & MOut.Depth <= d3Max;
i2 = MOut.Depth > d3Max & MOut.Depth <= d4Max;
i1 = MOut.Depth > d4Max;


% Create a new variable ncV and assign code
ncD = zeros(size(MOut.Velocitymagnitude));

ncD(i1) = 1;
ncD(i2) = 2;
ncD(i3) = 3;
ncD(i4) = 4;
ncD(i5) = 5;

clearvars d1Max d2Max d3Max d4Max
%% Assign Shear Codes (0 - 2) to each grid cell
% 0.00 - 2.0 = 0
% 2.01 - 20 = 1
% >20 = 2
t1Max = 1.5;
t2Max = 20;

i0 = MOut.ShearStressmagnitude > 0 & MOut.ShearStressmagnitude <= t1Max;
i1 = MOut.ShearStressmagnitude > t1Max & MOut.ShearStressmagnitude <= t2Max;
i2 = MOut.ShearStressmagnitude > t2Max;

% % Create a new variable NC.t to hold the shear stress codes
% NC.t = zeros(height(MOut),1);
% %NC.t(i0) = 0;
% NC.t(i1) = 1;
% NC.t(i2) = 2;

% Create a new variable ncV and assign code
ncS = zeros(size(MOut.Velocitymagnitude));

ncS(i0) = 0;
ncS(i1) = 1;
ncS(i2) = 2;

%% Calculate mH
%NC.mH = (NC.d + NC.v).*NC.t;
mH = (ncD + ncV).*ncS;

% Assign unit values based on the following:
% pool (MH 2-4), run (MH: 5-9), fast run (MH: 10-18), riffle (MH: 20), low energy
% (MH = 0)
% 1 = pool
% 2 = run
% 3 = riffle
% 4 = low energy (shallow or backwater)
%mh1Max = 
i1 = mH >= 2 & mH <= 6;
i2 = mH >= 7 & mH <= 10; 
i3 = mH >= 11 & mH <= 20;
%i4 = mH == 20;
% Specify the shallow/backwater threshold
dT = 1;
i4 = mH == 0 & MOut.Depth > 0 & MOut.Depth <= dT;
i5 = mH == 0 & MOut.Depth > 0 & MOut.Depth > dT;

mU = zeros(size(MOut.Depth));

mU(i1) = 1;
mU(i2) = 2;
mU(i3) = 3;
mU(i4) = 4;
mU(i5) = 5;

%% Plot the output


figure
scatter(MOut.X(i1),MOut.Y(i1),10,map(1,:),'filled')
hold
scatter(MOut.X(i2),MOut.Y(i2),10,map(2,:),'filled')
scatter(MOut.X(i3),MOut.Y(i3),10,map(3,:),'filled')
scatter(MOut.X(i4),MOut.Y(i4),10,map(4,:),'filled')
scatter(MOut.X(i5),MOut.Y(i5),10,map(5,:),'filled')


lgd = legend('Pool','Run','Riffle','Shallow','Backwater');
axis equal


%% Create a new table to hold MU values and other data
MU = table;
MU.x = MOut.X;
MU.y = MOut.Y;
MU.modU = MOut.Velocitymagnitude;
MU.modD = MOut.Depth;
MU.modS = MOut.ShearStressmagnitude;
MU.mU = mU;

%% Clean up the workspace
clearvars dT iFP iZ mH ncD ncS ncV t1Max t2Max

%% Calculate SN coordinates for the MU table
rMax = 200; %Maximum search distance
%Define the CL node spacing
ptSpace = 10;


clLength    =   max(CLProf.s);
nDiscr = floor(clLength/ptSpace);

transParam = [3 3 5 nDiscr rMax];
[snOut,~,clOut,~,~,iTrans] = xy2sn([CLProf.x,CLProf.y],[MU.x,MU.y],transParam);

MU.s = snOut(:,1);
MU.n = snOut(:,2);

figure
scatter(MU.x,MU.y,20,MU.s)

clearvars('filename','iTrans','nDiscr','rMax','transParam')

%% Bring in decimated depth map
%% Add image data to the MU table
% The morpho structure contains a number of images which are a culmination
% of the modeling analysis,they are:
% imgMU = morphologic unit map
% imgD  = model derived water depth
% imgNIR = band 34 (846 nm)
% imgSTD = variance of NIR band calculated using stdfilt function
% imgU  = modeled velocity
% imgS  = modeled shear stress
whos(rcM)

% next step:
% what can we do if we just have an image, field measurements and lidar
% calculate slope at all model points


%% Bring in image-derived depth
% This is the decimated depth map. The model output depth could be
% different depending on how well the model was fitting the actual water
% surface so lets stick to the image-derived depth.

ImgD = rcM.ImgD;
figure
scatter(ImgD.x,ImgD.y)
hold
scatter(MU.x,MU.y)
% just to make sure they are in fact different

%% Now create an interpolated surface from the decimated depth map

[X,Y] = meshgrid(imXdata,imYdata);

% Create a dmap surface from the decimated data
imgD = griddata(ImgD.x(:,1),ImgD.y(:,1),ImgD.depth(:,1),X,Y,'natural');
imgD(isnan(imgD)) = nan;

% Apply the water mask to remove the out of channel interpolation
 imgD  =   applyMask(imgD,rcI.waterMask);
 imgD(imgD == 0) = nan;
figure
imagesc(imgD)

clearvars X Y

%% Use impixel to extract depth values at model points
pix     =   impixel(imXdata,imYdata,flipud(imgD),MU.x(:,1),MU.y(:,1));
MU.imgD = pix(:,1);
figure
scatter(MU.x,MU.y,20,MU.imgD,'Filled')
figure
scatter(MU.x,MU.y,20,MU.modD,'Filled')
diff = MU.imgD - MU.modD;
figure
scatter(MU.x,MU.y,20,diff,'Filled')
colorbar


%% Add water surface slope to the mix
%% Fit a spline to the elevation profile
% Using the mean lidar at each centerline node calculated in
% SurfRfl_CASI2014_WSE.m
% This will allow us to calculate elevation as a function of s
[xData, yData] = prepareCurveData( CLProf.s,CLProf.meanZL);
% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.009602507575656453;

% Fit model to data.
[WSFit, WSFitGOF] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'WaterSurface' );
h = plot( WSFit, xData, yData );
legend( h, 'profZ vs. profS', 'WaterSurface', 'Location', 'NorthEast' );
% Label axes
xlabel profS
ylabel profZ
grid on    
    
clearvars('ft','WSFitGOF','h')   

%% Use the Model Fit to Calculate WSE at every depth point
   
   MU.wsZ = feval(WSFit,MU.s);

% Calculate slope    
  MU.wsSlope = -1*differentiate(WSFit,MU.s);
  
  
%% Bring in img data
whos(rcM)
% Use impixel to extract depth values at model points
pix     =   impixel(imXdata,imYdata,flipud(rcM.imgNIR),MU.x(:,1),MU.y(:,1));
MU.imgNIR = pix(:,1);

pix     =   impixel(imXdata,imYdata,flipud(rcM.imgSTD),MU.x(:,1),MU.y(:,1));
MU.imgSTD = pix(:,1);

%% Clean up the MU table
% remove model points where img depth is nan (outside the masked area)
iC = isnan(MU.imgD) == 1;
MU(iC,:) = [];

% looks like the best we can do with the casi image is 65% Pools are
% predicted at 80%, runs 48%, riffle 67%, 

%%% pulse_post3D processes data from MFiX simulations of 3D pulsating flow
%%% in volcanic eruptions.
%
% For a set of MFiX simulations, where <allruns> are the jobIDs of the
% simulations to be processed, pulse_post3D loads and manipulates MFiX
% output data to produce movies and plots of the time evolution of gas
% volume, gas temperature, gas pressure, and gas/solid momentum density
% over the model domain for the duration of each simulation. Output figures
% are saved to the same folder <run> that contains the relevant MFiX data
% files. 
%
% Special functions called: setCnsts3D; calc_inletFlow; volume3D;
%   entrainment3D; particleConc3D; gasTemperature3D; flowDensity3D
% Last edit: Taryn Black, 17 November 2015

clear all

%%% =================================================================== %%%
%%% ================= S E T  R U N  V A R I A B L E S ================= %%%

% ID numbers of MFiX runs to be processed:
  allruns = [888187];
  
% Set path. Must end in / & contain dirs titled by runIDs being processed.
  runpath = '/Users/taryn/OneDrive/Documents/MATLAB/MFIX_temp/';
%     runpath = '~/data2/rundata/';
  
% Set path for location of post-processing scripts (cannot change path file
% on Atlas cluster).
 postpath = '/Users/taryn/Documents/GitHub/pulse-post';
%    postpath = '~/data2/pulse-post';
  
% Choose whether to display ('on') or suppress ('off') figures.
% Note: vis must be 'off' when running remotely in -nodisplay mode.
  vis = 'on';

% Define isosurfaces for which gas volume fraction should be plotted
% <isoEPG>. Plumeedge is the gas volume fraction that defines the boundary
% of the plume. Set an RGB triple color (row in <colEPG>) and transparency
% <trnEPG> value for each isosurface. Note: the number of rows in <colEPG>
% and values in <trnEPG> must equal the number of isosurfaces.
  plumeedge = 1 - 1E-6;
  isoEPG = [plumeedge]; %,0.99995,0.9999,0.9995];
  colEPG = [0.5 0.5 0.5];   % gray
%           [1 0 0;
%             1 1 0;
%             0 1 1;
%             0 0 0];       
  trnEPG = [0.5];%[0.1,0.2,0.3,0.5];

% Number of labels to display on figures for each dimension
  xpoints = 5;      % horizontal axis 1
  ypoints = 9;      % vertical axis
  zpoints = 5;      % horizontal axis 2
  
% Scaling factor (meters*_fact) and unit for axes labels
  xfact = 1000;
  labelxunit = 'km';
  yfact = 1000;
  labelyunit = 'km';
  zfact = 1000;
  labelzunit = 'km';
  
% Entrainment animation colorbar limits (between -1 and 1)
  cmin = -0.5;
  cmax = 0.5;
  
% Viewing azimuth and elevation (in degrees) for 3D animations
  viewaz = -37.5;
  viewel = 20;
  
% Save animation timesteps as individual figures of specified filetype.
% Options: tif (*only if your photo viewer can read multipage
% tif, e.g. Windows Photo Viewer), jpeg, png, bmp, other image formats...
  imtype = 'tif';
  
% Other constants...
  Rgas = 461.5;     % gas constant
  
% Parameters for particle concentration slices figure. sdist* is the
% fraction along the *axis at which to cut the slice (between [] and 1).
  sdistX = 0.5; % scales to IMAX
  sdistY = [];  % scales to KMAX
  sdistZ = [];  % scales to JMAX (remember MFIX 'y' is MATLAB 'up'!)
  
% use in calc_inletFlow - calculate characteristic velocity, etc at which
% value for unsteady gas volume fraction?
  charEPG = 'mingas'; % 'maxgas' 'avggas'
  

%%% =================================================================== %%%
%%% =================================================================== %%%


% Display error if color/transparency doesn't match # of isosurfaces:
  [Crow,Ccol] = size(colEPG);
  if Crow ~= length(isoEPG)
      error('Error. \nNumber of rows in colEPG must match length of isoEPG.')
  elseif length(trnEPG) ~= length(isoEPG)
      error('Error. \nLength of trnEPG must match length of isoEPG.')
  end
  

for i = 1:length(allruns)
    
    run = allruns(i);
    sprintf('Now processing run #%d',run)

    dir = sprintf('%s%d',runpath,run);
        
  % Load and set constant simulation parameters
    ghostcells = 4;     % MFiX adds these to each domain dimension
    [IMAX,JMAX,KMAX,LENGTH,HEIGHT,WIDTH,RO_S1,RO_S2,RO_S3,NFR_S1,NFR_S2,...
        NFR_S3,PULSE,FREQ,MING,MAXG,VENT_R,DT,TSTOP,ATMOS,TROPO,BC_EPG,...
        BC_PG,BC_TG,BC_TS1,BC_TS2,BC_TS3] = setCnsts3D(run,dir,ghostcells);
    cd(postpath)

    timesteps = (TSTOP/DT)+1;
    
  % Define grid
%%% #TODO#: check to see if X/Y/Z are these even used
    X = LENGTH/(IMAX-ghostcells):LENGTH/(IMAX-ghostcells):LENGTH;
    Y = HEIGHT/(JMAX-ghostcells):HEIGHT/(JMAX-ghostcells):HEIGHT;
    Z = WIDTH/(KMAX-ghostcells):WIDTH/(KMAX-ghostcells):WIDTH;
  % Dimension resolution [meters/cell], excluding ghost cells
    XRES = LENGTH/(IMAX - ghostcells);
    YRES = HEIGHT/(JMAX - ghostcells);
    ZRES = WIDTH/(KMAX - ghostcells);
        
  % Set figure properties. Note: MFiX horizontal dimensions are LENGTH
  % (x; MATLAB X) and WIDTH (z; MATLAB Y), vertical dimension is HEIGHT
  % (y; MATLAB Z). tick* values are in meters.
    tickx = linspace(0,LENGTH,xpoints);
    labelx = tickx(2:end)/xfact;
    ticky = linspace(0,HEIGHT,ypoints);
    labely = ticky(2:end)/yfact;
    tickz = linspace(0,WIDTH,zpoints);
    labelz = tickz(2:end)/zfact;

    time = (0:timesteps-1)*DT;
      
  % Calculate characteristic inlet velocity for use in entrainment script.
    [XG,vel_char,MFR] = calc_inletFlow(charEPG,MING,MAXG,PULSE,BC_EPG,...
        BC_PG,BC_TG,Rgas,RO_S1,RO_S2,RO_S3,NFR_S1,NFR_S2,NFR_S3,BC_TS1,...
        BC_TS2,BC_TS3,VENT_R);
        
  % Manipulate data to time evolution over domain and save output
    vidEPG = volume3D(run,dir,vis,ghostcells,tickx,labelx,labelxunit,...
        ticky,labely,labelyunit,tickz,labelz,labelzunit,plumeedge,XRES,...
        YRES,ZRES,timesteps,IMAX,JMAX,KMAX,isoEPG,colEPG,trnEPG,time,...
        PULSE,FREQ,postpath);
    cd(postpath)

    vidEntr = entrainment3D(run,dir,vis,ghostcells,IMAX,JMAX,KMAX,...
        tickx,labelx,labelxunit,ticky,labely,labelyunit,tickz,labelz,...
        labelzunit,plumeedge,XRES,YRES,ZRES,postpath,PULSE,FREQ,time,...
        vel_char,cmin,cmax,viewaz,viewel,imtype);
    cd(postpath)

    vidPartConc = particleConc3D(run,dir,vis,IMAX,JMAX,KMAX,ghostcells,...
        tickx,labelx,labelxunit,ticky,labely,labelyunit,tickz,labelz,...
        labelzunit,XRES,YRES,ZRES,postpath,sdistX,sdistY,sdistZ);
    cd(postpath)
    
    vidGasTemp = gasTemperature3D(run,dir,vis,ghostcells,IMAX,JMAX,KMAX,...
        tickx,labelx,labelxunit,ticky,labely,labelyunit,tickz,labelz,...
        labelzunit,XRES,YRES,ZRES,sdistX,sdistY,sdistZ,postpath,ATMOS,...
        TROPO,Y);
    cd(postpath)
        
    vidFlowDens = flowDensity3D(run,dir,vis,IMAX,JMAX,KMAX,ghostcells,...
        postpath,Rgas,RO_S1,RO_S2,RO_S3,plumeedge,PULSE,FREQ,time,tickx,...
        labelx,labelxunit,ticky,labely,labelyunit,tickz,labelz,...
        labelzunit,XRES,YRES,ZRES,sdistX,sdistY,sdistZ);
    cd(postpath)

%     clearvars -except i allruns

end

%%% pulse_post3D processes output from simulations of 3D pulsating flow
%%% in volcanic eruptions generated from a modified MFiX code base.
%
% pulse_post3D loads and manipulates simulation output to produce movies,
% plots, and calculation output of the time evolution of gas volume,
% entrainment, particle concentration, gas temperature, flow density, and
% flow buoyancy over the model domain for the duration of each simulation.
% Figures are saved to the same folder <run> that contains the simulation
% output files.
%
% Functions called: setCnsts3D; calc_inletFlow; entrainment3D;
% particleConc3D; gasTemperature3D; flowDensity3D; velocity3D
%
% Last edit: Taryn Black, 11 April 2016

% clear all

%%% =================================================================== %%%
%%% ================== U S E R   P A R A M E T E R S ================== %%%

% ------------------- DEFINE SIMULATION IDS AND PATHS ------------------- %
% Name of run to be processed.
  run = 'testdata';
  
% Directory containing run data files to be processed.
%   runpath = sprintf('~/scratch/%s',run);
  runpath = 'C:/Users/taryn/Documents/GitHub/pulse-post/testdata';

% Directory where movies, images, and text files will be saved.
%   savepath = sprintf('~/data/ProductionRuns_Storage/%s/Figures',run);
  savepath = 'C:/Users/taryn/OneDrive/Documents/testdata_figs';
  
% Directory containing the suite of post-processing scripts.
%   postpath = '~/data/pulse-post';
  postpath = 'C:/Users/taryn/Documents/GitHub/pulse-post';
% ----------------------------------------------------------------------- %


% -------------- POST-PROCESSING DISPLAY AND SAVE SETTINGS -------------- %
% Choose whether to display ('on') or suppress ('off') figures.
% NOTE: must be 'off' when running remotely in -nodisplay mode.
  vis = 'on';

% Set end time (seconds) for movies. Use [] to process all timesteps.
  tstop = 300;
  
% Save animation timesteps as individual figures of specified filetype.
% Options: tif, jpeg, png, bmp, etc.
  imtype = 'tif';
% ----------------------------------------------------------------------- %
  

% --------------------- DEFINE 'PHYSICAL' CONSTANTS --------------------- %
% Gas volume fraction that defines the boundary of the plume
  plumeedge = 1 - 1E-6;
  
% Gas volume fraction at which to calculate characteristic (choked)
% velocity for the UNSTEADY case. Options: mingas, maxgas, avggas
  charEPG = 'mingas';
  
% Gas constant
  Rgas = 461.5;
  
% Gas viscosity [Pa.s]
  MU_G0 = 1E-5;
  
% Gravitational acceleration [m/s2]
  g = 9.81;
% ----------------------------------------------------------------------- %


% ----------------- PHYSICAL PROPERTY DISPLAY SETTINGS ------------------ %
% Define isosurfaces <isoEPG> for which gas volume fraction should be
% plotted. Set an RGB triple color (row in <colEPG>) and transparency
% <trnEPG> value for each isosurface. 
% NOTE: number of rows in <colEPG> and values in <trnEPG> must equal the
% number of isosurfaces.
  isoEPG = [plumeedge];
  colEPG = [0.7 0.7 0.7];   
  trnEPG = [0.7];
  
% Colorbar limits
  % Entrainment (must be between -1 and 1)
    entrainment_crange = [-0.5 0.5];
  % Particle concentration (log scale; must be less than 0)
    particleConc_crange = [-6 -3];
  % Flow density [kg/m3]
    flowDensity_crange = [0 4];
  % Relative density of flow (atmos_density - flow_density)
    flowBuoyancy_crange = [-1 1];
  % Gas temperature [K] (max is defined by vent inlet temperature)
    gasTemperature_cmin = 300; 
  % Velocity magnitude [m/s]
    velocity_crange = [0 300];
  % Vorticity
    vorticity_crange = [-20 20];
  % Mass flux [kg/m2.s]
    massflux_crange = [-1E5 1E5];
  
% Altitudes [meters] at which to calculate vertical mass flux
% To include script-calculated jet height, include one NaN in vector
  massflux_alts = [0 2000 4000 NaN];
    
% Slice distance and direction for 3D-slice figures. sdist* is the fraction
% along the *axis at which to cut the slice (between 0 and 1).
% NOTE: can only slice along one axis; use [] for other two axes.
  sdistX = 0.5; % scales to IMAX
  sdistY = [];  % scales to KMAX
  sdistZ = [];  % scales to JMAX (remember MFIX 'Y' is MATLAB 'z' (up)!)
% ----------------------------------------------------------------------- %


% --------------------- FIGURE PROPERTIES SETTINGS ---------------------- %
% Viewing azimuth and elevation (in degrees) for 3D animations
  viewaz = -37.5;
  viewel = 20;
  
% Number of tick labels to display on axes for each dimension
  Xpoints = 5;      % horizontal axis 1
  Ypoints = 11;      % vertical axis
  Zpoints = 5;      % horizontal axis 2

% Scaling factor (meters*_fact) and unit for axes labels
  Xfact = 1000;     % horizontal axis 1
  labelXunit = 'km';
  Yfact = 1000;     % vertical axis
  labelYunit = 'km';
  Zfact = 1000;     % horizontal axis 2
  labelZunit = 'km';
% ----------------------------------------------------------------------- %   


% ---------------------- FILE READ TYPES AND NAMES ---------------------- %
% readtype == 1 :  file was originally written in binary, and was converted
%   to ascii using convert_Pulse.f90. Name prefix is {*}_t%02d.txt
% readtype == 2 :  file was originally written in ascii and exists as a
%   single file instead of pre-separated timesteps. Name is {originalname}.
% fname* format : fname = {'prefix for readtype=1' 'name for readtype=2'}
  readEPG  = 1;
  fnameEPG = {'EP' 'EP_Ga'};
  
  readEPS1  = 1;
  fnameEPS1 = {'EP' ''};
  readEPS2  = 1;
  fnameEPS2 = {'EP' ''};
  readEPS3  = 1;
  fnameEPS3 = {'EP' ''};
  
  readTG   = 1;
  fnameTG  = {'T_G' 'T_Ga'};
  
  readROG  = 1;
  fnameROG = {'Current_Density' 'RO_Ga'};
  
  readUG  = 2;
  fnameUG = {'U_G' 'U_Ga'};
  readVG  = 1;
  fnameVG = {'U_G' 'V_Ga'};
  readWG  = 2;
  fnameWG = {'U_G' 'W_Ga'};
  
  readVS1  = 1;
  fnameVS1 = {'V_S' 'V_S1a'};
  readVS2  = 1;
  fnameVS2 = {'V_S' 'V_S2a'};
  readVS3  = 1;
  fnameVS3 = {'V_S' 'V_S3a'};
% ----------------------------------------------------------------------- %


%%% ============== E N D   U S E R   P A R A M E T E R S ============== %%%
%%% =================================================================== %%%


%%% CHECK FOR PROBLEMS BEFORE PROCESSING
% Display error if color/transparency doesn't match # of isosurfaces:
  [Crow,Ccol] = size(colEPG);
  if Crow ~= length(isoEPG)
      error('Error: Number of rows in colEPG must match length of isoEPG.')
  elseif length(trnEPG) ~= length(isoEPG)
      error('Error: Length of trnEPG must match length of isoEPG.')
  end
  
  
%%% GENERATE FIGURES AND DO CALCULATIONS FOR RUN    
  fprintf('Processing run %s.\n',run)
  titlerun = strrep(run,'_','\_');
        
% Define post-processing constants from simulation parameters
  ghostcells = 4;     % MFiX adds these to each domain dimension
  [IMAX,JMAX,KMAX,LENGTH,HEIGHT,WIDTH,RO_S1,RO_S2,RO_S3,NFR_S1,...
      NFR_S2,NFR_S3,PULSE,FREQ,MING,MAXG,VENT_R,DT,TSTOP,ATMOS,...
      TROPO,BC_EPG,BC_PG,BC_TG,BC_TS1,BC_TS2,BC_TS3,D_S1,D_S2,D_S3] ...
      = setCnsts3D(run,runpath,ghostcells,tstop);
  cd(postpath)
  timesteps = TSTOP/DT;
  time = (0:timesteps)*DT;
    
% Define dimensional spatial resolution [meters/cell] and domain grid
  XRES = LENGTH/(IMAX - ghostcells);
  YRES = HEIGHT/(JMAX - ghostcells);
  ZRES = WIDTH/(KMAX - ghostcells);
  XGRID = XRES:XRES:LENGTH;
  YGRID = YRES:YRES:HEIGHT;
  ZGRID = ZRES:ZRES:WIDTH;
    
% Define tick locations and labels
  tickx = linspace(0,LENGTH,Xpoints);
  labelx = tickx(2:end)/Xfact;
  ticky = linspace(0,HEIGHT,Ypoints);
  labely = ticky(2:end)/Yfact;
  tickz = linspace(0,WIDTH,Zpoints);
  labelz = tickz(2:end)/Zfact;
  
% Calculate characteristic inlet velocity
  [XG,vel_char,MFR,MASSFLUX,MFR_SOL,MASSFLUX_SOL,jetheight] = ...
      calc_inletFlow(charEPG,MING,MAXG,PULSE,BC_EPG,BC_PG,BC_TG,Rgas,...
      RO_S1,RO_S2,RO_S3,NFR_S1,NFR_S2,NFR_S3,BC_TS1,BC_TS2,BC_TS3,VENT_R,g);
  jetheight = jetheight/YRES;
  
% Calculate particle Stokes numbers
  [Stokes_S1,Stokes_S2,Stokes_S3] = calcStokes(RO_S1,RO_S2,RO_S3,D_S1,...
      D_S2,D_S3,vel_char,MU_G0,VENT_R,PULSE,FREQ);
  dlmwrite(fullfile(savepath,'particleStokes.txt'),[Stokes_S1 Stokes_S2 ...
      Stokes_S3],'delimiter','\t');
  
% Redefine mass flux altitudes as vector indices
  massflux_alts(massflux_alts==0) = YRES;
  massflux_alts(isnan(massflux_alts)) = round(jetheight)*YRES;
  massflux_alts = sort(massflux_alts);
  massflux_legend = strcat(strsplit(num2str(massflux_alts/Yfact)),...
      {' '},cellstr(labelYunit));
  massflux_alts = massflux_alts./YRES;
  
% Evolving volume fraction calculations and figure
  [EPG_vent,EPS_vent] = evolveVolumeFraction(MING,MAXG,FREQ,PULSE,DT,...
      TSTOP,vis,titlerun,imtype,time,savepath,postpath,run);
   
% Entrainment and gas volume fraction calculations and figures
  entrainment3D(run,runpath,vis,ghostcells,IMAX,JMAX,KMAX,tickx,labelx,...
      labelXunit,ticky,labely,labelYunit,tickz,labelz,labelZunit,...
      plumeedge,XRES,YRES,ZRES,postpath,PULSE,FREQ,time,vel_char,...
      jetheight,entrainment_crange,viewaz,viewel,imtype,titlerun,...
      timesteps,isoEPG,colEPG,trnEPG,DT,VENT_R,savepath,...
     readEPG,fnameEPG,readUG,fnameUG,readVG,fnameVG,readWG,fnameWG);
  cd(postpath)
  
% Particle concentration calculations and figures
  particleConc3D(run,runpath,vis,IMAX,JMAX,KMAX,ghostcells,tickx,labelx,...
      labelXunit,ticky,labely,labelYunit,tickz,labelz,labelZunit,...
      XRES,YRES,ZRES,postpath,sdistX,sdistY,sdistZ,...
      particleConc_crange,titlerun,timesteps,PULSE,FREQ,time,imtype,...
      plumeedge,D_S1,D_S2,D_S3,savepath,readEPS1,fnameEPS1,...
      readEPS2,fnameEPS2,readEPS3,fnameEPS3,readEPG,fnameEPG);
  cd(postpath)

% Gas temperature calculations and figures
  gasTemperature3D(run,runpath,vis,ghostcells,IMAX,JMAX,KMAX,tickx,...
      labelx,labelXunit,ticky,labely,labelYunit,tickz,labelz,...
      labelZunit,XRES,YRES,ZRES,sdistX,sdistY,sdistZ,postpath,ATMOS,...
      TROPO,YGRID,BC_TG,gasTemperature_cmin,PULSE,FREQ,time,titlerun,...
      timesteps,imtype,plumeedge,savepath,readTG,fnameTG,readEPG,fnameEPG);
  cd(postpath)
 
% Density and buoyancy calculations and figures  
  flowDensity3D(run,runpath,vis,IMAX,JMAX,KMAX,ghostcells,postpath,...
      RO_S1,RO_S2,RO_S3,plumeedge,PULSE,FREQ,time,tickx,labelx,...
      labelXunit,ticky,labely,labelYunit,tickz,labelz,labelZunit,...
      XRES,YRES,ZRES,sdistX,sdistY,sdistZ,flowDensity_crange,titlerun,...
      flowBuoyancy_crange,timesteps,imtype,savepath,readEPG,fnameEPG,...
      readROG,fnameROG,readEPS1,fnameEPS1,readEPS2,fnameEPS2,readEPS3,...
      fnameEPS3,jetheight);
  cd(postpath)

% Velocity magnitude calculations and figures
  velocity3D( runpath,sdistX,sdistY,sdistZ,vis,run,...
      timesteps,postpath,IMAX,JMAX,KMAX,ghostcells,velocity_crange,...
      PULSE,time,titlerun,FREQ,tickx,XRES,labelx,labelXunit,...
      ticky,YRES,labely,labelYunit,tickz,ZRES,labelz,labelZunit,...
      imtype,plumeedge,viewaz,viewel,YGRID,vorticity_crange,...
      savepath,readEPG,fnameEPG,readUG,fnameUG,readVG,fnameVG,...
      readWG,fnameWG);
  cd(postpath)
        
% Mass flux calculations and figures
  massFlux3D(runpath,vis,viewaz,viewel,ghostcells,IMAX,JMAX,KMAX,...
      tickx,ticky,tickz,XRES,YRES,ZRES,labelx,labely,labelz,...
      labelXunit,labelYunit,labelZunit,run,timesteps,postpath,...
      massflux_alts,RO_S1,RO_S2,RO_S3,plumeedge,massflux_crange,...
      PULSE,FREQ,time,titlerun,massflux_legend,imtype,savepath,readEPG,...
      fnameEPG,readVS1,fnameVS1,readVS2,fnameVS2,readVS3,fnameVS3,...
      readEPS1,fnameEPS1,readEPS2,fnameEPS2,readEPS3,fnameEPS3,...
      jetheight,MASSFLUX_SOL,sdistX,sdistY,sdistZ);
  cd(postpath)
  
  close all

  disp('* ================================================= *')
  disp('   P O S T - P R O C E S S I N G   C O M P L E T E')
  disp('* ================================================= *')

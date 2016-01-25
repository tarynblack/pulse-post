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
% Last edit: Taryn Black, 22 January 2015

clear all

%%% =================================================================== %%%
%%% ================== U S E R   P A R A M E T E R S ================== %%%

% ------------------- DEFINE SIMULATION IDS AND PATHS ------------------- %
% Names of runs to be processed.
%   runIDs = {'F_1_998_9999'};
  runIDs = {'testdata'};
  
% Path for directory containing post-processing data directories (runIDs).
%   runpath = '~/data2/rundata/';
  runpath = 'C:/Users/taryn/Documents/GitHub/pulse-post/';
  
% Path for location of post-processing scripts.
%   postpath = '~/data2/pulse-post';
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
    entrainment_cmin = -0.5;
    entrainment_cmax = 0.5;
  % Particle concentration
    particleConc_cmin = -5;
    particleConc_cmax = 10;
  % Flow density [kg/m3]
    flowDensity_cmin = 0;
    flowDensity_cmax = 8;
  % Relative density of flow (atmos_density - flow_density)
    flowBuoyancy_cmin = -3;
    flowBuoyancy_cmax = 1;
  % Gas temperature [K] (max is defined by vent inlet temperature)
    gasTemperature_cmin = 300; 
  % Velocity magnitude [m/s]
    velocity_cmin = 0;
    velocity_cmax = 300;
    
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
  Ypoints = 9;      % vertical axis
  Zpoints = 5;      % horizontal axis 2

% Scaling factor (meters*_fact) and unit for axes labels
  Xfact = 1000;     % horizontal axis 1
  labelXunit = 'km';
  Yfact = 1000;     % vertical axis
  labelYunit = 'km';
  Zfact = 1000;     % horizontal axis 2
  labelZunit = 'km';
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
  
  
%%% GENERATE FIGURES AND DO CALCULATIONS FOR EACH RUN
for i = 1:length(runIDs)
    
      run = runIDs{i};
      disp(sprintf('Now processing run %s',run))
      dir = sprintf('%s%s',runpath,run);
      titlerun = strrep(run,'_','\_');
        
    % Define post-processing constants from simulation parameters
      ghostcells = 4;     % MFiX adds these to each domain dimension
      [IMAX,JMAX,KMAX,LENGTH,HEIGHT,WIDTH,RO_S1,RO_S2,RO_S3,NFR_S1,...
          NFR_S2,NFR_S3,PULSE,FREQ,MING,MAXG,VENT_R,DT,TSTOP,ATMOS,...
          TROPO,BC_EPG,BC_PG,BC_TG,BC_TS1,BC_TS2,BC_TS3] ...
          = setCnsts3D(run,dir,ghostcells,tstop);
      cd(postpath)
      timesteps = TSTOP/DT;
      time = (0:timesteps)*DT;
    
    % Define domain grid and dimensional spatial resolution [meters/cell]
      X = LENGTH/(IMAX-ghostcells):LENGTH/(IMAX-ghostcells):LENGTH;
      Y = HEIGHT/(JMAX-ghostcells):HEIGHT/(JMAX-ghostcells):HEIGHT;
      Z = WIDTH/(KMAX-ghostcells):WIDTH/(KMAX-ghostcells):WIDTH;
      XRES = LENGTH/(IMAX - ghostcells);
      YRES = HEIGHT/(JMAX - ghostcells);
      ZRES = WIDTH/(KMAX - ghostcells);
        
    % Define tick locations and labels
      tickx = linspace(0,LENGTH,Xpoints);
      labelx = tickx(2:end)/Xfact;
      ticky = linspace(0,HEIGHT,Ypoints);
      labely = ticky(2:end)/Yfact;
      tickz = linspace(0,WIDTH,Zpoints);
      labelz = tickz(2:end)/Zfact;
      
    % Calculate characteristic inlet velocity
      [XG,vel_char,MFR] = calc_inletFlow(charEPG,MING,MAXG,PULSE,BC_EPG,...
          BC_PG,BC_TG,Rgas,RO_S1,RO_S2,RO_S3,NFR_S1,NFR_S2,NFR_S3,...
          BC_TS1,BC_TS2,BC_TS3,VENT_R);
        
    % Entrainment and gas volume fraction calculations and figures
      entrainment3D(run,dir,vis,ghostcells,IMAX,JMAX,KMAX,tickx,labelx,...
          labelXunit,ticky,labely,labelYunit,tickz,labelz,labelZunit,...
          plumeedge,XRES,YRES,ZRES,postpath,PULSE,FREQ,time,vel_char,...
          entrainment_cmin,entrainment_cmax,viewaz,viewel,imtype,...
          titlerun,timesteps,isoEPG,colEPG,trnEPG,DT,VENT_R);
      cd(postpath)
      
    % Particle concentration calculations and figures
      particleConc3D(run,dir,vis,IMAX,JMAX,KMAX,ghostcells,tickx,labelx,...
          labelXunit,ticky,labely,labelYunit,tickz,labelz,labelZunit,...
          XRES,YRES,ZRES,postpath,sdistX,sdistY,sdistZ,...
          particleConc_cmin,particleConc_cmax,titlerun,timesteps,PULSE,...
          FREQ,time,imtype,plumeedge);
      cd(postpath)
    
    % Gas temperature calculations and figures
      gasTemperature3D(run,dir,vis,ghostcells,IMAX,JMAX,KMAX,tickx,...
          labelx,labelXunit,ticky,labely,labelYunit,tickz,labelz,...
          labelZunit,XRES,YRES,ZRES,sdistX,sdistY,sdistZ,postpath,ATMOS,...
          TROPO,Y,BC_TG,gasTemperature_cmin,PULSE,FREQ,time,titlerun,...
          timesteps,imtype,plumeedge);
      cd(postpath)
     
    % Density and buoyancy calculations and figures  
      flowDensity3D(run,dir,vis,IMAX,JMAX,KMAX,ghostcells,postpath,...
          RO_S1,RO_S2,RO_S3,plumeedge,PULSE,FREQ,time,tickx,labelx,...
          labelXunit,ticky,labely,labelYunit,tickz,labelz,labelZunit,...
          XRES,YRES,ZRES,sdistX,sdistY,sdistZ,flowDensity_cmin,...
          flowDensity_cmax,titlerun,flowBuoyancy_cmin,flowBuoyancy_cmax,...
          timesteps,imtype);
      cd(postpath)

    % Velocity magnitude calculations and figures
      velocity3D( dir,sdistX,sdistY,sdistZ,vis,run,...
          timesteps,postpath,IMAX,JMAX,KMAX,ghostcells,velocity_cmin,...
          velocity_cmax,PULSE,time,titlerun,FREQ,tickx,XRES,labelx,labelXunit,...
          ticky,YRES,labely,labelYunit,tickz,ZRES,labelz,labelZunit,imtype,plumeedge);
      cd(postpath)
      
      close all
      clearvars -except i allruns
      disp('* ================================================= *')
      disp('   P O S T - P R O C E S S I N G   C O M P L E T E')
      disp('* ================================================= *')

end

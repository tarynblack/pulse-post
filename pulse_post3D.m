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
% Special functions called: mfixData3D; setCnsts3D; volume3D
% Last edit: Taryn Black, 12 October 2015

clear all

%%% =================================================================== %%%
%%% ================= S E T  R U N  V A R I A B L E S ================= %%%

% ID numbers of MFiX runs to be processed:
  allruns = [734662 734663 734664];
  
% Set path. Must end in / & contain dirs titled by runIDs being processed.
%   runpath = '/Users/Taryn/Documents/MATLAB/MFIX_temp/'; 
  runpath = '~/data2/rundata/';
  
% Set path for location of post-processing files (cannot change path file
% on Atlas cluster).
%   postpath = '/Users/Taryn/Documents/GitHub/pulse-post';
  postpath = '~/data2/pulse-post';
  
% Choose whether to display ('on') or suppress ('off') figures.
% Note: vis must be 'off' when running remotely in -nojvm mode.
  vis = 'off';
  
% Choose whether to load and process a variable (1) or skip it (0).
  onE = 1;      % EP_G
  onP = 0;      % P_G
  onT = 0;      % T_G, T_S1, T_S2
  onV = 0;      % U_G, U_S1, U_S2, V_G, V_S1, V_S2, W_G, W_S2, W_S2
  onR = 0;      % RO_G, ROP_S1, ROP_S2
  onX = 0;      % X_G2

% Define isosurfaces for which gas volume fraction should be plotted
% <isoEPG>. Set an RGB triple color (row in <colEPG>) and transparency
% <trnEPG> value for each isosurface. Note: the number of rows in <colEPG>
% and values in <trnEPG> must equal the number of isosurfaces.
  isoEPG = [0.99999,0.99995,0.9999,0.9995];
  colEPG = [1 0 0;
            1 1 0;
            0 1 1;
            0 0 0];       
  trnEPG = [0.1,0.2,0.3,0.5];

% Number of labels to display on figures for each dimension
  xpoints = 5;      % horizontal axis 1
  ypoints = 9;      % vertical axis
  zpoints = 5;      % horizontal axis 2

% TODO: check that these are necessary and automate in mfixconst
ventwidth = 400;    % meters
T_range = [300 1100];
P_range = [5E4 5E5];
BD_range = [0 7];
V_range = [0 400];

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
   
    % Import data generated in MFiX
%       [  EP_G,...
%           P_G,...
%           T_G,T_S1,T_S2,...
%           U_G,U_S1,U_S2,...
%           V_G,V_S1,V_S2,...
%           W_G,W_S1,W_S2,...
%           RO_G,ROP_S1,ROP_S2,...
%           X_G2  ] = mfixData3D(run,dir,onE,onP,onT,onV,onR,onX);
      cd(postpath)
        
    % Load and set constant simulation parameters
      ghostcells = 4;     % MFiX adds these to each domain dimension
      [IMAX,JMAX,KMAX,LENGTH,HEIGHT,WIDTH,RO_S1,RO_S2,PULSE,FREQ,DT] ...
          = setCnsts3D(run,dir);
      cd(postpath)
%       timesteps = length(EP_G)/(IMAX*JMAX*KMAX);
    % TODO: add TSTOP to mfixconst
    timesteps = (600/DT)+1;
      % Cell edges? [meters]
      % TODO: check to see if X/Y/Z need ghostcells or not in calcs, where
      % are these even used?
        X = LENGTH/IMAX:LENGTH/IMAX:LENGTH;
        Y = HEIGHT/JMAX:HEIGHT/JMAX:HEIGHT;
        Z = WIDTH/KMAX:WIDTH/KMAX:WIDTH;
      % Dimension resolution [meters/cell], excluding ghost cells
        XRES = LENGTH/(IMAX - ghostcells);
        YRES = HEIGHT/(JMAX - ghostcells);
        ZRES = WIDTH/(KMAX - ghostcells);
        
    % Set figure properties. Note: MFiX horizontal dimensions are LENGTH
    % (x; MATLAB x) and WIDTH (z; MATLAB y), vertical dimension is HEIGHT
    % (y; MATLAB z).
      xkilolabels = linspace(0,LENGTH/1000,xpoints);
      ykilolabels = linspace(0,HEIGHT/1000,ypoints);
      zkilolabels = linspace(0,WIDTH/1000,zpoints);
      time = (0:timesteps-1)*DT;
        
    % Manipulate data to time evolution over domain and save output
%       if onE == 1
%           vidEPG = volume3D(run,dir,vis,ghostcells,xkilolabels,...
%               ykilolabels,zkilolabels,timesteps,EP_G,IMAX,JMAX,KMAX,...
%               isoEPG,colEPG,trnEPG,time,PULSE,FREQ,postpath);
%       end
    vidEPG = volume3D(run,dir,vis,ghostcells,xkilolabels,...
              ykilolabels,zkilolabels,timesteps,IMAX,JMAX,KMAX,...
              isoEPG,colEPG,trnEPG,time,PULSE,FREQ,postpath);
      cd(postpath)

        
%     clearvars -except i allruns

end

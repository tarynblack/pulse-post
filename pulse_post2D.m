% pulse_post2D processes data from MFIX pulsating flow simulations.
%   For a specified MFIX simulation, where <run> is the simulation jobID,
%   pulse_post loads and manipulates MFIX output data to produce movies and
%   plots of the time evolution of gas volume, gas temperature, gas
%   pressure, and gas/solid momentum density over the model domain for the
%   duration of the simulation. Output figures are saved to the same folder
%   <run> that contains the relevant MFIX data files.
%   Functions called: mfixData, setCnsts, volume, temperature, pressure,
%   momentum
%   Taryn Black, last edit 25 April 2015

% clear all

allruns = 169101;
% allruns = [169065 169066 169067 169068 169069 169070 169071 169072 169076 169079 169080 169081 169082 169083 169084 169085 169086 169087 169092 169093 169094 169096 169097 169098 169099 169100 169101];

for i = 1:length(allruns)

%%% ================ SET RUN VARIABLES ================ %%%
run = allruns(i)
h = [0 100 500];    % heights to sample momentum, in m above vent
ventwidth = 100;
dt = 3;             % write-out timestep [s]
T_range = [300 900];
P_range = [5E4 5E5];
BD_range = [0 7];
V_range = [0 400];
%%% =================================================== %%%


% cd(sprintf('%d',run))

%%% Import data generated in MFIX
    [EP_G,RO_G,ROP_S1,ROP_S2,V_G,U_G,V_S1,V_S2,T_G,P_G] = mfixData(run);

%%% Load and set simulations properties
    [IMAX,JMAX,HEIGHT,LENGTH,RO_S1,RO_S2,PULSE,FREQ] = setCnsts(run);
    X = LENGTH/IMAX:LENGTH/IMAX:LENGTH;
    Y = HEIGHT/JMAX:HEIGHT/JMAX:HEIGHT;
    XRES = LENGTH/(IMAX-2);     % x-resolution, m per cell (exc. ghost cells)
    YRES = HEIGHT/(JMAX-2);     % y-resolution, m per cell (exc. ghost cells)
    timesteps = length(EP_G)/(IMAX*JMAX);
    
%%% Set figure properties
    distance = linspace(0,LENGTH/1000,5);
    altitude = linspace(0,HEIGHT/1000,9);
    time = (0:timesteps-1)*dt;
    
%%% Manipulate data to time evolution over domain and save output
    mGV = volume2D(EP_G,timesteps,IMAX,JMAX,X,Y,distance,altitude,run,PULSE,FREQ,time);
        movie2avi(mGV,sprintf('gVol_%d.avi',run));
        cd ..
%     mGT = temperature2D(T_G,timesteps,IMAX,JMAX,X,Y,distance,altitude,run,T_range,PULSE,FREQ,time);
%         movie2avi(mGT,sprintf('gTemp_%d.avi',run));
%         cd ..
%     mGP = pressure2D(P_G,timesteps,IMAX,JMAX,X,Y,distance,altitude,run,P_range,PULSE,FREQ,time);
%         movie2avi(mGP,sprintf('gPres_%d.avi',run));
%         cd ..
%     [mBD,ROB] = bulkdens2D(EP_G,P_G,T_G,RO_S1,RO_S2,timesteps,IMAX,JMAX,X,Y,distance,altitude,run,BD_range,PULSE,FREQ,time);
%         movie2avi(mBD,sprintf('bulkDens_%d.avi',run));
%         cd ..
%     mGVel = velocity2D(V_G,U_G,timesteps,IMAX,JMAX,X,Y,distance,altitude,run,V_range,PULSE,FREQ,time);
%         movie2avi(mGVel,sprintf('gVelo_%d.avi',run));
%         cd ..
%     [mEntr,hFig2,hFig3,avg_entrain,plumevolume,vol_inlet] = entrainment2D(EP_G,V_G,U_G,run,timesteps,distance,altitude,XRES,YRES,IMAX,JMAX,PULSE,FREQ,time);
%         movie2avi(mEntr,sprintf('entrain_%d.avi',run),'Compression','None');
%         saveas(hFig2,sprintf('AvgEntrCoeff_%d',run),'png');
%         saveas(hFig3,sprintf('VolumeEntrd_%d',run),'png');
%         save(sprintf('avg_entrain_%d.txt',run),'avg_entrain','-ascii')
%         save(sprintf('plumevolume_%d.txt',run),'plumevolume','-ascii')
%         save(sprintf('vol_inlet_%d.txt',run),'vol_inlet','-ascii')
%         close all
%         cd ..
%     [h1,h2,h3] = momentum2D(h,EP_G,RO_G,ROP_S1,ROP_S2,V_G,V_S1,V_S2,IMAX,JMAX,LENGTH,XRES,YRES,ventwidth,timesteps,run);
%         saveas(h1,sprintf('gMomSep_%d',run),'jpg'); 
%         saveas(h2,sprintf('sMom_%d',run),'jpg');
%         saveas(h3,sprintf('gMomStack_%d',run),'jpg')
clearvars -except i allruns
end

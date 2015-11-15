function [ vidFlowDens ] = flowDensity3D( run,dir,vis,IMAX,JMAX,KMAX,...
    ghostcells,postpath,RO_S1,RO_S2,RO_S3 )
%flowDensity3D calculates the net density of the flow from gas and particle
%densities and volume fractions.
%   Detailed explanation goes here
%
%   Special functions called: varchunk3D
%   Last edit: Taryn Black, 14 November 2015

    varname = 'Flow density';
    
%%% Initialize figure frames
    % #TODO# decide how to display (isosurf or slice?)
    
%%% Initialize video
    vidFlowDens = VideoWriter(sprintf('vidFlowDens_%d.avi',run));
    vidFlowDens.Quality = 100;
    vidFlowDens.FrameRate = 10;
    open(vidFlowDens);
    set(gcf,'Visible',vis);

%%% Calculate flow density at each timestep, plot, and save video.
    fID_EPG = fopen(sprintf('EP_G_%d',run));
    fID_PG  = fopen(sprintf('P_G_%d',run));
    fID_TG  = fopen(sprintf('T_G_%d',run));
    fID_RS1 = fopen(sprintf('ROP_S1_%d',run));
    fID_RS2 = fopen(sprintf('ROP_S2_%d',run));
    fID_RS3 = fopen(sprintf('ROP_S3_%d',run));
    
    t = 0;
    while ~feof(fID_EPG)
        
        t = t+1;
        
        cd(postpath)
        try
            EPG = varchunk3D(fID_EPG,IMAX,JMAX,KMAX,ghostcells);
        catch ME
            warning('Error in varchunk3D at t=%d s:\n%s\nContinuing to next simulation.',time(t-1),ME.identifier)
            break
        end
        PG  = varchunk(fID_PG,IMAX,JMAX,KMAX,ghostcells);
        TG  = varchunk(fID_TG,IMAX,JMAX,KMAX,ghostcells);
        RS1 = varchunk(fID_RS1,IMAX,JMAX,KMAX,ghostcells);
        RS2 = varchunk(fID_RS2,IMAX,JMAX,KMAX,ghostcells);
        RS3 = varchunk(fID_RS3,IMAX,JMAX,KMAX,ghostcells);
        
      % Skip processing for initial timestep
        if t==1;
            continue
        end        
        cla;
        
        flowdensity = (EPG*PG/(R*TG)) + ((1-EPG)*(RS1*RO_S1 + RS2*RO_S2 + RS3*RO_S3));
    
    
    end

    cd(dir)
    close(vidFlowDens);
    cd(postpath)
    sprintf('Flow density processing complete. \nvidFlowDens_%d has been saved to %s',run,dir)
    
end


function [ vidFlowDens ] = flowDensity3D( run,dir,vis,IMAX,JMAX,KMAX,...
    ghostcells,postpath,RO_S1,RO_S2,RO_S3,plumeedge,PULSE,FREQ,time,...
    tickx,labelx,labelxunit,ticky,labely,labelyunit,tickz,labelz,...
    labelzunit,XRES,YRES,ZRES,sdistX,sdistY,sdistZ )
%flowDensity3D calculates the net density of the flow from gas and particle
%densities and volume fractions.
%   Detailed explanation goes here
%
%   Special functions called: varchunk3D, pulsetitle
%   Last edit: Taryn Black, 14 November 2015

    varname = 'Flow density';

%%% #TODO# decide whether or not to keep this - just using it temporarly to
%%% visualize buoyancy slice.
%%% Ensure that 'no slice' directions are empty and determine figure
%%% viewing angle based on slice direction
    if sdistX==0
        sdistX = [];
    end
    if sdistY==0
        sdistY = [];
    end
    if sdistZ==0
        sdistZ = [];
    end
    
    if isempty(sdistX) && isempty(sdistY)
        saz = 0;
        sel = 0;
    elseif isempty(sdistY) && isempty(sdistZ)
        saz = 90;
        sel = 0;
    elseif isempty(sdistX) && isempty(sdistZ)
        saz = 0;
        sel = 90;
    else [saz,sel] = view(3);
    end
    
%%% Initialize figure frames 
%%% #TODO# display as density slice. Below is temporary to visualize buoyancy slice.
    cd(dir)
    fig = figure('Name','Plume Buoyancy','visible',vis);
    hold on
    view(saz,sel)
    axis equal
    axis([ghostcells-1,IMAX-(ghostcells/2),ghostcells-1,...
        KMAX-(ghostcells/2),ghostcells-1,JMAX-(ghostcells/2)]);
    set(gca,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,'FontSize',12)
    xlabel(sprintf('\\bf Distance_x (%s)',labelxunit),'FontSize',12)
    set(gca,'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,'FontSize',12)
    ylabel(sprintf('\\bf Distance_z (%s)',labelzunit),'FontSize',12)
    set(gca,'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely,'FontSize',12)
    zlabel(sprintf('\\bf Altitude (%s)',labelyunit),'FontSize',12)
    
    
%%% Initialize video
    vidFlowDens = VideoWriter(sprintf('vidFlowDens_%d.avi',run));
    vidFlowDens.Quality = 100;
    vidFlowDens.FrameRate = 10;
    open(vidFlowDens);
    set(gcf,'Visible',vis);

%%% Calculate flow density at each timestep, plot, and save video.
    fID_EPG = fopen(sprintf('EP_G_%d',run));
    fID_ROG  = fopen(sprintf('RO_G_%d',run));
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
        ROG  = varchunk3D(fID_ROG,IMAX,JMAX,KMAX,ghostcells);
        RS1 = varchunk3D(fID_RS1,IMAX,JMAX,KMAX,ghostcells);
        RS2 = varchunk3D(fID_RS2,IMAX,JMAX,KMAX,ghostcells);
        RS3 = varchunk3D(fID_RS3,IMAX,JMAX,KMAX,ghostcells);
        
      % Skip processing for initial timestep
        if t==1;
            continue
        end        
        cla;
        
      % Density at every point in domain: sum of gas and particle densities
%%% #TODO# change calculation: when using Mary's data, RS1 --> RS1*RO_S1
        domaindensity = (EPG.*ROG) + (RS1 + RS2 + RS3);
        
      % Determine whether flow is buoyant (less dense than atmosphere)
      % using mean flow/atmospheric densities at each altitude slice.
        inplume = EPG <= plumeedge;
        inatmos = EPG > plumeedge;
        
        flowdensity = domaindensity.*inplume;
        atmsdensity = domaindensity.*inatmos;
        
        for i = 1:(JMAX-ghostcells)
            z_flowd = flowdensity(:,:,i);
            z_atmsd = atmsdensity(:,:,i);
            avgflowdens_altitude(i) = mean(z_flowd(z_flowd~=0));
            avgatmsdens_altitude(i) = mean(z_atmsd(z_atmsd~=0));
        end
        
        % Relative density of atmosphere to flow. Positive = flow is
        % positively buoyant (less dense than atmosphere). Negative = flow
        % is negatively buoyant (more dense than atmosphere). Value =
        % difference in density
          reldensity = avgatmsdens_altitude - avgflowdens_altitude;
          reldens_3D = repmat(reldensity',1,KMAX-ghostcells,IMAX-ghostcells);
          reldens_3D = permute(reldens_3D,[3 2 1]);
          flowbuoyancy = reldens_3D.*inplume;
          
%%% #TODO# change cmin/cmax/(numcolors?) to realistic once testing on
%%% cluster. Or put in pulse_post so they aren't hard coded.
        % Define buoyancy colormap: red = rise, blue = collapse.
          numcolors = 256;
          zeropoint = abs(min(reldensity))/(abs(max(reldensity))+abs(min(reldensity)));
          cmaplims = [1 0 0;    % red
                      1 1 1;    % white
                      0 0 1];   % blue
          cmin = -9;
          cmax = 1;
          fixcpts = [numcolors-1 numcolors*(1-(abs(cmax)/(abs(cmin)+abs(cmax)))) 0];
          cmap = interp1(fixcpts/numcolors,cmaplims,linspace(0,1,numcolors));
          colormap(cmap)
          
          view(saz,sel)
          hFB = slice(flowbuoyancy,sdistX*IMAX,sdistY*KMAX,sdistZ*JMAX);
          hFB.FaceColor = 'interp';
          hFB.EdgeColor = 'none';
          colorbar
          caxis([cmin cmax]);
          tL = pulsetitle(varname,PULSE,time,t,run,FREQ);
          title(tL,'FontSize',12,'FontWeight','bold');
          
%%% #TODO# clean up this section
          cd(dir)
            vidfig = 'buoyancy.jpg';
            saveas(fig,vidfig)
            img = imread(vidfig);
            writeVideo(vidFlowDens,img);
    
    end

    cd(dir)
    close(vidFlowDens);
    cd(postpath)
    sprintf('Flow density processing complete. \nvidFlowDens_%d has been saved to %s',run,dir)
    
end


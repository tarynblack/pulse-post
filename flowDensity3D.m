function [ vidFlowDens ] = flowDensity3D( run,dir,vis,IMAX,JMAX,KMAX,...
    ghostcells,postpath,RO_S1,RO_S2,RO_S3,plumeedge,PULSE,FREQ,time,...
    tickx,labelx,labelxunit,ticky,labely,labelyunit,tickz,labelz,...
    labelzunit,XRES,YRES,ZRES,sdistX,sdistY,sdistZ,flowDensity_cmin,...
    flowDensity_cmax,titlerun,flowBuoyancy_cmin,flowBuoyancy_cmax,...
    timesteps,imtype )
%flowDensity3D calculates the net density of the flow from gas and particle
%densities and volume fractions.
%   Detailed explanation goes here
%
%   Special functions called: varchunk3D, pulsetitle
%   Last edit: Taryn Black, 18 November 2015

    varD = 'Flow density';
    varB = 'Flow buoyancy';

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
    cd(dir)
    
    figDens = figure('Name','Plume Density','visible',vis,'units','normalized','outerposition',[0 0 0.5 1]);
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
    
    figBuoy = figure('Name','Plume Buoyancy','visible',vis,'units','normalized','outerposition',[0.5 0 0.5 1]);
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
    vidFlowDens = VideoWriter(sprintf('vidFlowDens_%s.avi',run));
    vidFlowDens.Quality = 100;
    vidFlowDens.FrameRate = 10;
    open(vidFlowDens);
    set(gcf,'Visible',vis);
    
    vidFlowBuoy = VideoWriter(sprintf('vidFlowBuoy_%s.avi',run));
    vidFlowBuoy.Quality = 100;
    vidFlowBuoy.FrameRate = 10;
    open(vidFlowBuoy);
    set(gcf,'Visible',vis);

%%% Calculate flow density at each timestep, plot, and save video.
    fID_EPG = fopen(sprintf('EP_G_%s',run));
    fID_ROG = fopen(sprintf('RO_G_%s',run));
    fID_RS1 = fopen(sprintf('ROP_S1_%s',run));
    fID_RS2 = fopen(sprintf('ROP_S2_%s',run));
    fID_RS3 = fopen(sprintf('ROP_S3_%s',run));
    
    t = 0;
    while t <= timesteps %~feof(fID_EPG)
        
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
        
      % Plot density as a slice through the domain.
        figure(figDens)
        view(saz,sel)
        hD = slice(domaindensity,sdistX*IMAX,sdistY*KMAX,sdistZ*JMAX);
        hD.FaceColor = 'interp';
        hD.EdgeColor = 'none';
        colormap(figDens,'default')
        colorbar
        caxis([flowDensity_cmin flowDensity_cmax]);
        tLD = pulsetitle(varD,PULSE,time,t,titlerun,FREQ);
        title(tLD,'FontSize',12,'FontWeight','bold');
        hold on
        
      % Determine whether flow is buoyant (less dense than atmosphere)
      % using mean flow/atmospheric densities at each altitude slice.
        inplume = EPG <= plumeedge;
        inatmos = EPG > plumeedge;
        
        flowdensity = domaindensity.*inplume;
        atmsdensity = domaindensity.*inatmos;
        
        for i = 1:(JMAX-ghostcells)
            z_atmsd = atmsdensity(:,:,i);
            avgatmsdens_altitude(i) = mean(z_atmsd(z_atmsd~=0));
        end
        
        % Relative density of atmosphere to flow. Positive = flow is
        % positively buoyant (less dense than atmosphere). Negative = flow
        % is negatively buoyant (more dense than atmosphere). Value =
        % difference in density
          avgatmsdens_3D = repmat(avgatmsdens_altitude',1,KMAX-ghostcells,IMAX-ghostcells);
          avgatmsdens_3D = permute(avgatmsdens_3D,[3 2 1]);
          avgatmsdens_3D_inplume = avgatmsdens_3D.*inplume;
          flowbuoyancy = avgatmsdens_3D_inplume - flowdensity;
            mindens = min(min(min(flowbuoyancy(flowbuoyancy~=0))));
            maxdens = max(max(max(flowbuoyancy(flowbuoyancy~=0))));
          
        % Define buoyancy colormap: red = rise, blue = collapse.
          numcolors = 256;
          zeropoint = abs(mindens)/(abs(maxdens)+abs(mindens));
          cmaplims = [1 0 0;    % red
                      1 1 1;    % white
                      0 0 1];   % blue
          fixcpts = [numcolors-1 numcolors*(1-(abs(flowBuoyancy_cmax)/(abs(flowBuoyancy_cmin)+abs(flowBuoyancy_cmax)))) 0];
          cmapB = interp1(fixcpts/numcolors,cmaplims,linspace(0,1,numcolors));
          colormap(figBuoy,cmapB)
          
          figure(figBuoy)
          view(saz,sel)
          hB = slice(flowbuoyancy,sdistX*IMAX,sdistY*KMAX,sdistZ*JMAX);
          hB.FaceColor = 'interp';
          hB.EdgeColor = 'none';
          colorbar
          caxis([flowBuoyancy_cmin flowBuoyancy_cmax]);
          tLB = pulsetitle(varB,PULSE,time,t,titlerun,FREQ);
          title(tLB,'FontSize',12,'FontWeight','bold');
          
          
          cd(dir)
            vidfigD = 'FlowDensity.jpg';
            saveas(figDens,vidfigD)
            imgD = imread(vidfigD);
            writeVideo(vidFlowDens,imgD);
            
          cd(dir)
            vidfigB = 'FlowBuoyancy.jpg';
            saveas(figBuoy,vidfigB)
            imgB = imread(vidfigB);
            writeVideo(vidFlowBuoy,imgB);
            
        % Save each timestep as an individual figure in either a
        % multipage tif file or other image filetype (user-specified).
          if strcmp(imtype,'tif') == 1 || strcmp(imtype,'tiff') == 1
              imwrite(imgD,sprintf('FlowDens_tsteps_%s.tif',run),'WriteMode','append')
          else saveas(figDens,sprintf('FlowDens_%03ds_%s.%s',time(t),run,imtype));
          end
          
        % Save each timestep as an individual figure in either a
        % multipage tif file or other image filetype (user-specified).
          if strcmp(imtype,'tif') == 1 || strcmp(imtype,'tiff') == 1
              imwrite(imgB,sprintf('FlowBuoy_tsteps_%s.tif',run),'WriteMode','append')
          else saveas(figBuoy,sprintf('FlowBuoy_%03ds_%s.%s',time(t),run,imtype));
          end
    
    end

    cd(dir)
    close(vidFlowDens);
    close(vidFlowBuoy);
    cd(postpath)
    sprintf('Flow density processing complete.\nvidFlowDens_%s has been saved to %s.\nvidFlowBuoy_%s has been saved to %s',run,dir,run,dir)
    
end


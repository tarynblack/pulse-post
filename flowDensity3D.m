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
%   Last edit: Taryn Black, 9 December 2015

  % Clear directory of appending files from previous processing attempts
    cd(dir)
    delete('flowdensity_*','atmsdensity_*','flowbuoyancy_*',...
        'FlowDens_*','FlowBuoy_*');
    
    
  % ----------------------- FIGURE INITIALIZATION ----------------------- %
  % Define variable names for figures
    varD = 'Flow density';
    varB = 'Flow buoyancy';
    cd(dir)

  % Ensure that 'no slice' directions are empty and determine figure
  % viewing angle based on slice direction
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
    
  % Flow density slice: figure and axes properties     
    figDens = figure('Name','Plume Density','visible',vis,'units',...
        'normalized','outerposition',[0 0 0.5 1]);
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
    
  % Flow buoyancy slice: figure and axes properties
    figBuoy = figure('Name','Plume Buoyancy','visible',vis,'units',...
        'normalized','outerposition',[0.5 0 0.5 1]);
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
    
  % Flow density slice: video
    vidFlowDens = VideoWriter(sprintf('vidFlowDens_%s.avi',run));
    vidFlowDens.Quality = 100;
    vidFlowDens.FrameRate = 10;
    open(vidFlowDens);
    set(gcf,'Visible',vis);
    
  % Flow buoyancy: video
    vidFlowBuoy = VideoWriter(sprintf('vidFlowBuoy_%s.avi',run));
    vidFlowBuoy.Quality = 100;
    vidFlowBuoy.FrameRate = 10;
    open(vidFlowBuoy);
    set(gcf,'Visible',vis);
  % ===================================================================== %

  
  % File import specifications: columns to read or skip for each variable
    EGimport   = '%f%*f%*f%*f%*f%*f%*f';
    EPS1import = '%*f%f%*f%*f%*f%*f%*f';
    EPS2import = '%*f%*f%f%*f%*f%*f%*f';
    EPS3import = '%*f%*f%*f%f%*f%*f%*f';
    ROGimport  = '%*f%f%*f%*f%*f%*f%*f';
  
    
  % =================== B E G I N   T I M E   L O O P =================== %
    t = 0;
    while t <= timesteps 
        
        t = t+1;
       
      % Queue up current timestep files
        cd(dir)
        fclose('all');
        clear fID*;
        fID_EPG = fopen(sprintf('EP_t%02d.txt',t));
        fID_EPS1 = fopen(sprintf('EP_t%02d.txt',t));
        fID_EPS2 = fopen(sprintf('EP_t%02d.txt',t));
        fID_EPS3 = fopen(sprintf('EP_t%02d.txt',t));
        fID_ROG = fopen(sprintf('Current_Density_t%02d.txt',t));
        cd(postpath)
        
      % Prepare vol fracs and gas dens for full domain at current timestep
        try
            EPG = varchunk3D(fID_EPG,EGimport,IMAX,JMAX,KMAX,ghostcells);
        catch ME
            warning('Error in varchunk3D at t=%d s:\n%s\nContinuing to next simulation.',...
                time(t),ME.identifier)
            break
        end
        ROG  = varchunk3D(fID_ROG,ROGimport,IMAX,JMAX,KMAX,ghostcells);
        EPS1 = varchunk3D(fID_EPS1,EPS1import,IMAX,JMAX,KMAX,ghostcells);
        EPS2 = varchunk3D(fID_EPS2,EPS2import,IMAX,JMAX,KMAX,ghostcells);
        EPS3 = varchunk3D(fID_EPS3,EPS3import,IMAX,JMAX,KMAX,ghostcells);
        
      % Skip processing for first timestep when there is no plume.
        if t==1;
            continue
        end        
        cla;
        
      % Calculate flow density at every point in domain
        domaindensity = (EPG.*ROG)+(EPS1*RO_S1)+(EPS2*RO_S2)+(EPS3*RO_S3);
        
        
      % ------------------- FLOW DENSITY SLICE FIGURE ------------------- %
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
      % ================================================================= %
      
      
      % --------------------- OVERLAY PLUME OUTLINE --------------------- %
        figure(figDens)
        hEPD = contourslice(EPG,sdistX*IMAX,sdistY*KMAX,0,[plumeedge plumeedge]);
        set(hEPD,'EdgeColor',[1 1 1],'LineWidth',0.5);
      % ================================================================= %
        
      
      % -------------------- CALCULATE FLOW BUOYANCY -------------------- %
      % Determine whether domain elements are inside or outside plume
        inplume = EPG <= plumeedge;
        inatmos = EPG > plumeedge;
        
      % Separate domain density matrix into flow and atmosphere densities
        flowdensity = domaindensity.*inplume;
        atmsdensity = domaindensity.*inatmos;
        
      % Calculate average atmospheric density at each altitude gridpoint
        for i = 1:(JMAX-ghostcells)
            z_atmsd = atmsdensity(:,:,i);
            avgatmsdens_altitude(i) = mean(z_atmsd(z_atmsd~=0));
        end
        
      % Calculate density of flow relative to atmosphere.
      % Positive: flow is positively buoyant (less dense than atmosphere)
      % Negative: flow is negatively buoyant (more dense than atmosphere)
        avgatmsdens_3D = repmat(avgatmsdens_altitude',1,KMAX-ghostcells,...
            IMAX-ghostcells);
        avgatmsdens_3D = permute(avgatmsdens_3D,[3 2 1]);
        avgatmsdens_3D_inplume = avgatmsdens_3D.*inplume;
        flowbuoyancy = avgatmsdens_3D_inplume - flowdensity;
        mindens = min(min(min(flowbuoyancy(flowbuoyancy~=0))));
        maxdens = max(max(max(flowbuoyancy(flowbuoyancy~=0))));
      % ================================================================= %
      
      
      % --------------------- BUOYANCY SLICE FIGURES -------------------- %
      % Define buoyancy colormap: red = rise, blue = collapse.
        numcolors = 256;
        zeropoint = abs(mindens)/(abs(maxdens)+abs(mindens));
        cmaplims = [1 0 0;    % red
                    1 1 1;    % white
                    0 0 1];   % blue
        fixcpts = [numcolors-1 numcolors*(1-(abs(flowBuoyancy_cmax)/...
            (abs(flowBuoyancy_cmin)+abs(flowBuoyancy_cmax)))) 0];
        cmapB = interp1(fixcpts/numcolors,cmaplims,linspace(0,1,numcolors));
        colormap(figBuoy,cmapB)
          
      % Plot slice of flow buoyancy ~ relative density
        figure(figBuoy)
        view(saz,sel)
        hB = slice(flowbuoyancy,sdistX*IMAX,sdistY*KMAX,sdistZ*JMAX);
          hB.FaceColor = 'interp';
          hB.EdgeColor = 'none';
        colorbar
          caxis([flowBuoyancy_cmin flowBuoyancy_cmax]);
        tLB = pulsetitle(varB,PULSE,time,t,titlerun,FREQ);
        title(tLB,'FontSize',12,'FontWeight','bold');
      % ================================================================= %
      
      
      % --------------------- OVERLAY PLUME OUTLINE --------------------- %
        figure(figBuoy)
        hEPB = contourslice(EPG,sdistX*IMAX,sdistY*KMAX,0,[plumeedge plumeedge]);
        set(hEPB,'EdgeColor',[0 0 0],'LineWidth',0.5);
      % ================================================================= %
        
        
      % --------- SAVE CURRENT FRAMES TO VIDEOS AND IMAGE FILES --------- %
        cd(dir)
        
      % Append current flow density frame to vidFlowDens
        vidfigD = 'FlowDensity.jpg';
        saveas(figDens,vidfigD)
        imgD = imread(vidfigD);
        writeVideo(vidFlowDens,imgD);
          
      % Append current flow buoyancy frame to vidFlowBuoy
        vidfigB = 'FlowBuoyancy.jpg';
        saveas(figBuoy,vidfigB)
        imgB = imread(vidfigB);
        writeVideo(vidFlowBuoy,imgB);
          
      % Save density and buoyancy calculations at each timestep
        dlmwrite(fullfile(sprintf('%s',dir),sprintf('flowdensity_t%03d.txt',...
            time(t))),flowdensity,'delimiter','\t','precision','%g');
        dlmwrite(fullfile(sprintf('%s',dir),sprintf('atmsdensity_t%03d.txt',...
            time(t))),atmsdensity,'delimiter','\t','precision','%g');
        dlmwrite(fullfile(sprintf('%s',dir),sprintf('flowbuoyancy_t%03d.txt',...
            time(t))),flowbuoyancy,'delimiter','\t','precision','%g');
          
      % If user-specified image filetype is tif, append current timestep
      % frame to multipage tif file. Otherwise, save frame as independent
      % image named by timestep.
        if strcmp(imtype,'tif') == 1 || strcmp(imtype,'tiff') == 1
            imwrite(imgD,sprintf('FlowDens_tsteps_%s.tif',run),'WriteMode','append')
            imwrite(imgB,sprintf('FlowBuoy_tsteps_%s.tif',run),'WriteMode','append')
        else
            saveas(figDens,sprintf('FlowDens_%03ds_%s.%s',time(t),run,imtype));
            saveas(figBuoy,sprintf('FlowBuoy_%03ds_%s.%s',time(t),run,imtype));
        end
      % ================================================================= %

    end
  % ===================== E N D   T I M E   L O O P ===================== %
  
  
  % End video write and finish video files
    cd(dir)
    close(vidFlowDens);
    close(vidFlowBuoy);
    
    cd(postpath)
    disp('Flow density processing complete.')
    fprintf('vidFlowDens_%s has been saved to %s.\n',run,dir)
    fprintf('vidFlowBuoy_%s has been saved to %s.\n',run,dir)
  
end


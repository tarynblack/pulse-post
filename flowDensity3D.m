function [ vidFlowDens ] = flowDensity3D( run,runpath,vis,IMAX,JMAX,KMAX,...
    ghostcells,postpath,RO_S1,RO_S2,RO_S3,plumeedge,PULSE,FREQ,time,...
    tickx,labelx,labelXunit,ticky,labely,labelYunit,tickz,labelz,...
    labelZunit,XRES,YRES,ZRES,sdistX,sdistY,sdistZ,flowDensity_cmin,...
    flowDensity_cmax,titlerun,flowBuoyancy_cmin,flowBuoyancy_cmax,...
    timesteps,imtype,savepath )
%flowDensity3D calculates the net density of the flow from gas and particle
%densities and volume fractions.
%   Detailed explanation goes here
%
%   Special functions called: varchunk3D, pulsetitle
%   Last edit: Taryn Black, 2 March 2016

  % Clear directory of appending files from previous processing attempts
    cd(savepath)
    delete('flowdensity_*','atmsdensity_*','flowreldens_*',...
        'FlowDens_*','FlowRelD_*');
    
    
  % ----------------------- FIGURE INITIALIZATION ----------------------- %
  % Define variable names for figures
    varD = 'Flow density';
    varB = 'Relative density';
    cd(runpath)

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
        'normalized','outerposition',[0 0 0.4 1],'PaperPositionMode',...
        'auto','color','w');
    axDens = axes('Parent',figDens,'box','on','TickDir','in','FontSize',12);
    hold on
    grid(axDens,'on');axDens.Layer = 'top';
    view(axDens,saz,sel)
    axis(axDens,'equal',[0,IMAX-ghostcells,0,KMAX-ghostcells,0,...
        JMAX-ghostcells]);
    set(axDens,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,...
        'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,...
        'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely)
    xlabel(axDens,sprintf('\\bf Distance_x (%s)',labelXunit))
    ylabel(axDens,sprintf('\\bf Distance_z (%s)',labelZunit))
    zlabel(axDens,sprintf('\\bf Altitude (%s)',labelYunit))
    colormap(figDens,'default')
    cbDens = colorbar(axDens,'AxisLocation','in','FontSize',12);
    cbDens.Label.String = '\bfFlow Density (kg/m^3)';
        
    
  % Flow relative density slice: figure and axes properties
    figRelD = figure('Name','Relative density','visible',vis,'units',...
        'normalized','outerposition',[0.5 0 0.4 1],'PaperPositionMode',...
        'auto','color','w');
    axRelD = axes('Parent',figRelD,'box','on','TickDir','in','FontSize',12);
    hold on
    grid(axRelD,'on');axRelD.Layer = 'top';
    view(axRelD,saz,sel)
    axis(axRelD,'equal',[0,IMAX-ghostcells,0,KMAX-ghostcells,0,...
        JMAX-ghostcells]);
    set(axRelD,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,...
        'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,...
        'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely)
    xlabel(axRelD,sprintf('\\bf Distance_x (%s)',labelXunit))
    ylabel(axRelD,sprintf('\\bf Distance_z (%s)',labelZunit))
    zlabel(axRelD,sprintf('\\bf Altitude (%s)',labelYunit))
    cbRelD = colorbar(axRelD,'AxisLocation','in','FontSize',12);
    cbRelD.Label.String = '\bfFlow Density relative to Atmosphere (kg/m^3)';
    
  % Flow density slice: video
    cd(savepath)
    vidFlowDens = VideoWriter(sprintf('vidFlowDens_%s.avi',run));
    vidFlowDens.Quality = 100;
    vidFlowDens.FrameRate = 10;
    open(vidFlowDens);
    set(gcf,'Visible',vis);
    
  % Flow relative density: video
    vidFlowReld = VideoWriter(sprintf('vidFlowRelD_%s.avi',run));
    vidFlowReld.Quality = 100;
    vidFlowReld.FrameRate = 10;
    open(vidFlowReld);
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
        cd(runpath)
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
        cla(axDens);
        hD = slice(0.5:(IMAX-ghostcells-0.5),0.5:(KMAX-ghostcells-0.5),...
            0.5:(JMAX-ghostcells-0.5),domaindensity,...
            sdistX*(IMAX-ghostcells),sdistY*(KMAX-ghostcells),...
            sdistZ*(JMAX-ghostcells));
        hD.FaceColor = 'interp';
        hD.EdgeColor = 'none';
        caxis(axDens,[flowDensity_cmin flowDensity_cmax]);
        tLD = pulsetitle(varD,PULSE,time,t,titlerun,FREQ);
        title(tLD,'FontSize',12,'FontWeight','bold');
      % ================================================================= %
      
      
      % --------------------- OVERLAY PLUME OUTLINE --------------------- %
        figure(figDens)
        hEPD = contourslice(EPG,sdistX*(IMAX-ghostcells),...
            sdistY*(KMAX-ghostcells),0,[plumeedge plumeedge]);
        set(hEPD,'EdgeColor',[1 1 1],'LineWidth',0.5);
      % ================================================================= %
        
      
      % ---------------- CALCULATE FLOW RELATIVE DENSITY ---------------- %
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
        flowreldens = avgatmsdens_3D_inplume - flowdensity;
%         mindens = min(flowreldens(flowreldens(:)~=0));
%         maxdens = max(flowreldens(flowreldens(:)~=0));
      % ================================================================= %
      
      
      % ----------------- RELATIVE DENSITY SLICE FIGURES ---------------- %
      % Define relative density colormap: red = rise, blue = collapse.
        numcolors = 256;
%         zeropoint = abs(mindens)/(abs(maxdens)+abs(mindens));
        cmaplims = [1 0 0;    % red
                    1 1 1;    % white
                    0 0 1];   % blue
        fixcpts = [numcolors-1 numcolors*(1-(abs(flowBuoyancy_cmax)/...
            (abs(flowBuoyancy_cmin)+abs(flowBuoyancy_cmax)))) 0];
        cmapB = interp1(fixcpts/numcolors,cmaplims,linspace(0,1,numcolors));
        colormap(figRelD,cmapB)
          
      % Plot slice of flow relative density
        figure(figRelD)
        cla(axRelD);
        hB = slice(0.5:(IMAX-ghostcells-0.5),0.5:(KMAX-ghostcells-0.5),...
            0.5:(JMAX-ghostcells-0.5),flowreldens,...
            sdistX*(IMAX-ghostcells),sdistY*(KMAX-ghostcells),...
            sdistZ*(JMAX-ghostcells));
        hB.FaceColor = 'interp';
        hB.EdgeColor = 'none';
        caxis(axRelD,[flowBuoyancy_cmin flowBuoyancy_cmax]);
        tLB = pulsetitle(varB,PULSE,time,t,titlerun,FREQ);
        title(tLB,'FontSize',12,'FontWeight','bold');
      % ================================================================= %
      
      
      % --------------------- OVERLAY PLUME OUTLINE --------------------- %
        figure(figRelD)
        hEPB = contourslice(EPG,sdistX*(IMAX-ghostcells),...
            sdistY*(KMAX-ghostcells),0,[plumeedge plumeedge]);
        set(hEPB,'EdgeColor',[0 0 0],'LineWidth',0.5);
      % ================================================================= %
        
        
      % --------- SAVE CURRENT FRAMES TO VIDEOS AND IMAGE FILES --------- %
        cd(savepath)
        
      % Append current flow density frame to vidFlowDens
        vidfigD = 'FlowDensity.jpg';
        saveas(figDens,fullfile(savepath,vidfigD))
        imgD = imread(vidfigD);
        writeVideo(vidFlowDens,imgD);
          
      % Append current flow relative density frame to vidFlowRelD
        vidfigRD = 'FlowRelDensity.jpg';
        saveas(figRelD,fullfile(savepath,vidfigRD))
        imgRD = imread(vidfigRD);
        writeVideo(vidFlowReld,imgRD);
          
      % Save density and relative density calculations at each timestep
        dlmwrite(fullfile(savepath,sprintf('flowdensity_t%03d.txt',...
            time(t))),flowdensity,'delimiter','\t','precision','%g');
        dlmwrite(fullfile(savepath,sprintf('atmsdensity_t%03d.txt',...
            time(t))),atmsdensity,'delimiter','\t','precision','%g');
        dlmwrite(fullfile(savepath,sprintf('flowreldens_t%03d.txt',...
            time(t))),flowreldens,'delimiter','\t','precision','%g');
          
      % If user-specified image filetype is tif, append current timestep
      % frame to multipage tif file. Otherwise, save frame as independent
      % image named by timestep.
        if strcmp(imtype,'tif') == 1 || strcmp(imtype,'tiff') == 1
            imwrite(imgD,fullfile(savepath,sprintf('FlowDens_tsteps_%s.tif',...
                run)),'tif','WriteMode','append')
            imwrite(imgRD,fullfile(savepath,sprintf('FlowRelD_tsteps_%s.tif',...
                run)),'tif','WriteMode','append')
        else
            saveas(figDens,fullfile(savepath,...
                sprintf('FlowDens_%03ds_%s.%s',time(t),run,imtype)));
            saveas(figRelD,fullfile(savepath,...
                sprintf('FlowRelD_%03ds_%s.%s',time(t),run,imtype)));
        end
      % ================================================================= %

    end
  % ===================== E N D   T I M E   L O O P ===================== %
  
  
  % End video write and finish video files
    cd(savepath)
    close(vidFlowDens);
    close(vidFlowReld);
    
    cd(postpath)
    disp('Flow density processing complete.')
    fprintf('vidFlowDens_%s has been saved to %s.\n',run,savepath)
    fprintf('vidFlowRelD_%s has been saved to %s.\n',run,savepath)
  
end

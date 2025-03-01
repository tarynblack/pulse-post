function [ vidFlowDens ] = flowDensity3D( run,runpath,vis,IMAX,JMAX,KMAX,...
    ghostcells,postpath,RO_S1,RO_S2,RO_S3,plumeedge,PULSE,FREQ,time,...
    tickx,labelx,labelXunit,ticky,labely,labelYunit,tickz,labelz,...
    labelZunit,XRES,YRES,ZRES,sdistX,sdistY,sdistZ,flowDensity_crange,...
    titlerun,flowBuoyancy_crange,timesteps,imtype,savepath,readEPG,...
    fnameEPG,readROG,fnameROG,readEPS1,fnameEPS1,readEPS2,fnameEPS2,...
    readEPS3,fnameEPS3,jetheight )
%flowDensity3D calculates the net density of the flow from gas and particle
%densities and volume fractions.
%   Detailed explanation goes here
%
%   Special functions called: loadTimestep3D, pulsetitle
%   Last edit: Taryn Black, 14 April 2016

  % Clear directory of appending files from previous processing attempts
    cd(savepath)
    delete('flowdensity_*','atmsdensity_*','flowreldens_*',...
        'FlowDens*','FlowRelD*','avgRelDens_*');
    
    
  % ----------------------- FIGURE INITIALIZATION ----------------------- %
  % Define variable names for figures
    varD = 'Bulk density';
    varB = 'Density contrast';
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
        sel = 90;
    elseif isempty(sdistY) && isempty(sdistZ)
        saz = 90;
        sel = 0;
    elseif isempty(sdistX) && isempty(sdistZ)
        saz = 0;
        sel = 0;
    else [saz,sel] = view(3);
    end
    
  % Bulk density slice: figure and axes properties     
    figDens = figure('Name','Bulk density','visible',vis,'units',...
        'centimeters','outerposition',[0 0 19 18.75],'PaperPositionMode',...
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
    cbDens.Label.String = '\bfBulk density (kg/m^3)';
        
    
  % Density contrast slice: figure and axes properties
    figRelD = figure('Name','Density contrast','visible',vis,'units',...
        'centimeters','outerposition',[16.66 0 19 18.75],'PaperPositionMode',...
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
    cbRelD.Label.String = '\bf\rho_{atmosphere} - \rho_{mixture} (kg/m^3)';
  % Define relative density colormap: red = rise, blue = collapse.
    numcolors = 256;
    cmaplims = [1 0 0;    % red
                1 1 1;    % white
                0 0 1];   % blue
    fixcpts = [numcolors-1 numcolors*(1-(abs(flowBuoyancy_crange(2))/...
        (abs(flowBuoyancy_crange(1))+abs(flowBuoyancy_crange(2))))) 0];
    cmapB = interp1(fixcpts/numcolors,cmaplims,linspace(0,1,numcolors));
    colormap(figRelD,cmapB)
    
  % Bulk density slice: video
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
   
    
  % Preallocate vectors
    avgRDJH = zeros(1,timesteps);
  
    
  % =================== B E G I N   T I M E   L O O P =================== %
    t = 0;
    while t <= timesteps 
        
        t = t+1;
       
      % Queue up current timestep files
        cd(runpath)
        fclose('all');
        clear fID*;
        
        cd(postpath)
        fID_EPG = fileReadType(fnameEPG,readEPG,t,runpath,postpath);
        fID_EPS1 = fileReadType(fnameEPS1,readEPS1,t,runpath,postpath);
        fID_EPS2 = fileReadType(fnameEPS2,readEPS2,t,runpath,postpath);
        fID_EPS3 = fileReadType(fnameEPS3,readEPS3,t,runpath,postpath);
        fID_ROG = fileReadType(fnameROG,readROG,t,runpath,postpath);
        
      % Prepare vol fracs and gas dens for full domain at current timestep
        try
            EPG = loadTimestep3D(fID_EPG,EGimport,readEPG,IMAX,JMAX,KMAX,ghostcells);
        catch ME
            warning('Error in loadTimestep3D at t=%d s:\n%s\nContinuing to next simulation.',...
                time(t),ME.identifier)
            break
        end
        ROG  = loadTimestep3D(fID_ROG,ROGimport,readROG,IMAX,JMAX,KMAX,ghostcells);
        EPS1 = loadTimestep3D(fID_EPS1,EPS1import,readEPS1,IMAX,JMAX,KMAX,ghostcells);
        EPS2 = loadTimestep3D(fID_EPS2,EPS2import,readEPS2,IMAX,JMAX,KMAX,ghostcells);
        EPS3 = loadTimestep3D(fID_EPS3,EPS3import,readEPS3,IMAX,JMAX,KMAX,ghostcells);
        
      % Skip processing for first timestep when there is no plume.
        if t==1;
            continue
        end        
        cla;
        
      % Calculate bulk density at every point in domain
        domaindensity = (EPG.*ROG)+(EPS1*RO_S1)+(EPS2*RO_S2)+(EPS3*RO_S3);
        
        
      % ------------------- BULK DENSITY SLICE FIGURE ------------------- %
        figure(figDens)
        cla(axDens);
        hD = slice(0.5:(IMAX-ghostcells-0.5),0.5:(KMAX-ghostcells-0.5),...
            0.5:(JMAX-ghostcells-0.5),domaindensity,...
            sdistX*(IMAX-ghostcells),sdistY*(KMAX-ghostcells),...
            sdistZ*(JMAX-ghostcells));
        hD.FaceColor = 'interp';
        hD.EdgeColor = 'none';
        caxis(axDens,flowDensity_crange);
        tLD = pulsetitle(varD,PULSE,time,t,titlerun,FREQ);
        title(tLD,'FontSize',12,'FontWeight','bold');
        hEPD = contourslice(EPG,sdistX*(IMAX-ghostcells),...
            sdistY*(KMAX-ghostcells),0,[plumeedge plumeedge]);
        set(hEPD,'EdgeColor',[1 1 1],'LineWidth',0.5);
        set(figDens,'Visible',vis);
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

      % Calculate average relative density of flow at jet height
        inplumeJH  = inplume(:,:,round(jetheight));
        reldensJH  = flowreldens(:,:,round(jetheight));
        avgRDJH(t) = mean(reldensJH(inplumeJH));
        dlmwrite(fullfile(savepath,sprintf('avgRelDens_JetHeight_%s.txt',...
            run)),[time(t) avgRDJH(t)],'-append','delimiter','\t');
      % ================================================================= %
      
      
      % ----------------- RELATIVE DENSITY SLICE FIGURES ---------------- %
        figure(figRelD)
        cla(axRelD);
        hB = slice(0.5:(IMAX-ghostcells-0.5),0.5:(KMAX-ghostcells-0.5),...
            0.5:(JMAX-ghostcells-0.5),flowreldens,...
            sdistX*(IMAX-ghostcells),sdistY*(KMAX-ghostcells),...
            sdistZ*(JMAX-ghostcells));
        hB.FaceColor = 'interp';
        hB.EdgeColor = 'none';
        caxis(axRelD,flowBuoyancy_crange);
        tLB = pulsetitle(varB,PULSE,time,t,titlerun,FREQ);
        title(tLB,'FontSize',12,'FontWeight','bold');
        hEPB = contourslice(EPG,sdistX*(IMAX-ghostcells),...
            sdistY*(KMAX-ghostcells),0,[plumeedge plumeedge]);
        set(hEPB,'EdgeColor',[0 0 0],'LineWidth',0.5);
        set(figRelD,'Visible',vis);
      % ================================================================= %
        
        
      % --------- SAVE CURRENT FRAMES TO VIDEOS AND IMAGE FILES --------- %
        cd(savepath)
        
      % Append current bulk density frame to vidFlowDens
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
    
    
    if strcmp(PULSE,'T') == 1
      str = sprintf('%s: Unsteady flow %g Hz',titlerun,FREQ);
    elseif strcmp(PULSE,'F') == 1
      str = sprintf('%s: Steady flow',titlerun);
    end

  % --------- PLOT TIME SERIES OF RELATIVE DENSITY AT JET HEIGHT -------- %
    figRDJH = figure('Name','Jet height relative density','visible',vis,...
        'units','centimeters','outerposition',[0 0 33.33 15],...
        'PaperPositionMode','auto','color','w');
    axRDJH = axes('Parent',figRDJH,'box','on','TickDir','in','Fontsize',12);
    grid(axRDJH,'on');
    hold on
    allRDJH = load(sprintf('avgRelDens_JetHeight_%s.txt',run));
    allRDJH = allRDJH(:,2);
    negidx = allRDJH<=0;
    negRDJH = NaN(length(allRDJH),1);
    negRDJH(negidx) = allRDJH(negidx);
    hRDJH = plot(time(2:end),allRDJH,'r',time(2:end),negRDJH,'b','LineWidth',2);
    xlim([0 time(end)]);
    ylim(flowBuoyancy_crange)
    xlabel('\bfTime (s)')
    ylabel('\bf\rho_{atmosphere} - \rho_{mixture} (kg/m^3)')
    title(axRDJH,{sprintf('Mixture density contrast at jet height (%.3f km)',...
        jetheight*YRES/1000);sprintf('%s',str)});
    saveas(figRDJH,fullfile(savepath,sprintf('JetHeightBuoyancy_%s.jpg',run)));
  % ===================================================================== %

  
    cd(postpath)
    disp('Bulk density processing complete.')
    fprintf('vidFlowDens_%s has been saved to %s.\n',run,savepath)
    fprintf('vidFlowRelD_%s has been saved to %s.\n',run,savepath)
  
end

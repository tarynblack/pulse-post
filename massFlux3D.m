function [ vidMFlux ] = massFlux3D( runpath,vis,viewaz,viewel,ghostcells,...
    IMAX,JMAX,KMAX,tickx,ticky,tickz,XRES,YRES,ZRES,labelx,labely,labelz,...
    labelXunit,labelYunit,labelZunit,run,timesteps,postpath,massflux_alts,...
    RO_S1,RO_S2,RO_S3,plumeedge,massflux_crange,PULSE,FREQ,time,...
    titlerun,massflux_legend,imtype,savepath,readEPG,fnameEPG,...
    readROG,fnameROG,readVG,fnameVG,readVS1,fnameVS1,readVS2,fnameVS2,...
    readVS3,fnameVS3,readEPS1,fnameEPS1,readEPS2,fnameEPS2,...
    readEPS3,fnameEPS3,jetheight,MASSFLUX_SOL,VENT_R,MING,sdistX,sdistY,...
    sdistZ )
%massFlux3D Summary of this function goes here
%   Detailed explanation goes here
%
%   Functions called: loadTimestep3D; pulsetitle
%   Last edit: Taryn Black, 6 April 2016

  % Clear directory of appending files from previous processing attempts
    cd(savepath)
    delete('*massflux*','*MFlux*','*AvgNetMF*','Collapse*','*Ongaro*',...
        'NetMassFlux_*','EP*vent')
    
    
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
    
  % Subtightplot properties
    gap = [0.01 0.01];
    ht = 0.15;
    wd = 0.15;


  % ----------------------- FIGURE INITIALIZATION ----------------------- %
  % Define variable names for figures
    varMF = 'Solid mass flux';
    varMFZ = 'Horizontally averaged solid mass flux';
    varVFE = 'Volume fraction of gas and solids';
    cd(runpath)
    
  % Mass flux slice: figure and axes properties
    figMFlux = figure('Name','Mass flux','units','centimeters',...
        'outerposition',[0 0 28 18.75],'visible',vis,'PaperPositionMode',...
        'auto','color','w');
    cd(postpath)
    axMFlux1 = subtightplot(1,2,1,gap,ht,wd);
        hold on
        axis(axMFlux1,'equal',[0,IMAX-ghostcells,0,KMAX-ghostcells,0,...
            JMAX-ghostcells]);
        set(axMFlux1,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,...
            'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,...
            'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely);
        xlabel(axMFlux1,sprintf('\\bfDistance_x (%s)',labelXunit))
        ylabel(axMFlux1,sprintf('\\bfDistance_z (%s)',labelZunit))
        zlabel(axMFlux1,sprintf('\\bfAltitude (%s)',labelYunit))
    axMFlux2 = subtightplot(1,2,2,gap,ht,wd);
        hold on
        axis(axMFlux2,'equal',[0,IMAX-ghostcells,0,KMAX-ghostcells,0,...
            JMAX-ghostcells]);
        set(axMFlux2,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,...
          'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,...
          'ZTick',ticky(2:end)/YRES,'ZTickLabel','')
        xlabel(axMFlux2,sprintf('\\bf Distance_x (%s)',labelXunit))
        ylabel(axMFlux2,sprintf('\\bf Distance_z (%s)',labelZunit))
%         zlabel(axMFlux2,sprintf('\\bfAltitude (%s)',labelYunit))
    set([axMFlux1 axMFlux2],'box','on','TickDir','in','FontSize',12)
    grid(axMFlux1,'on');grid(axMFlux2,'on');
    axMFlux1.Layer = 'top';axMFlux2.Layer = 'top';
    view(axMFlux1,saz,sel);view(axMFlux2,viewaz,viewel);
    cmapMFlux = [winter;[0.9 0.9 0.9];flipud(autumn)];
    colormap(axMFlux1,cmapMFlux);colormap(axMFlux2,cmapMFlux);
    cbMFlux = colorbar(axMFlux2,'AxisLocation','in','FontSize',12);
    cbMFlux.Label.String = '\bfMass Flux (kg/m^2s)';
    cbMFlux.Ticks = -log10(abs(massflux_crange(1))):log10(massflux_crange(2));
    cbMFlux.TickLabels = cellstr(num2str([-10.^(abs(-log10(abs...
        (massflux_crange(1))):-1)) 0 10.^(1:log10(massflux_crange(2)))]','%.0e\n'));
    
  % Average mass flux time series: figure and axes properties
    figAvgMFZ = figure('Name','Spatially averaged mass flux with altitude',...
        'units','centimeters','outerposition',[0 0 33.33 18.75],'visible',...
        vis,'PaperPositionMode','auto','color','w');
    axAvgMFZ = axes('Parent',figAvgMFZ,'box','on','TickDir','in','FontSize',12);
    hold on
    grid(axAvgMFZ,'on');
    axis(axAvgMFZ,[2*massflux_crange(1),2*massflux_crange(2),0,JMAX-ghostcells]);
    set(axAvgMFZ,'YTick',ticky(2:end)/YRES,'YTickLabel',labely);
    xlabel(axAvgMFZ,'\bfNet mass flux (kg/m^2s)')
    ylabel(axAvgMFZ,sprintf('\\bfAltitude (%s)',labelYunit))
    hMFZ = plot(0,0,'DisplayName','Previous profiles');
    hBlk = plot(0,0,'k','LineWidth',2,'DisplayName',...
        'Current profile');
    hJet = plot(2*massflux_crange(1):1E3:2*massflux_crange(2),...
        jetheight*ones(1,length(2*massflux_crange(1):1E3:2*massflux_crange(2))),...
        '--','Color',[0.2 0.5 0.2],'LineWidth',1.5,'DisplayName',...
        sprintf('Jet height (%.3f km)',jetheight*YRES/1000));
    hMFZleg = legend(axAvgMFZ,[hBlk hMFZ hJet]);
    set(hMFZleg,'FontSize',12,'Location','Northwest')
    
        
  % Initialize video
    cd(savepath)
    vidMFlux = VideoWriter(sprintf('vidMFlux_%s.avi',run));
    vidMFlux.Quality = 100;
    vidMFlux.FrameRate = 10;
    open(vidMFlux);
    set(gcf,'Visible',vis);
  % ===================================================================== %
  
  
  % File import specifications: columns to read (%) or skip (%*) for each variable
    EPGimport = '%f%*f%*f%*f%*f%*f%*f';
%     ROGimport = '%*f%f%*f%*f%*f%*f%*f';
%     VGimport  = '%*f%f%*f%*f%*f%*f';
    EPS1import = '%*f%f%*f%*f%*f%*f%*f';
    EPS2import = '%*f%*f%f%*f%*f%*f%*f';
    EPS3import = '%*f%*f%*f%f%*f%*f%*f';
    VS1import = '%f%*f%*f';
    VS2import = '%*f%f%*f';
    VS3import = '%*f%*f%f';
    
  % Preallocate vectors
    netMF_alts = zeros(length(massflux_alts),timesteps);
    collapse_Ongaro = zeros(1,timesteps);
    
  
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
%         fID_ROG = fileReadType(fnameROG,readROG,t,runpath,postpath);
%         fID_VG = fileReadType(fnameVG,readVG,t,runpath,postpath);
        fID_EPS1 = fileReadType(fnameEPS1,readEPS1,t,runpath,postpath);
        fID_EPS2 = fileReadType(fnameEPS2,readEPS2,t,runpath,postpath);
        fID_EPS3 = fileReadType(fnameEPS3,readEPS3,t,runpath,postpath);
        fID_VS1 = fileReadType(fnameVS1,readVS1,t,runpath,postpath);
        fID_VS2 = fileReadType(fnameVS2,readVS2,t,runpath,postpath);
        fID_VS3 = fileReadType(fnameVS3,readVS3,t,runpath,postpath);
        
      % Prepare velocities for full domain at current timestep
        try
            EPG = loadTimestep3D(fID_EPG,EPGimport,readEPG,IMAX,JMAX,KMAX,ghostcells);
        catch ME
            warning('Error in loadTimestep3D at t=%d s:\n%s\nContinuing to next simulation.',...
                time(t),ME.identifier)
            break
        end
%         ROG = loadTimestep3D(fID_ROG,ROGimport,readROG,IMAX,JMAX,KMAX,ghostcells);
%         V_G  = loadTimestep3D(fID_VG,VGimport,readVG,IMAX,JMAX,KMAX,ghostcells);
        EPS1 = loadTimestep3D(fID_EPS1,EPS1import,readEPS1,IMAX,JMAX,KMAX,ghostcells);
        EPS2 = loadTimestep3D(fID_EPS2,EPS2import,readEPS2,IMAX,JMAX,KMAX,ghostcells);
        EPS3 = loadTimestep3D(fID_EPS3,EPS3import,readEPS3,IMAX,JMAX,KMAX,ghostcells);
        V_S1 = loadTimestep3D(fID_VS1,VS1import,readVS1,IMAX,JMAX,KMAX,ghostcells);
        V_S2 = loadTimestep3D(fID_VS2,VS2import,readVS2,IMAX,JMAX,KMAX,ghostcells);
        V_S3 = loadTimestep3D(fID_VS3,VS3import,readVS3,IMAX,JMAX,KMAX,ghostcells);
        
      % Skip processing for first timestep when there is no plume
        if t==1;
            continue
        end
        
      % Calculate vertical solid mass flux at specified altitudes
%         massflux = ROG.*V_G + RO_S1*V_S1 + RO_S2*V_S2 + RO_S3*V_S3;
        massflux = RO_S1*EPS1.*V_S1 + RO_S2*EPS2.*V_S2 + RO_S3*EPS3.*V_S3;
        netmassflux = squeeze(sum(sum(massflux)));
        netMF_alts(:,t) = netmassflux(massflux_alts);
    
      % Save calculated mass fluxes at each timestep
        dlmwrite(fullfile(savepath,sprintf('massflux_all_t%03d.txt',...
            time(t))),massflux,'delimiter','\t');
        dlmwrite(fullfile(savepath,sprintf('netmassflux_%s.txt',...
            run)),[time(t) netmassflux'],'-append','delimiter','\t');
        
      % Separate positive and negative mass fluxes for logarithmic plotting
        logMF = massflux;
        logMF(logMF>0) = log10(logMF(logMF>0));
        logMF(logMF<0) = -log10(abs(logMF(logMF<0)));
        
      % Calculate most-negative flux and vent flux for Ongaro criterion
%         ventflux = zeros(IMAX-ghostcells,KMAX-ghostcells);
        negMF  = abs(min(netmassflux(netmassflux<0)));
        if isempty(negMF)
            negMF = 0;
        end
%         for i = ((IMAX-ghostcells)/2)-(VENT_R/XRES):((IMAX-ghostcells)/2)+(VENT_R/XRES)
%             for k = ((KMAX-ghostcells)/2)-(VENT_R/ZRES):((KMAX-ghostcells)/2)+(VENT_R/ZRES)
%                 if (i - ((IMAX-ghostcells)/2))^2 + (k - ((KMAX-ghostcells)/2))^2 <= (VENT_R/XRES)^2
%                     ventflux(i,k) = massflux(i,k,1);
%                 end
%             end
%         end
%         netventflux = squeeze(sum(sum(ventflux)));
        collapse_Ongaro(t) = negMF/MASSFLUX_SOL;%netventflux;
        dlmwrite(fullfile(savepath,sprintf('collapseOngaro_%s.txt',run)),...
            [time(t) negMF MASSFLUX_SOL collapse_Ongaro(t)],...
            '-append','delimiter','\t');

        
        
      % --------------------- MASS FLUX SLICE FIGURE -------------------- %
        figure(figMFlux)
        cla(axMFlux1);cla(axMFlux2);
        hMF1 = slice(axMFlux1,0.5:(IMAX-ghostcells-0.5),0.5:(KMAX-ghostcells-0.5),...
            0.5:(JMAX-ghostcells-0.5),logMF,sdistX*(IMAX-ghostcells),...
            sdistY*(KMAX-ghostcells),sdistZ*(JMAX-ghostcells));
          hMF1.FaceColor = 'interp';
          hMF1.EdgeColor = 'none';
          hEPZ1 = contourslice(axMFlux1,EPG,sdistX*(IMAX-ghostcells),...
              sdistY*(KMAX-ghostcells),0,[plumeedge plumeedge]);
          set(hEPZ1,'EdgeColor',[1 1 1],'LineWidth',0.5);
          tMF1 = pulsetitle(varMF,PULSE,time,t,titlerun,FREQ);
          title(axMFlux1,tMF1,'FontWeight','bold');
          caxis(axMFlux1,[-log10(abs(massflux_crange(1))) log10(massflux_crange(2))]);
        hMF2 = slice(axMFlux2,1:IMAX-ghostcells,1:KMAX-ghostcells,...
            1:JMAX-ghostcells,logMF,[],[],massflux_alts);
          set(hMF2,'FaceColor','interp','EdgeColor','none')
          hEPZ2 = contourslice(axMFlux2,EPG,0,0,massflux_alts,...
              [plumeedge plumeedge]);
          set(hEPZ2,'EdgeColor',[1 1 1],'LineWidth',0.5);
          caxis(axMFlux2,[-log10(abs(massflux_crange(1))) log10(massflux_crange(2))]);
%         tMF = pulsetitle(varMF,PULSE,time,t,titlerun,FREQ);
        tMF2 = sprintf('Jet height: %.3f km',jetheight*YRES/1000);
        title(axMFlux2,[tMF2],'FontWeight','bold');
        PosMF1 = get(axMFlux1,'position');
        PosMF2 = get(axMFlux2,'position');
        PosMF2(3:4) = PosMF1(3:4);
        set(axMFlux2,'position',PosMF2);
      % ================================================================= %
      
      
      % ----------- SPATIALLY AVERAGED MASS FLUX WITH HEIGHT ------------ %
        figure(figAvgMFZ)
        hold on
        set(hMFZ,'Color',[0.55 0.55 0.55],'LineWidth',0.5);
        hMFZ = plot(netmassflux,1:JMAX-ghostcells,'k','LineWidth',2);        
        tMFZ = pulsetitle(varMFZ,PULSE,time,t,titlerun,FREQ);
        title(axAvgMFZ,tMFZ,'FontWeight','bold');
      % ================================================================= %
      
      
      % --------- SAVE CURRENT FRAMES TO VIDEOS AND IMAGE FILES --------- %
        cd(savepath)
        
      % Append current mass flux frame to vidMFlux
        vidfigMF = 'MFluxCurrent.jpg';
        saveas(figMFlux,fullfile(savepath,vidfigMF));
        imgMF = imread(vidfigMF);
        writeVideo(vidMFlux,imgMF);
        
      % Save current frame of average net mass flux figure
        figMFZ = 'AvgNetMFCurrent.jpg';
        saveas(figAvgMFZ,fullfile(savepath,figMFZ));
        imgMFZ = imread(figMFZ);
        
      % If user-specified image filetype is tif, append current timestep
      % frame to multipage tif file. Otherwise, save frame as independent
      % image named by timestep.
        if strcmp(imtype,'tif') == 1 || strcmp(imtype,'tiff') == 1
            imwrite(imgMF,fullfile(savepath,sprintf('MFlux_tsteps_%s.tif',...
                run)),'tif','WriteMode','append');
            imwrite(imgMFZ,fullfile(savepath,sprintf('AvgNetMF_tsteps_%s.tif',...
                run)),'tif','WriteMode','append');
        else
            saveas(figMF,fullfile(savepath,sprintf('MFlux_t%03d_%s.%s',...
                time(t),run,imtype)));
            saveas(figMFZ,fullfile(savepath,sprintf('AvgNetMF_t%03d_%s.%s',...
                time(t),run,imtype)));
        end
      % ================================================================= %
      
    end
  % ===================== E N D   T I M E   L O O P ===================== %
  
  
  % End video write and finish video files
    cd(savepath)
    close(vidMFlux);
    
    
    if strcmp(PULSE,'T') == 1
      str = sprintf('%s: Unsteady flow %g Hz',titlerun,FREQ);
    elseif strcmp(PULSE,'F') == 1
      str = sprintf('%s: Steady flow',titlerun);
    end
    
   % ------------------ NET MASS FLUX TIME SERIES PLOTS ------------------ % 
    figNetMF = figure('Name','Net Mass Flux','visible',vis,'units',...
        'centimeters','outerposition',[0 0 33.33 18.75],'PaperPositionMode',...
        'auto','color','w');
    axNetMF = axes('Parent',figNetMF,'box','on','FontSize',12);
    grid(axNetMF,'on');
    axis(axNetMF,[0,time(end),2*massflux_crange(1),2*massflux_crange(2)]);
    hold on
    plot(time,netMF_alts);
    legend(axNetMF,massflux_legend)
    title(axNetMF,sprintf('%s: Net solid mass flux through specified altitudes',str))
    xlabel(axNetMF,'\bfTime (s)')
    ylabel(axNetMF,'\bfNet mass flux (kg/m^2s)')
    saveas(figNetMF,fullfile(savepath,sprintf('NetMassFlux_tseries_%s.jpg',run)))
  % ===================================================================== %
  
  
  % ---------------- TIME-AVERAGED MASS FLUX WITH HEIGHT ---------------- %
    figure(figAvgMFZ)
    hold on
    set(hMFZ,'Color',[0.55 0.55 0.55],'LineWidth',0.5);
    allNMF = load(sprintf('netmassflux_%s.txt',run));
    avgNMF = mean(allNMF(:,2:end),1);
    hNMF = plot(avgNMF,1:JMAX-ghostcells,'-.','Color',[0 0.4 0.7],'LineWidth',3);
    set(hNMF,'DisplayName','Time-averaged profile')
    set(hMFZ,'DisplayName','Individual timestep profiles')
    title(axAvgMFZ,sprintf('%s: Time-averaged mass flux with height',str))
    hMFZleg = legend(axAvgMFZ,[hNMF hMFZ hJet]);
    set(hMFZleg,'FontSize',12,'Location','Northwest')
    saveas(figAvgMFZ,fullfile(savepath,sprintf('TimeAvgNetMF_%s.jpg',run)));
    
    avgNegMF = abs(min(avgNMF(avgNMF<0)));
    if isempty(avgNegMF)
        avgNegMF = 0;
    end
    avg_Ongaro = avgNegMF/MASSFLUX_SOL;
    dlmwrite(fullfile(savepath,sprintf('avgOngaroCrit_%s.txt',run)),...
        [avgNegMF MASSFLUX_SOL avg_Ongaro],'delimiter','\t');
  % ===================================================================== %
  
  
  % ---------------- COLLAPSE CRITERION TIME SERIES PLOT ---------------- %
    figCollapse = figure('Name','Collapse criterion','units','centimeters',...
        'outerposition',[0 0 33.33 18.75],'visible',vis,'PaperPositionMode',...
        'auto','color','w');
    axCollapse = axes('Parent',figCollapse,'box','on','TickDir','in',...
        'FontSize',12);
    grid(axCollapse,'on');
    axis(axCollapse,[0,time(end),0,1]);
    hold on
    plot(time,collapse_Ongaro,'.-',time,0.9*ones(1,length(time)),'k--',time,...
        0.65*ones(1,length(time)),'k-.',time,0.5*ones(1,length(time)),'k:');
    xlabel(axCollapse,'\bfTime (s)');
    ylabel(axCollapse,'\bfCollapse criterion ratio');
    title(axCollapse,sprintf('%s: Collapse criterion ratio...',str));
    legend('SUPERRATIO!','Near-total collapse','Partial collapse','Incipient collapse');
    saveas(figCollapse,fullfile(savepath,sprintf('CollapseCriterion_%s.jpg',run)));
  % ===================================================================== %
  
  
  cd(postpath)
  disp('Mass flux processing complete.')
  fprintf('vidMFlux_%s has been saved to %s.\n',run,savepath)

end

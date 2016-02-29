function [ vidMFlux ] = massFlux3D( dir,vis,viewaz,viewel,ghostcells,...
    IMAX,JMAX,KMAX,tickx,ticky,tickz,XRES,YRES,ZRES,labelx,labely,labelz,...
    labelXunit,labelYunit,labelZunit,run,timesteps,postpath,massflux_alts,...
    RO_S1,RO_S2,RO_S3,plumeedge,massflux_cmin,massflux_cmax,PULSE,FREQ,...
    time,titlerun,massflux_legend,imtype )
%massFlux3D Summary of this function goes here
%   Detailed explanation goes here
%
%   Functions called: varchunk3D; pulsetitle
%   Last edit: Taryn Black, 18 February 2016

  % Clear directory of appending files from previous processing attempts
    cd(dir)
    delete('massflux*','netmassflux*','*MFlux*')


  % ----------------------- FIGURE INITIALIZATION ----------------------- %
  % Define variable names for figures
    varMF = 'Mass flux';
    
  % Figure and axes properties
    figMFlux = figure('Name','Mass flux','units','normalized',...
        'outerposition',[0 0 0.4 1],'visible',vis,'PaperPositionMode',...
        'auto','color','w');
    axMFlux = axes('Parent',figMFlux,'box','on','TickDir','in','FontSize',12);
    hold on
    grid(axMFlux,'on');axMFlux.Layer = 'top';
    view(axMFlux,viewaz,viewel)
    axis(axMFlux,'equal',[0,IMAX-ghostcells,0,KMAX-ghostcells,0,...
        JMAX-ghostcells]);
    set(axMFlux,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,...
        'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,...
        'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely);
    xlabel(axMFlux,sprintf('\\bfDistance_x (%s)',labelXunit))
    ylabel(axMFlux,sprintf('\\bfDistance_z (%s)',labelZunit))
    zlabel(axMFlux,sprintf('\\bfAltitude (%s)',labelYunit))
    cbMFlux = colorbar(axMFlux,'AxisLocation','out','FontSize',12);
    cbMFlux.Label.String = '\bfMass Flux (kg/m^2s)';
    
  % Initialize video
    vidMFlux = VideoWriter(sprintf('vidMFlux_%s.avi',run));
    vidMFlux.Quality = 100;
    vidMFlux.FrameRate = 10;
    open(vidMFlux);
    set(gcf,'Visible',vis);
  % ===================================================================== %
  
  
  % File import specifications: columns to read (%) or skip (%*) for each variable
    EPGimport = '%f%*f%*f%*f%*f%*f%*f';
    ROGimport = '%*f%f%*f%*f%*f%*f%*f';
    VGimport  = '%*f%f%*f%*f%*f%*f';
%     VS1import = '%f%*f%*f';
%     VS2import = '%*f%f%*f';
%     VS3import = '%*f%*f%f';
    
    VS1import = '%*f%f%*f%*f%*f%*f%*f';
    VS2import = '%*f%*f%f%*f%*f%*f%*f';
    VS3import = '%*f%*f%*f%f%*f%*f%*f';
    
  
  % =================== B E G I N   T I M E   L O O P =================== %
    t = 0;
    while t <= timesteps
        
        t = t+1;
        
      % Queue up current timestep files
        cd(dir)
        fclose('all');
        clear fID*;
        fID_EPG = fopen(sprintf('EP_t%02d.txt',t));
        fID_ROG = fopen(sprintf('EP_t%02d.txt',t));
        fID_VG  = fopen(sprintf('U_G_t%02d.txt',t));
%         fID_VS1 = fopen(sprintf('V_S_t%02d.txt',t));
%         fID_VS2 = fopen(sprintf('V_S_t%02d.txt',t));
%         fID_VS3 = fopen(sprintf('V_S_t%02d.txt',t));

        fID_VS1 = fopen(sprintf('EP_t%02d.txt',t));
        fID_VS2 = fopen(sprintf('EP_t%02d.txt',t));
        fID_VS3 = fopen(sprintf('EP_t%02d.txt',t));
        cd(postpath)
        
      % Prepare velocities for full domain at current timestep
        try
            EPG = varchunk3D(fID_EPG,EPGimport,IMAX,JMAX,KMAX,ghostcells);
        catch ME
            warning('Error in varchunk3D at t=%d s:\n%s\nContinuing to next simulation.',...
                time(t),ME.identifier)
            break
        end
        ROG = varchunk3D(fID_ROG,ROGimport,IMAX,JMAX,KMAX,ghostcells);
        V_G  = varchunk3D(fID_VG,VGimport,IMAX,JMAX,KMAX,ghostcells);
        V_S1 = varchunk3D(fID_VS1,VS1import,IMAX,JMAX,KMAX,ghostcells);
        V_S2 = varchunk3D(fID_VS2,VS2import,IMAX,JMAX,KMAX,ghostcells);
        V_S3 = varchunk3D(fID_VS3,VS3import,IMAX,JMAX,KMAX,ghostcells);
        
      % Skip processing for first timestep when there is no plume
        if t==1;
            continue
        end
        
      % Calculate vertical mass flux (space-varied and net) at specified altitudes
        massflux = ROG.*V_G + RO_S1*V_S1 + RO_S2*V_S2 + RO_S3*V_S3;
        netmassflux = squeeze(sum(sum(massflux)));
        netMF_alts = zeros(length(massflux_alts),timesteps);
        netMF_alts(:,t) = netmassflux(massflux_alts);
    
      % Save all calculated mass fluxes (incl. net) at each timestep
        dlmwrite(fullfile(sprintf('%s',dir),sprintf('massflux_all_t%03d.txt',...
            time(t))),massflux,'delimiter','\t');
        dlmwrite(fullfile(sprintf('%s',dir),sprintf('netmassflux_%s.txt',...
            run)),[time(t) netmassflux'],'-append','delimiter','\t');
        
        
      % --------------------- MASS FLUX SLICE FIGURE -------------------- %
        figure(figMFlux)
        cla(axMFlux);
        hMF = slice(axMFlux,1:IMAX-ghostcells,1:KMAX-ghostcells,...
            1:JMAX-ghostcells,massflux,[],[],massflux_alts);
        set(hMF,'FaceColor','interp','EdgeColor','none')
        hEPZ = contourslice(axMFlux,EPG,0,0,massflux_alts,...
            [plumeedge plumeedge]);
        set(hEPZ,'EdgeColor',[1 1 1],'LineWidth',0.5);
        caxis(axMFlux,[massflux_cmin massflux_cmax]);
        tMF = pulsetitle(varMF,PULSE,time,t,titlerun,FREQ);
        title(axMFlux,tMF,'FontWeight','bold');
      % ================================================================= %
      
      
      % --------- SAVE CURRENT FRAMES TO VIDEOS AND IMAGE FILES --------- %
        cd(dir)
        
      % Append current mass flux frame to vidMFlux
        vidfigMF = 'MFluxCurrent.jpg';
        saveas(figMFlux,vidfigMF);
        imgMF = imread(vidfigMF);
        writeVideo(vidMFlux,imgMF);
        
      % If user-specified image filetype is tif, append current timestep
      % frame to multipage tif file. Otherwise, save frame as independent
      % image named by timestep.
        if strcmp(imtype,'tif') == 1 || strcmp(imtype,'tiff') == 1
            imwrite(imgMF,sprintf('MFlux_tsteps_%s.tif',run),'WriteMode','append')
        else
            saveas(figMF,sprintf('MFlux_t%03d_%s.%s',time(t),run,imtype));
        end
      % ================================================================= %
      
    end
  % ===================== E N D   T I M E   L O O P ===================== %
  
  
  % End video write and finish video files
    cd(dir)
    close(vidMFlux);
    
    
  % ------------------ NET MASS FLUX TIME SERIES PLOTS ------------------ %
    if strcmp(PULSE,'T') == 1
      str = sprintf('%s: Unsteady flow %g Hz',titlerun,FREQ);
    elseif strcmp(PULSE,'F') == 1
      str = sprintf('%s: Steady flow',titlerun);
    end
    
    figNetMF = figure('Name','Net Mass Flux','visible',vis,'units',...
        'normalized','outerposition',[0 0 1 1],'PaperPositionMode',...
        'auto','color','w');
    axNetMF = axes('Parent',figNetMF,'box','on','FontSize',12);
    grid(axNetMF,'on');
    hold on
    plot(time,netMF_alts);
    legend(axNetMF,massflux_legend)
    title(axNetMF,sprintf('%s: Net mass flux through specified altitudes',str))
    xlabel(axNetMF,'\bfTime (s)')
    ylabel(axNetMF,'\bfNet mass flux (kg/m^2s)')
    saveas(figNetMF,sprintf('NetMassFlux_tseries_%s.jpg',run))
  % ===================================================================== %
  
  
  cd(postpath)
  disp('Mass flux processing complete.')
  fprintf('vidMFlux_%s has been saved to %s.\n',run,dir)

end


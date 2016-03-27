function [ vidMFlux ] = massFlux3D( runpath,vis,viewaz,viewel,ghostcells,...
    IMAX,JMAX,KMAX,tickx,ticky,tickz,XRES,YRES,ZRES,labelx,labely,labelz,...
    labelXunit,labelYunit,labelZunit,run,timesteps,postpath,massflux_alts,...
    RO_S1,RO_S2,RO_S3,plumeedge,massflux_cmin,massflux_cmax,PULSE,FREQ,...
    time,titlerun,massflux_legend,imtype,savepath,readEPG,fnameEPG,...
    readROG,fnameROG,readVG,fnameVG,readVS1,fnameVS1,readVS2,fnameVS2,...
    readVS3,fnameVS3 )
%massFlux3D Summary of this function goes here
%   Detailed explanation goes here
%
%   Functions called: loadTimestep3D; pulsetitle
%   Last edit: Taryn Black, 21 March 2016

  % Clear directory of appending files from previous processing attempts
    cd(savepath)
    delete('massflux*','netmassflux*','*MFlux*')


  % ----------------------- FIGURE INITIALIZATION ----------------------- %
  % Define variable names for figures
    varMF = 'Mass flux';
    
  % Figure and axes properties
    figMFlux = figure('Name','Mass flux','units','centimeters',...
        'outerposition',[0 0 22 18.75],'visible',vis,'PaperPositionMode',...
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
    cbMFlux = colorbar(axMFlux,'AxisLocation','in','FontSize',12);
    cbMFlux.Label.String = '\bfMass Flux (kg/m^2s)';
    cbMFlux.Ticks = -log10(abs(massflux_cmin)):log10(massflux_cmax);
    cbMFlux.TickLabels = cellstr(num2str([-10.^(abs(-log10(abs...
        (massflux_cmin)):-1)) 0 10.^(1:log10(massflux_cmax))]','%.0e\n'));
    
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
    ROGimport = '%*f%f%*f%*f%*f%*f%*f';
    VGimport  = '%*f%f%*f%*f%*f%*f';
    VS1import = '%f%*f%*f';
    VS2import = '%*f%f%*f';
    VS3import = '%*f%*f%f';
    
  
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
        fID_ROG = fileReadType(fnameROG,readROG,t,runpath,postpath);
        fID_VG = fileReadType(fnameVG,readVG,t,runpath,postpath);
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
        ROG = loadTimestep3D(fID_ROG,ROGimport,readROG,IMAX,JMAX,KMAX,ghostcells);
        V_G  = loadTimestep3D(fID_VG,VGimport,readVG,IMAX,JMAX,KMAX,ghostcells);
        V_S1 = loadTimestep3D(fID_VS1,VS1import,readVS1,IMAX,JMAX,KMAX,ghostcells);
        V_S2 = loadTimestep3D(fID_VS2,VS2import,readVS2,IMAX,JMAX,KMAX,ghostcells);
        V_S3 = loadTimestep3D(fID_VS3,VS3import,readVS3,IMAX,JMAX,KMAX,ghostcells);
        
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
        dlmwrite(fullfile(savepath,sprintf('massflux_all_t%03d.txt',...
            time(t))),massflux,'delimiter','\t');
        dlmwrite(fullfile(savepath,sprintf('netmassflux_%s.txt',...
            run)),[time(t) netmassflux'],'-append','delimiter','\t');
        
      % Separate positive and negative mass fluxes for logarithmic plotting
        logMF = massflux;
        logMF(logMF>0) = log10(logMF(logMF>0));
        logMF(logMF<0) = -log10(abs(logMF(logMF<0)));
        
        
      % --------------------- MASS FLUX SLICE FIGURE -------------------- %
        figure(figMFlux)
        cla(axMFlux);
        hMF = slice(axMFlux,1:IMAX-ghostcells,1:KMAX-ghostcells,...
            1:JMAX-ghostcells,logMF,[],[],massflux_alts);
        set(hMF,'FaceColor','interp','EdgeColor','none')
        hEPZ = contourslice(axMFlux,EPG,0,0,massflux_alts,...
            [plumeedge plumeedge]);
        set(hEPZ,'EdgeColor',[1 1 1],'LineWidth',0.5);
        caxis(axMFlux,[-log10(abs(massflux_cmin)) log10(massflux_cmax)]);
        tMF = pulsetitle(varMF,PULSE,time,t,titlerun,FREQ);
        title(axMFlux,tMF,'FontWeight','bold');
      % ================================================================= %
      
      
      % --------- SAVE CURRENT FRAMES TO VIDEOS AND IMAGE FILES --------- %
        cd(savepath)
        
      % Append current mass flux frame to vidMFlux
        vidfigMF = 'MFluxCurrent.jpg';
        saveas(figMFlux,fullfile(savepath,vidfigMF));
        imgMF = imread(vidfigMF);
        writeVideo(vidMFlux,imgMF);
        
      % If user-specified image filetype is tif, append current timestep
      % frame to multipage tif file. Otherwise, save frame as independent
      % image named by timestep.
        if strcmp(imtype,'tif') == 1 || strcmp(imtype,'tiff') == 1
            imwrite(imgMF,fullfile(savepath,sprintf('MFlux_tsteps_%s.tif',...
                run)),'tif','WriteMode','append')
        else
            saveas(figMF,fullfile(savepath,sprintf('MFlux_t%03d_%s.%s',...
                time(t),run,imtype)));
        end
      % ================================================================= %
      
    end
  % ===================== E N D   T I M E   L O O P ===================== %
  
  
  % End video write and finish video files
    cd(savepath)
    close(vidMFlux);
    
    
  % ------------------ NET MASS FLUX TIME SERIES PLOTS ------------------ %
    if strcmp(PULSE,'T') == 1
      str = sprintf('%s: Unsteady flow %g Hz',titlerun,FREQ);
    elseif strcmp(PULSE,'F') == 1
      str = sprintf('%s: Steady flow',titlerun);
    end
    
    figNetMF = figure('Name','Net Mass Flux','visible',vis,'units',...
        'centimeters','outerposition',[0 0 33.33 18.75],'PaperPositionMode',...
        'auto','color','w');
    axNetMF = axes('Parent',figNetMF,'box','on','FontSize',12);
    grid(axNetMF,'on');
    hold on
    plot(time,netMF_alts);
    legend(axNetMF,massflux_legend)
    title(axNetMF,sprintf('%s: Net mass flux through specified altitudes',str))
    xlabel(axNetMF,'\bfTime (s)')
    ylabel(axNetMF,'\bfNet mass flux (kg/m^2s)')
    saveas(figNetMF,fullfile(savepath,sprintf('NetMassFlux_tseries_%s.jpg',run)))
  % ===================================================================== %
  
  
  cd(postpath)
  disp('Mass flux processing complete.')
  fprintf('vidMFlux_%s has been saved to %s.\n',run,savepath)

end


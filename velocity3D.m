 function [ vidVelo ] = velocity3D( dir,sdistX,sdistY,sdistZ,vis,run,...
     timesteps,postpath,IMAX,JMAX,KMAX,ghostcells,velocity_cmin,...
     velocity_cmax,PULSE,time,titlerun,FREQ,tickx,XRES,labelx,labelXunit,...
     ticky,YRES,labely,labelYunit,tickz,ZRES,labelz,labelZunit,imtype,...
     plumeedge,viewaz,viewel,YGRID,vorticity_cmin,vorticity_cmax,zvort_alts )
%velocity3D calculates the magnitude of gas velocity and plots as a slice
%over time. Also plots vorticity.
%   Detailed explanation goes here
%
%   Functions called: varchunk3D; pulsetitle
%   Last edit: Taryn Black, 18 February 2016

  % Clear directory of appending files from previous processing attempts
    cd(dir)   
    delete('FlowSpeed*','Vorticity*')


  % ----------------------- FIGURE INITIALIZATION ----------------------- %
  % Define variable names for figures
    varU = 'Flow speed';
    varVX = 'Vorticity, \omega_x';
    varVY = 'Vorticity, \omega_y';
    varVZ = 'Vorticity, \omega_z';

  % Subtightplot properties
    gap = [0.03 0.03];
    ht = 0.15;
    wd = 0.15;
    
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
    
  % Figure and axes properties
    figVelo = figure('Name','Flow speed','units','normalized',...
        'outerposition',[0 0 0.4 1],'visible',vis,'PaperPositionMode',...
        'auto','color','w');
    axVelo = axes('Parent',figVelo,'box','on','TickDir','in','FontSize',12);
    hold on
    grid(axVelo,'on');axVelo.Layer = 'top';
    view(axVelo,saz,sel)
    axis(axVelo,'equal',[0,IMAX-ghostcells,0,KMAX-ghostcells,0,...
        JMAX-ghostcells]);
    set(axVelo,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,...
        'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,...
        'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely);
    xlabel(axVelo,sprintf('\\bf Distance_x (%s)',labelXunit))
    ylabel(axVelo,sprintf('\\bf Distance_z (%s)',labelZunit))
    zlabel(axVelo,sprintf('\\bf Altitude (%s)',labelYunit))
    cbVelo = colorbar(axVelo,'AxisLocation','in','FontSize',12);
    cbVelo.Label.String = '\bfFlow Speed (m/s)';

    figVort = figure('Name','Vorticity','units','normalized','visible',...
        vis,'outerposition',[0 0 1 1],'PaperPositionMode','auto',...
        'color','w');
    cd(postpath)
    axVortX = subtightplot(1,3,1,gap,ht,wd);
        hold on
        axis(axVortX,'equal',[0,IMAX-ghostcells,0,KMAX-ghostcells,0,...
            JMAX-ghostcells]);
        set(axVortX,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,...
            'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,...
            'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely);
%         xlabel(axVortX,sprintf('\\bf Distance_x (%s)',labelXunit))
        ylabel(axVortX,sprintf('\\bf Distance_z (%s)',labelZunit))
        zlabel(axVortX,sprintf('\\bf Altitude (%s)',labelYunit))
    axVortY = subtightplot(1,3,2,gap,ht,wd);
        hold on
        axis(axVortY,'equal',[0,IMAX-ghostcells,0,KMAX-ghostcells,0,...
            JMAX-ghostcells]);
        set(axVortY,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,...
            'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,...
            'ZTick',ticky(2:end)/YRES,'ZTickLabel','')
%         xlabel(axVortY,sprintf('\\bf Distance_x (%s)',labelXunit))
%         ylabel(axVortY,sprintf('\\bf Distance_z (%s)',labelZunit))
%         zlabel(axQV,sprintf('\\bf Altitude (%s)',labelYunit))
    axVortZ = subtightplot(1,3,3,gap,ht,wd);
        hold on
        axis(axVortZ,'equal',[0,IMAX-ghostcells,0,KMAX-ghostcells,0,...
            JMAX-ghostcells]);
        set(axVortZ,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,...
            'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,...
            'ZTick',ticky(2:end)/YRES,'ZTickLabel','')
        xlabel(axVortZ,sprintf('\\bf Distance_x (%s)',labelXunit))
%         ylabel(axVortZ,sprintf('\\bf Distance_z (%s)',labelZunit))
        cbVort = colorbar(axVortZ,'Location','eastoutside','AxisLocation',...
           'out','FontSize',12);
        cbVort.Label.String = '\bfFlow Vorticity';
    set([axVortX axVortY axVortZ],'box','on','TickDir','in','FontSize',12)
    grid(axVortX,'on'); grid(axVortY,'on'); grid(axVortZ,'on');
    view(axVortX,viewaz,viewel); view(axVortY,viewaz,viewel); 
        view(axVortZ,viewaz,viewel)
    cd(dir)
    
  % Initialize video
    vidVelo = VideoWriter(sprintf('vidVelo_%s.avi',run));
    vidVelo.Quality = 100;
    vidVelo.FrameRate = 10;
    open(vidVelo);
    set(gcf,'Visible',vis);
    
    vidVort = VideoWriter(sprintf('vidVort_%s.avi',run));
    vidVort.Quality = 100;
    vidVort.FrameRate = 10;
    open(vidVort);
    set(gcf,'Visible',vis);
  % ===================================================================== %
  
  
  % File import specifications: columns to read or skip for each variable
    EGimport = '%f%*f%*f%*f%*f%*f%*f';
    UGimport = '%f%*f%*f%*f%*f%*f';
    VGimport = '%*f%f%*f%*f%*f%*f';
    WGimport = '%*f%*f%f%*f%*f%*f';
    
    
  % =================== B E G I N   T I M E   L O O P =================== %
    t = 0;
    while t <= timesteps
        
        t = t+1;
        
      % Queue up current timestep files
        cd(dir)
        fclose('all');
        clear fID*;
        fID_EPG = fopen(sprintf('EP_t%02d.txt',t));
        fID_UG = fopen(sprintf('U_G_t%02d.txt',t));
        fID_VG = fopen(sprintf('U_G_t%02d.txt',t));
        fID_WG = fopen(sprintf('U_G_t%02d.txt',t));
        cd(postpath)
        
      % Prepare velocities for full domain at current timestep
        try
            EPG = varchunk3D(fID_EPG,EGimport,IMAX,JMAX,KMAX,ghostcells);
        catch ME
            warning('Error in varchunk3D at t=%d s:\n%s\nContinuing to next simulation.',...
                time(t),ME.identifier)
            break
        end
        U_G = varchunk3D(fID_UG,UGimport,IMAX,JMAX,KMAX,ghostcells);
        V_G = varchunk3D(fID_VG,VGimport,IMAX,JMAX,KMAX,ghostcells);
        W_G = varchunk3D(fID_WG,WGimport,IMAX,JMAX,KMAX,ghostcells);
        
      % Skip processing for first timestep when there is no plume
        if t==1;
            continue
        end
        
      % Calculate magnitude of velocity everywhere
        flowspeed = sqrt(U_G.^2 + V_G.^2 + W_G.^2);
        
      
      % -------------------- FLOW SPEED SLICE FIGURE -------------------- %
        figure(figVelo)
        cla(axVelo);
        hFS = slice(0.5:(IMAX-ghostcells-0.5),0.5:(KMAX-ghostcells-0.5),...
            0.5:(JMAX-ghostcells-0.5),flowspeed,...
            sdistX*(IMAX-ghostcells),sdistY*(KMAX-ghostcells),...
            sdistZ*(JMAX-ghostcells));
        hFS.FaceColor = 'interp';
        hFS.EdgeColor = 'none';
        caxis(axVelo,[velocity_cmin velocity_cmax]);
        tFS = pulsetitle(varU,PULSE,time,t,titlerun,FREQ);
        title(tFS,'FontSize',12,'FontWeight','bold');
      % ================================================================= %
      
      
      % --------------------- OVERLAY PLUME OUTLINE --------------------- %
        hEP = contourslice(EPG,sdistX*(IMAX-ghostcells),...
            sdistY*(KMAX-ghostcells),0,[plumeedge plumeedge]);
        set(hEP,'EdgeColor',[1 1 1],'LineWidth',0.5);
      % ================================================================= %
      
      
      % ----------------- CALCULATE AND PLOT VORTICITY ------------------ %
        [curlx,curly,curlz] = curl(U_G,V_G,W_G);
        
        figure(figVort)
        cla(axVortX);cla(axVortY);cla(axVortZ);
        
        hVX = slice(axVortX,0.5:(IMAX-ghostcells-0.5),...
            0.5:(KMAX-ghostcells-0.5),0.5:(JMAX-ghostcells-0.5),curlx,...
            0.5*(IMAX-ghostcells),0.5*(KMAX-ghostcells),[]);
        set(hVX,'FaceColor','interp','EdgeColor','none') 
        hEPX = contourslice(axVortX,EPG,0.5*(IMAX-ghostcells),...
            0.5*(KMAX-ghostcells),0,[plumeedge plumeedge]);
        set(hEPX,'EdgeColor',[1 1 1],'LineWidth',0.5);   
        hVXv = plot3(axVortX,0.5*(IMAX-ghostcells)*ones(1,length(YGRID)),...
            0.5*(KMAX-ghostcells)*ones(1,length(YGRID)),...
            0.5:(JMAX-ghostcells-0.5),'k');
        caxis(axVortX,[vorticity_cmin vorticity_cmax]);
        tVX = pulsetitle(varVX,PULSE,time,t,titlerun,FREQ);
        title(axVortX,tVX,'FontSize',12,'FontWeight','bold');
          
        hVY = slice(axVortY,0.5:(IMAX-ghostcells-0.5),...
            0.5:(KMAX-ghostcells-0.5),0.5:(JMAX-ghostcells-0.5),curly,...
            0.5*(IMAX-ghostcells),0.5*(KMAX-ghostcells),[]);
        set(hVY,'FaceColor','interp','EdgeColor','none')
        hEPY = contourslice(axVortY,EPG,0.5*(IMAX-ghostcells),...
            0.5*(KMAX-ghostcells),0,[plumeedge plumeedge]);
        set(hEPY,'EdgeColor',[1 1 1],'LineWidth',0.5); 
        hVYv = plot3(axVortY,0.5*(IMAX-ghostcells)*ones(1,length(YGRID)),...
            0.5*(KMAX-ghostcells)*ones(1,length(YGRID)),...
            0.5:(JMAX-ghostcells-0.5),'k');
        caxis(axVortY,[vorticity_cmin vorticity_cmax]);
        tVY = pulsetitle(varVY,PULSE,time,t,titlerun,FREQ);
        title(axVortY,tVY,'FontWeight','bold');
        
        hVZ = slice(axVortZ,1:IMAX-ghostcells,1:KMAX-ghostcells,...
            1:JMAX-ghostcells,curlz,[],[],zvort_alts);
        set(hVZ,'FaceColor','interp','EdgeColor','none')
        hEPZ = contourslice(axVortZ,EPG,0,0,zvort_alts,...
            [plumeedge plumeedge]);
        set(hEPZ,'EdgeColor',[1 1 1],'LineWidth',0.5); 
        caxis(axVortZ,[vorticity_cmin vorticity_cmax]);
        tVZ = pulsetitle(varVZ,PULSE,time,t,titlerun,FREQ);
        title(axVortZ,tVZ,'FontWeight','bold');
        
        PosS1 = get(axVortX,'position');
        PosS2 = get(axVortY,'position');
        PosS3 = get(axVortZ,'position');
        PosS2(3:4) = PosS1(3:4);
        PosS3(3:4) = PosS1(3:4);
        set(axVortY,'position',PosS2);
        set(axVortZ,'position',PosS3);
      % ================================================================= %
      
      
      % --------- SAVE CURRENT FRAMES TO VIDEOS AND IMAGE FILES --------- %
        cd(dir)
        
      % Append current flow speed frame to vidVelo, vidvort
        vidfigvelo = 'FlowSpeedCurrent.jpg';
        saveas(figVelo,vidfigvelo);
        imgvelo = imread(vidfigvelo);
        writeVideo(vidVelo,imgvelo);
        
        vidfigvort = 'VorticityCurrent.jpg';
        saveas(figVort,vidfigvort);
        imgvort = imread(vidfigvort);
        writeVideo(vidVort,imgvort);
        
      % If user-specified image filetype is tif, append current timestep
      % frame to multipage tif file. Otherwise, save frame as independent
      % image named by timestep.
        if strcmp(imtype,'tif') == 1 || strcmp(imtype,'tiff') == 1
            imwrite(imgvelo,sprintf('FlowSpeed_tsteps_%s.tif',run),'WriteMode','append')
            imwrite(imgvort,sprintf('Vorticity_tsteps_%s.tif',run),'WriteMode','append')
        else
            saveas(figVelo,sprintf('FlowSpeed_%03ds_%s.%s',time(t),run,imtype));
            saveas(figVort,sprintf('Vorticity_%03ds_%s.%s',time(t),run,imtype));
        end 
      % ================================================================= %
                    
    end
  % ===================== E N D   T I M E   L O O P ===================== %
    
  
  % End video write and finish video files
    cd(dir)
    close(vidVelo)
    close(vidVort)
    
    cd(postpath)
    disp('Flow speed processing complete.')
    fprintf('vidVelo_%s has been saved to %s.\n',run,dir)
    fprintf('vidVort_%s has been saved to %s.\n',run,dir)

end

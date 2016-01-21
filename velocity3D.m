 function [ vidVelo ] = velocity3D( dir,sdistX,sdistY,sdistZ,vis,run,...
     timesteps,postpath,IMAX,JMAX,KMAX,ghostcells,velocity_cmin,...
     velocity_cmax,PULSE,time,titlerun,FREQ,tickx,XRES,labelx,labelXunit,...
     ticky,YRES,labely,labelYunit,tickz,ZRES,labelz,labelZunit,imtype,plumeedge )
%velocity3D calculates the magnitude of gas velocity and plots as a slice
%over time.
%   Detailed explanation goes here
%
%   Functions called:
%   Last edit: Taryn Black, 21 January 2016

  % Clear directory of appending files from previous processing attempts
    cd(dir)
    delete('FlowSpeed_*')


  % ----------------------- FIGURE INITIALIZATION ----------------------- %
  % Define variable names for figures
    varU = 'Flow speed';
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
    
  % Figure and axes properties
    figVelo = figure('Name','Flow speed','units','normalized',...
        'outerposition',[0 0 0.5 1],'visible',vis);
    set(figVelo,'color','w')
    axVelo = axes('Parent',figVelo);
    hold on
    set(axVelo,'box','on','FontSize',12)
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
    
  % Initialize video
    vidVelo = VideoWriter(sprintf('vidVelo_%s.avi',run));
    vidVelo.Quality = 100;
    vidVelo.FrameRate = 10;
    open(vidVelo);
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
        hc = colorbar;
        caxis([velocity_cmin velocity_cmax]);
        ylabel(hc,'\bf Flow Speed [m/s]','FontSize',12)
        tL = pulsetitle(varU,PULSE,time,t,titlerun,FREQ);
        title(tL,'FontSize',12,'FontWeight','bold');
      % ================================================================= %
      
      
      % --------------------- OVERLAY PLUME OUTLINE --------------------- %
        hEP = contourslice(EPG,sdistX*(IMAX-ghostcells),sdistY*(KMAX-ghostcells),0,[plumeedge plumeedge]);
        set(hEP,'EdgeColor',[1 1 1],'LineWidth',0.5);
      % ================================================================= %
      
      
      % --------- SAVE CURRENT FRAMES TO VIDEOS AND IMAGE FILES --------- %
        cd(dir)
        
      % Append current flow speed frame to vidVelo
        vidfig = 'FlowSpeedCurrent.jpg';
        saveas(figVelo,vidfig);
        img = imread(vidfig);
        writeVideo(vidVelo,img);
        
      % If user-specified image filetype is tif, append current timestep
      % frame to multipage tif file. Otherwise, save frame as independent
      % image named by timestep.
        if strcmp(imtype,'tif') == 1 || strcmp(imtype,'tiff') == 1
            imwrite(img,sprintf('FlowSpeed_tsteps_%s.tif',run),'WriteMode','append')
        else
            saveas(figVelo,sprintf('FlowSpeed_%03ds_%s.%s',time(t),run,imtype));
        end 
      % ================================================================= %
                    
    end
  % ===================== E N D   T I M E   L O O P ===================== %
    
  
  % End video write and finish video files
    cd(dir)
    close(vidVelo)
    
    cd(postpath)
    disp('Flow speed processing complete.')
    fprintf('vidVelo_%s has been saved to %s.\n',run,dir)

    end

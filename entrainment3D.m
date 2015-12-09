function [ vidEntr ] = entrainment3D( run,dir,vis,ghostcells,IMAX,JMAX,...
    KMAX,tickx,labelx,labelxunit,ticky,labely,labelyunit,tickz,labelz,...
    labelzunit,plumeedge,XRES,YRES,ZRES,postpath,PULSE,FREQ,time,...
    vel_char,entrainment_cmin,entrainment_cmax,viewaz,viewel,imtype,...
    titlerun,timesteps,isoEPG,colEPG,trnEPG,DT,VENT_R )
%entrainment3D Summary of this function goes here
%   entrainment3D ---does things---
%
%   Functions called: varchunk3D; pulsetitle
%   Last edit: Taryn Black, 9 December 2015
    
  % Clear directory of appending files from previous processing attempts
    cd(dir)
    delete('coeff_avg-std*','entr_avg-std*','expn_avg-std*',...
        'plot_time*','ecoeff_all*','plumevolume_*','plumeheight_*',...
        'e_MortonConic_*','Entr_t*','EPG_t*','Quiver_t*');
    
    
  % ----------------------- FIGURE INITIALIZATION ----------------------- %
  % Define variable names for figures
    varEP = 'Gas volume fraction';
    varEn = 'Entrainment';
    cd(dir)

  % Gas volume fraction isosurface: figure and axes properties    
    figEP = figure('Name','Gas Volume Fraction','visible',vis,'units',...
        'normalized','outerposition',[0 0 0.5 1]);
    hold on
    view(viewaz,viewel)
    axis equal
    axis([ghostcells-1,IMAX-(ghostcells/2),ghostcells-1,...
        KMAX-(ghostcells/2),ghostcells-1,JMAX-(ghostcells/2)]);
    set(gca,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,'FontSize',12)
      xlabel(sprintf('\\bf Distance (%s)',labelxunit),'FontSize',12)
    set(gca,'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,'FontSize',12)
      ylabel(sprintf('\\bf Distance (%s)',labelzunit),'FontSize',12)
    set(gca,'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely,'FontSize',12)
      zlabel(sprintf('\\bf Altitude (%s)',labelyunit),'FontSize',12)
    grid on
    box on
    
  % Isonormal/velocity quivers: figure and axes properties
    figQ = figure('Name','Isonormals and Velocities','visible',vis,...
        'units','normalized','outerposition',[0 0 1 1]);
    hold on
    subfigQN = subplot(1,2,1);
      hold on
      view(viewaz,viewel)
      axis equal
      axis([ghostcells-1,IMAX-(ghostcells/2),ghostcells-1,...
          KMAX-(ghostcells/2),ghostcells-1,JMAX-(ghostcells/2)]);
      set(gca,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,'FontSize',12)
        xlabel(sprintf('\\bf Distance_x (%s)',labelxunit),'FontSize',12)
      set(gca,'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,'FontSize',12)
        ylabel(sprintf('\\bf Distance_z (%s)',labelzunit),'FontSize',12)
      set(gca,'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely,'FontSize',12)
        zlabel(sprintf('\\bf Altitude (%s)',labelyunit),'FontSize',12)
      grid on
      box on
    subfigQV = subplot(1,2,2);
      hold on
      view(viewaz,viewel)
      axis equal
      axis([ghostcells-1,IMAX-(ghostcells/2),ghostcells-1,...
          KMAX-(ghostcells/2),ghostcells-1,JMAX-(ghostcells/2)]);
      set(gca,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,'FontSize',12)
        xlabel(sprintf('\\bf Distance_x (%s)',labelxunit),'FontSize',12)
      set(gca,'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,'FontSize',12)
        ylabel(sprintf('\\bf Distance_z (%s)',labelzunit),'FontSize',12)
      set(gca,'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely,'FontSize',12)
        zlabel(sprintf('\\bf Altitude (%s)',labelyunit),'FontSize',12)
      grid on
      box on
        
  % Entrainment isosurface: figure and axes properties
    figEn = figure('Name','Entrainment','visible',vis,'units',...
        'normalized','outerposition',[0.5 0 0.5 1]);
    hold on
    view(viewaz,viewel)
    axis equal
    axis([ghostcells-1,IMAX-(ghostcells/2),ghostcells-1,...
        KMAX-(ghostcells/2),ghostcells-1,JMAX-(ghostcells/2)]);
    set(gca,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,'FontSize',12)
      xlabel(sprintf('\\bf Distance (%s)',labelxunit),'FontSize',12)
    set(gca,'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,'FontSize',12)
      ylabel(sprintf('\\bf Distance (%s)',labelzunit),'FontSize',12)
    set(gca,'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely,'FontSize',12)
      zlabel(sprintf('\\bf Altitude (%s)',labelyunit),'FontSize',12)
    grid on
    box on
    
  % Gas volume fraction isosurface: video
    vidEPG = VideoWriter(sprintf('vidEPG_%s.avi',run));
    vidEPG.Quality = 100;
    vidEPG.FrameRate = 10;
    open(vidEPG);
    set(gcf,'Visible',vis);
    
  % Isonormals/velocities quiver: video
    vidQ = VideoWriter(sprintf('vidQuiver_%s.avi',run));
    vidQ.Quality = 100;
    vidQ.FrameRate = 10;
    open(vidQ)
    set(gcf,'Visible',vis);
    
  % Entrainment isosurface: video
    vidEntr = VideoWriter(sprintf('vidEntr_%s.avi',run));
    vidEntr.Quality = 100;
    vidEntr.FrameRate = 10;
    open(vidEntr);
    set(gcf,'Visible',vis);
  % ===================================================================== %

  
  % File import specifications: columns to read or skip for each variable
    EGimport = '%f%*f%*f%*f%*f%*f%*f';
    UGimport = '%f%*f%*f%*f%*f%*f';
    VGimport = '%*f%f%*f%*f%*f%*f';
    WGimport = '%*f%*f%f%*f%*f%*f';
    
  % Preallocate vectors for text files
    avg_coeff = zeros(1,timesteps);
    std_coeff = zeros(1,timesteps);
    avg_entr  = zeros(1,timesteps);
    std_entr  = zeros(1,timesteps);
    avg_expn  = zeros(1,timesteps);
    std_expn  = zeros(1,timesteps);
    plumevolume = zeros(1,timesteps);
    plumeheight = zeros(1,timesteps);
    plumetoprad = zeros(1,timesteps);
    e_Mconic = zeros(1,timesteps);
     
    
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
        
      % Prepare EPG and velocities for full domain at current timestep
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
                
        
      % ------------ GAS VOLUME FRACTION ISOSURFACES FIGURES ------------ %
      % Find all specified isosurfaces and plot
        figure(figEP)
        cla;
        for j = 1:length(isoEPG)
            surf(j) = patch(isosurface(EPG,isoEPG(j)));
            set(surf(j),'FaceColor',colEPG(j,:),'EdgeColor','none',...
                'FaceAlpha',trnEPG(j));
            hold on
        end
        camlight('right')
        camlight('left')
        lighting gouraud
        tLEP = pulsetitle(varEP,PULSE,time,t,titlerun,FREQ);
        title(tLEP,'FontSize',12,'FontWeight','bold');
      % ================================================================= %

      
      % ---------- FIND PLUME ISONORMALS AND SURFACE VELOCITIES --------- %
      % Find plume boundary (isosurface)
        plumesurf = isosurface(EPG,plumeedge);
        
      % Extract coordinates of plume isosurface vertices
        plumeX = plumesurf.vertices(:,1);
        plumeY = plumesurf.vertices(:,2);
        plumeZ = plumesurf.vertices(:,3);
        
      % Calculate normals to plume surface (*-1 to point out of plume)
        plumenorm = -isonormals(EPG,plumesurf.vertices)';
        
      % Convert plume normals to unit normals
        unitnorm = zeros(3,length(plumenorm));
        for i = 1:length(plumenorm)
            unitnorm(:,i) = plumenorm(:,i)./norm(plumenorm(:,i));
        end
        
      % Interpolate for velocities at plume surface vertices
        plumeU = interp3(U_G,plumeX,plumeY,plumeZ);
        plumeV = interp3(V_G,plumeX,plumeY,plumeZ);
        plumeW = interp3(W_G,plumeX,plumeY,plumeZ);
        plumevelocity = [plumeU plumeW plumeV]';
        
      % Calculate magnitudes and directional components of plume unit
      % normal velocities (PUNV)
        PUNV_mag = dot(plumevelocity,unitnorm); 
        PUNV_X = PUNV_mag.*unitnorm(1,:);
        PUNV_Y = PUNV_mag.*unitnorm(2,:);
        PUNV_Z = PUNV_mag.*unitnorm(3,:); 
      % ================================================================= %
      
      
      % --------------- ISONORMAL/VELOCITY QUIVER FIGURES --------------- %  
      % Quiver plot of isonormals and velocities
        figure(figQ)
        q = 50; % reducement factor for quiver plot
        subfigQN = subplot(1,2,1);
          cla;
          hold on
          view(viewaz,viewel)
          axis equal
          axis([ghostcells-1,IMAX-(ghostcells/2),ghostcells-1,...
              KMAX-(ghostcells/2),ghostcells-1,JMAX-(ghostcells/2)]);
          set(gca,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,'FontSize',12)
            xlabel(sprintf('\\bf Distance_x (%s)',labelxunit),'FontSize',12)
          set(gca,'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,'FontSize',12)
            ylabel(sprintf('\\bf Distance_z (%s)',labelzunit),'FontSize',12)
          set(gca,'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely,'FontSize',12)
            zlabel(sprintf('\\bf Altitude (%s)',labelyunit),'FontSize',12)
          grid on
          box on
          quiver3(plumeX(1:q:length(plumeX))',plumeY(1:q:length(plumeY))',...
              plumeZ(1:q:length(plumeZ))',unitnorm(1,1:q:length(unitnorm)),...
              unitnorm(2,1:q:length(unitnorm)),unitnorm(3,1:q:length(unitnorm)),...
              'MaxHeadSize',10,'AutoScaleFactor',1,'LineWidth',0.1);
          title('Plume surface isonormals')
        subfigQV = subplot(1,2,2);
          cla;
          hold on
          view(viewaz,viewel)
          axis equal
          axis([ghostcells-1,IMAX-(ghostcells/2),ghostcells-1,...
              KMAX-(ghostcells/2),ghostcells-1,JMAX-(ghostcells/2)]);
          set(gca,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,'FontSize',12)
            xlabel(sprintf('\\bf Distance_x (%s)',labelxunit),'FontSize',12)
          set(gca,'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,'FontSize',12)
            ylabel(sprintf('\\bf Distance_z (%s)',labelzunit),'FontSize',12)
          set(gca,'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely,'FontSize',12)
            zlabel(sprintf('\\bf Altitude (%s)',labelyunit),'FontSize',12)
          grid on
          box on
          quiver3(plumeX(1:q:length(plumeX))',plumeY(1:q:length(plumeY))',...
              plumeZ(1:q:length(plumeZ))',PUNV_X(1:q:length(PUNV_X)),...
              PUNV_Y(1:q:length(PUNV_Y)),PUNV_Z(1:q:length(PUNV_Z)),...
              'MaxHeadSize',20,'AutoScaleFactor',5,'LineWidth',0.1);
          title('Velocity magnitudes')
        PosQN = get(subfigQN,'position');
        PosQV = get(subfigQV,'position');
        PosQV(3:4) = PosQN(3:4);
        set(subfigQV,'position',PosQV);
      % ================================================================= %
      
      
      % ----------- CALCULATE ENTRAINMENT AND SPATIAL AVG/STD ----------- %
      % Calculate entrainment(-)/expansion(+) coefficient from velocities
        e_coeff = PUNV_mag/vel_char;
        entr = e_coeff(e_coeff<0);
        expn = e_coeff(e_coeff>0);
        
      % Calculate spatially-averaged total coefficient, entrainment, and
      % expansion, and their standard deviations at each timestep
        avg_coeff(t) = mean(e_coeff);
        std_coeff(t) = std(e_coeff);
        avg_entr(t) = mean(entr);
        std_entr(t) = std(entr);
        avg_expn(t) = mean(expn);
        std_expn(t) = std(expn);
        
      % Save all calculated entrainment values at each timestep
        dlmwrite(fullfile(sprintf('%s',dir),sprintf('coeff_avg-std_%s.txt',...
            run)),[avg_coeff(t) std_coeff(t)],'-append','delimiter','\t',...
            'precision','%0.6f');
        dlmwrite(fullfile(sprintf('%s',dir),sprintf('entr_avg-std_%s.txt',...
            run)),[avg_entr(t) std_entr(t)],'-append','delimiter','\t',...
            'precision','%0.6f');
        dlmwrite(fullfile(sprintf('%s',dir),sprintf('expn_avg-std_%s.txt',...
            run)),[avg_expn(t) std_expn(t)],'-append','delimiter','\t',...
            'precision','%0.6f');
        dlmwrite(fullfile(sprintf('%s',dir),sprintf('plot_time_%s.txt',...
            run)),time(t)','-append','delimiter','\t','precision','%g');
        dlmwrite(fullfile(sprintf('%s',dir),sprintf('ecoeff_all_t%03d.txt',...
            time(t))),e_coeff,'delimiter','\t','precision','%0.6f');
      % ================================================================= %
      
      
      % ---------------- ENTRAINMENT ISOSURFACE FIGURES ----------------- %
      % Define concatenated colormap of entrainment so cmaps meet at e=0
        nmap = 256;
        colormap([winter(round(nmap*-entrainment_cmin/...
            (entrainment_cmax-entrainment_cmin)));...
            flipud(autumn(nmap-round(nmap*-entrainment_cmin/...
            (entrainment_cmax-entrainment_cmin))))]);
        caxis([entrainment_cmin entrainment_cmax])
        cmap = colormap;
        emap = linspace(entrainment_cmin,entrainment_cmax,nmap);
        e_color = zeros(length(e_coeff),3);   
        for k = 1:length(e_coeff)
            if e_coeff(k) > entrainment_cmax
                e_coeff(k) = entrainment_cmax;
            elseif e_coeff(k) < entrainment_cmin
                e_coeff(k) = entrainment_cmin;
            end
        end
        e_round = interp1(emap,emap,e_coeff,'nearest');
        emap = [emap; 1:nmap];
        for j = 1:length(e_coeff)
            e_color(j,:) = cmap(emap(2,emap(1,:) == e_round(j)),:);
        end
        colormap(figEn,cmap)

      % Plot entrainment on plume surface
        figure(figEn)
        cla;
        patch('Vertices',[plumeX plumeY plumeZ],'Faces',1:length(plumeX),...
            'FaceVertexCData',e_color,'FaceColor','none','EdgeColor',...
            'none','Marker','o','MarkerFaceColor','flat')
        colorbar
        caxis([entrainment_cmin entrainment_cmax])
        tLEn = pulsetitle(varEn,PULSE,time,t,titlerun,FREQ);
        tlEn2 = sprintf('Characteristic velocity: %g m/s',vel_char);
        title([tLEn;tlEn2],'FontSize',12,'FontWeight','bold');
        camlight('right')
        camlight('left')
        lighting gouraud
      % ================================================================= %
      
      
      % -------------- ENTRAINMENT COMPARISON CALCULATIONS -------------- %
      % Calculate plume volume [m3], height [m], and conic max radius [m]
        numcells = sum(sum(sum(EPG <= plumeedge)));
        plumevolume(t) = numcells*XRES*YRES*ZRES;
        plumeheight(t) = max(plumeZ)*ZRES;
        plumetoprad(t) = sqrt(3*plumevolume(t)/(pi*plumeheight(t)))-VENT_R;
        dlmwrite(fullfile(sprintf('%s',dir),sprintf('plumevolume_%s.txt',...
            run)),plumevolume(t),'-append','delimiter','\t','precision','%g');
        dlmwrite(fullfile(sprintf('%s',dir),sprintf('plumeheight_%s.txt',...
            run)),plumeheight(t),'-append','delimiter','\t','precision','%g');
        
      % Calculate entrainment coefficient using Morton linear assumption
      % and using conic plume simplification
        e_Morton = PUNV_mag./plumeV';
        e_Mconic(t) = (plumetoprad(t) - plumetoprad(t-1))/(DT*vel_char);
        dlmwrite(fullfile(sprintf('%s',dir),sprintf('e_MortonConic_%s.txt',...
            run)),e_Mconic(t),'-append','delimiter','\t','precision','%g');
      % ================================================================= %
      
      
      % --------- SAVE CURRENT FRAMES TO VIDEOS AND IMAGE FILES --------- %
        cd(dir)
        
      % Append current gas volume fraction frame to vidEPG
        vidfigEP = 'EPGCurrent.jpg';
        saveas(figEP,vidfigEP);
        imgEP = imread(vidfigEP);
        writeVideo(vidEPG,imgEP);

      % Append current isonormal/velocity quiver frame to vidQ
        vidfigQ = 'QuiverCurrent.jpg';
        saveas(figQ,vidfigQ);
        imgQ = imread(vidfigQ);
        writeVideo(vidQ,imgQ);
        
      % Append current entrainment frame to vidEntr
        vidfigEn = 'EntrCurrent.jpg';
        saveas(figEn,vidfigEn);
        imgEn = imread(vidfigEn);
        writeVideo(vidEntr,imgEn);
            
      % If user-specified image filetype is tif, append current timestep
      % frame to multipage tif file. Otherwise, save frame as independent
      % image named by timestep.
        if strcmp(imtype,'tif') == 1 || strcmp(imtype,'tiff') == 1
            imwrite(imgEn,sprintf('Entr_tsteps_%s.tif',run),'WriteMode','append')
            imwrite(imgQ,sprintf('Quiver_tsteps_%s.tif',run),'WriteMode','append')
            imwrite(imgEP,sprintf('EPG_tsteps_%s.tif',run),'WriteMode','append')
        else
            saveas(figEn,sprintf('Entr_t%03d_%s.%s',time(t),run,imtype));
            saveas(figQ,sprintf('Quiver_t%03d_%s.%s',time(t),run,imtype));
            saveas(figEP,sprintf('EPG_t%03d_%s.%s',time(t),run,imtype));
        end
      % ================================================================= %
            
    end
  % ===================== E N D   T I M E   L O O P ===================== %
    
    
  % End video write and finish video files
    cd(dir)
    close(vidEPG);
    close(vidQ);
    close(vidEntr);
      
      
  % ------------------- ENTRAINMENT TIME SERIES PLOTS ------------------- %
  % Total plume volume and change in plume volume
    if strcmp(PULSE,'T') == 1
      str = sprintf('%s: Unsteady flow %g Hz',titlerun,FREQ);
    elseif strcmp(PULSE,'F') == 1
      str = sprintf('%s: Steady flow',titlerun);
    end
    fig_plumevol = figure('Name','Plume Volume','visible',vis);
    hvol1 = subplot(2,1,1);
      plot(time,plumevolume)
      title(hvol1,{sprintf('%s: Total plume volume',str)},...
          'FontWeight','bold','FontSize',10)
      xlabel(hvol1,{'Time (s)'},'FontWeight','bold','FontSize',10)
      ylabel(hvol1,{'Volume (m^3)'},'FontWeight','bold','FontSize',10)
    hvol2 = subplot(2,1,2);
      plot(time(2:length(plumevolume)),diff(plumevolume))
      title(hvol2,{sprintf('%s: Change in plume volume',str)},...
          'FontWeight','bold','FontSize',10)
      xlabel(hvol2,{'Time (s)'},'FontWeight','bold','FontSize',10)
      ylabel(hvol2,{'\DeltaVolume (m^3)'},'FontWeight','bold','FontSize',10)
    saveas(fig_plumevol,sprintf('PlumeVolume_%s.jpg',run));
    
  % Plume-averaged entrainment/expansion
    fig_coeff = figure('Name','Entrainment Coefficients','visible',vis);
    hold on
    errorbar(time,avg_coeff,std_coeff,'k')
    errorbar(time,avg_entr,std_entr,'b')
    errorbar(time,avg_expn,std_expn,'r')
    title(sprintf('%s: Plume-averaged coefficients',str),'FontWeight',...
        'bold','FontSize',10)
    xlabel('Time (s)','FontWeight','bold','FontSize',10)
    ylim([-0.5 0.5])
    legend({'Total coefficient','Entrainment','Expansion'},'Box','on',...
        'Location','EastOutside','FontWeight','bold','FontSize',10)
    saveas(fig_coeff,sprintf('Coefficients_%s.jpg',run));
    
  % Comparison conic coefficient
    fig_Morton = figure('Name','Conic entrainment coefficient','visible',vis);
    plot(time,e_Mconic)
    title(sprintf('%s: Morton conic entrainment coefficient',str),...
        'FontWeight','bold','FontSize',10)
    xlabel('Time (s)','FontWeight','bold','FontSize',10)
    ylabel('Coefficient','FontWeight','bold','FontSize',10)
    ylim([0 0.5])
    saveas(fig_Morton,sprintf('MortonConic_%s.jpg',run));
  % =================================================================== %
  
  
  cd(postpath)
  disp('Entrainment processing complete.')
  fprintf('vidEPG_%s has been saved to %s.\n',run,dir)
  fprintf('vidQ_%s has been saved to %s.\n',run,dir)
  fprintf('vidEntr_%s has been saved to %s.\n',run,dir)

end

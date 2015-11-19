function [ vidEntr ] = entrainment3D( run,dir,vis,ghostcells,IMAX,...
    JMAX,KMAX,tickx,labelx,labelxunit,ticky,labely,labelyunit,tickz,...
    labelz,labelzunit,plumeedge,XRES,YRES,ZRES,postpath,PULSE,FREQ,...
    time,vel_char,entrainment_cmin,entrainment_cmax,viewaz,viewel,...
    imtype,titlerun,timesteps,isoEPG,colEPG,trnEPG,DT,VENT_R )
%entrainment3D Summary of this function goes here
%   entrainment3D ---does things---
%
%   Special functions called: varchunk3D; pulsetitle
%   Last edit: Taryn Black, 18 November 2015

    varEP = 'Gas volume fraction';
    varEn = 'Entrainment';
    
%%% Clear directory of appending files from previous processing attempts
    cd(dir)
    delete('*.txt','Entr_t*','EPG_t*');
    
%%% Initialize figure frames
    cd(dir)
    
    figEP = figure('Name','Gas Volume Fraction','visible',vis,'units','normalized','outerposition',[0 0 0.5 1]);
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
    
    figQ = figure('Name','Isonormals and Velocities','visible',vis,'units','normalized','outerposition',[0 0 1 1]);
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
            
    figEn = figure('Name','Entrainment','visible',vis,'units','normalized','outerposition',[0.5 0 0.5 1]);
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
    
%%% Initialize videos
    vidEPG = VideoWriter(sprintf('vidEPG_%s.avi',run));
    vidEPG.Quality = 100;
    vidEPG.FrameRate = 10;
    open(vidEPG);
    set(gcf,'Visible',vis);
    
    vidQ = VideoWriter(sprintf('vidQuiver_%s.avi',run));
    vidQ.Quality = 100;
    vidQ.FrameRate = 10;
    open(vidQ)
    set(gcf,'Visible',vis);
    
    vidEntr = VideoWriter(sprintf('vidEntr_%s.avi',run));
    vidEntr.Quality = 100;
    vidEntr.FrameRate = 10;
    open(vidEntr);
    set(gcf,'Visible',vis);

%%% Calculate entrainment at each timestep, plot, and save video.
    fID_EPG = fopen(sprintf('EP_G_%s',run));
    fID_UG  = fopen(sprintf('U_G_%s',run));
    fID_VG  = fopen(sprintf('V_G_%s',run));
    fID_WG  = fopen(sprintf('W_G_%s',run));
        
    t = 0;
    while t <= timesteps 
        
        t = t+1;
        
        cd(postpath)
      % Find gas volume fraction and velocities for entire domain at current timestep
        try
            EPG = varchunk3D(fID_EPG,IMAX,JMAX,KMAX,ghostcells);
        catch ME
            warning('Error in varchunk3D at t=%d s:\n%s\nContinuing to next simulation.',time(t),ME.identifier)
            break
        end
        U_G = varchunk3D(fID_UG,IMAX,JMAX,KMAX,ghostcells);
        V_G = varchunk3D(fID_VG,IMAX,JMAX,KMAX,ghostcells);
        W_G = varchunk3D(fID_WG,IMAX,JMAX,KMAX,ghostcells);
        
        % Skip processing for first timestep when there is no plume.
          if t==1;
              continue
          end
                
      % Find all specified gas volume fraction isosurfaces and plot.
        figure(figEP)
        cla;
        for j = 1:length(isoEPG)
            surf(j) = patch(isosurface(EPG,isoEPG(j)));
            set(surf(j),'FaceColor',colEPG(j,:),'EdgeColor','none','FaceAlpha',trnEPG(j));
            hold on
        end
        camlight('right')
        camlight('left')
        lighting gouraud
        
        tLEP = pulsetitle(varEP,PULSE,time,t,titlerun,FREQ);
        title(tLEP,'FontSize',12,'FontWeight','bold');

      % Find plume boundary (isosurface) and plot. 
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
        
      % Quiver plot of isonormals and velocities
        figure(figQ)
        q = 20; % reducement factor for quiver plot
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
        
      % Calculate entrainment(-)/expansion(+) coefficient from plume velocities
        e_coeff = PUNV_mag/vel_char;
        entr = e_coeff(e_coeff<0);
        expn = e_coeff(e_coeff>0);
        
      % Calculate plume-averaged total coefficient, entrainment, and
      % expansion and their standard deviations at each timestep
        avg_coeff(t) = mean(e_coeff);
        std_coeff(t) = std(e_coeff);
        avg_entr(t) = mean(entr);
        std_entr(t) = std(entr);
        avg_expn(t) = mean(expn);
        std_expn(t) = std(expn);
        
        dlmwrite(fullfile(sprintf('%s',dir),sprintf('coeff_avg-std_%s.txt',run)),[avg_coeff(t) std_coeff(t)],'-append','delimiter','\n','precision','%0.6f');
        dlmwrite(fullfile(sprintf('%s',dir),sprintf('entr_avg-std_%s.txt',run)),[avg_entr(t) std_entr(t)],'-append','delimiter','\n','precision','%0.6f');
        dlmwrite(fullfile(sprintf('%s',dir),sprintf('expn_avg-std_%s.txt',run)),[avg_expn(t) std_expn(t)],'-append','delimiter','\n','precision','%0.6f');
        dlmwrite(fullfile(sprintf('%s',dir),sprintf('plot_time_%s.txt',run)),time(t)','-append','delimiter','\n','precision','%g');
        dlmwrite(fullfile(sprintf('%s',dir),sprintf('ecoeff_all_t%d.txt',time(t))),e_coeff,'delimiter','\n','precision','%0.6f');
        
  % Plot plume surface color-coded by entrainment coefficient
    nmap = 256;
%      colormap(jet(nmap));
    colormap([winter(round(nmap*-entrainment_cmin/(entrainment_cmax-entrainment_cmin)));flipud(autumn(nmap-round(nmap*-entrainment_cmin/(entrainment_cmax-entrainment_cmin))))]);
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
    title(tLEn,'FontSize',12,'FontWeight','bold');
    camlight('right')
    camlight('left')
    lighting gouraud

  % Calculate plume volume [m3] and height [m], and max radius of conic
  % plume
    numcells = sum(sum(sum(EPG <= plumeedge)));
    plumevolume(t) = numcells*XRES*YRES*ZRES;
    plumeheight(t) = max(plumeZ)*ZRES;
    dlmwrite(fullfile(sprintf('%s',dir),sprintf('plumevolume_%s.txt',run)),plumevolume(t),'-append','delimiter','\n','precision','%g');
    dlmwrite(fullfile(sprintf('%s',dir),sprintf('plumeheight_%s.txt',run)),plumeheight(t),'-append','delimiter','\n','precision','%g');
    plumetoprad(t) = sqrt(3*plumevolume(t)/(pi*plumeheight(t))) - VENT_R;
        
  % Calculate entrainment coefficient using Morton linear assumption and
  % using conic plume simplification
    e_Morton = PUNV_mag./plumeV';
    e_Mconic(t) = (plumetoprad(t) - plumetoprad(t-1))/(DT*vel_char);
    dlmwrite(fullfile(sprintf('%s',dir),sprintf('e_MortonConic_%s.txt',run)),e_Mconic(t),'-append','delimiter','\n','precision','%g');
    
    cd(dir)
      vidfigEn = 'EntrCurrent.jpg';
      saveas(figEn,vidfigEn);
      imgEn = imread(vidfigEn);
      writeVideo(vidEntr,imgEn);

      vidfigQ = 'QuiverCurrent.jpg';
      saveas(figQ,vidfigQ);
      imgQ = imread(vidfigQ);
      writeVideo(vidQ,imgQ);
      
      vidfigEP = 'EPGCurrent.jpg';
      saveas(figEP,vidfigEP);
      imgEP = imread(vidfigEP);
      writeVideo(vidEPG,imgEP);
            
      % Save each timestep as an individual figure in either a
      % multipage tif file or other image filetype (user-specified).
        if strcmp(imtype,'tif') == 1 || strcmp(imtype,'tiff') == 1
            imwrite(imgEn,sprintf('Entr_tsteps_%s.tif',run),'WriteMode','append')
            imwrite(imgEP,sprintf('EPG_tsteps_%s.tif',run),'WriteMode','append')
        else
            saveas(figEn,sprintf('Entr_t%03d_%s.%s',time(t),run,imtype));
            saveas(figEP,sprintf('EPG_t%03d_%s.%s',time(t),run,imtype));
        end
            
    end
    
    cd(dir)
    close(vidEPG);
    close(vidQ);
    close(vidEntr);

    
    % Plot total plume volume and change in plume volume over time
      if strcmp(PULSE,'T') == 1
        str = sprintf('%s: Unsteady flow %.1f Hz',titlerun,FREQ);
      elseif strcmp(PULSE,'F') == 1
        str = sprintf('%s: Steady flow',titlerun);
      end
      fig_plumevol = figure('Name','Plume Volume','visible',vis);
        hvol1 = subplot(2,1,1);
            plot(time,plumevolume)
            title(hvol1,{sprintf('%s: Total plume volume',str)},'FontWeight','bold','FontSize',10)
            xlabel(hvol1,{'Time (s)'},'FontWeight','bold','FontSize',10)
            ylabel(hvol1,{'Volume (m^3)'},'FontWeight','bold','FontSize',10)
        hvol2 = subplot(2,1,2);
            plot(time(2:end),diff(plumevolume))
            title(hvol2,{sprintf('%s: Change in plume volume',str)},'FontWeight','bold','FontSize',10)
            xlabel(hvol2,{'Time (s)'},'FontWeight','bold','FontSize',10)
            ylabel(hvol2,{'\DeltaVolume (m^3)'},'FontWeight','bold','FontSize',10)
      saveas(fig_plumevol,sprintf('PlumeVolume_%s.jpg',run));
      
    % Plot plume-averaged entrainment/expansion over time
      fig_coeff = figure('Name','Entrainment Coefficients','visible',vis);
        hold on
        errorbar(time,avg_coeff,std_coeff,'k')
        errorbar(time,avg_entr,std_entr,'b')
        errorbar(time,avg_expn,std_expn,'r')
        title(sprintf('%s: Plume-averaged coefficients',str),'FontWeight','bold','FontSize',10)
        xlabel('Time (s)','FontWeight','bold','FontSize',10)
        ylim([-1 1])
        legend({'Total coefficient','Entrainment','Expansion'},'Box','on','Location','EastOutside','FontWeight','bold','FontSize',10)
      saveas(fig_coeff,sprintf('Coefficients_%s.jpg',run));
      
    % Plot comparison Morton conic coefficient
      fig_Morton = figure('Name','Morton conic entrainment coefficient','visible',vis);
        plot(time,e_Mconic)
        title(sprintf('%s: Morton conic entrainment coefficient',str),'FontWeight','bold','FontSize',10)
        xlabel('Time (s)','FontWeight','bold','FontSize',10)
        ylabel('Coefficient','FontWeight','bold','FontSize',10)
        ylim([-1 1])
      saveas(fig_Morton,sprintf('MortonConic_%s.jpg',run));
            
    cd(postpath)
    
    sprintf('Entrainment processing complete.\nvidEntr_%s has been saved to %s.\nvidEPG_%s has been saved to %s',run,dir,run,dir)

end


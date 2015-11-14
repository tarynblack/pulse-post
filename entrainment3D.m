% function [ vidEntr ] = entrainment3D( run,dir,vis,ghostcells,IMAX,...
%     JMAX,KMAX,tickx,labelx,labelxunit,ticky,labely,labelyunit,tickz,...
%     labelz,labelzunit,plumeedge,XRES,YRES,ZRES,postpath,PULSE,FREQ,...
%     time,vel_char,cmin,cmax,viewaz,viewel,imtype )
%entrainment3D Summary of this function goes here
%   entrainment3D ---does things---
%
%   Special functions called: varchunk3D; pulsetitle
%   Last edit: Taryn Black, 13 November 2015

    varname = 'Entrainment';
    
%%% Initialize figure frames
    cd(dir)
    fig = figure('Name','Entrainment','visible',vis);
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
    
%%% Initialize video
    vidEntr = VideoWriter(sprintf('vidEntr_%d.avi',run));
    vidEntr.Quality = 100;
    vidEntr.FrameRate = 10;
    open(vidEntr);
    set(gcf,'Visible',vis);

%%% Calculate entrainment at each timestep, plot, and save video.
    fID_EPG = fopen(sprintf('EP_G_%d',run));
    fID_UG  = fopen(sprintf('U_G_%d',run));
    fID_VG  = fopen(sprintf('V_G_%d',run));
    fID_WG  = fopen(sprintf('W_G_%d',run));
        
    t = 0;
    while ~feof(fID_EPG)
            
        t = t+1;
        
        cd(postpath)
        try
            EPG = varchunk3D(fID_EPG,IMAX,JMAX,KMAX,ghostcells);
        catch ME
            warning('Error in varchunk3D at t=%d s:\n%s\nContinuing to next simulation.',time(t-1),ME.identifier)
            break
        end
        
        % Skip processing for first timestep - no isosurface at step one
          if t==1;
              continue
          end
        
        cla;

      % Find plume boundary (isosurface) and plot. 
        plumesurf = isosurface(EPG,plumeedge);
%         plotplume = patch(plumesurf,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.5);
%         hold on
        
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
        
      % Find velocities for entire domain at current timestep
        U_G = varchunk3D(fID_UG,IMAX,JMAX,KMAX,ghostcells);
        V_G = varchunk3D(fID_VG,IMAX,JMAX,KMAX,ghostcells);
        W_G = varchunk3D(fID_WG,IMAX,JMAX,KMAX,ghostcells);
        
      % Interpolate for velocities at plume surface vertices
        plumeU = interp3(U_G,plumeX,plumeY,plumeZ);
        plumeV = interp3(V_G,plumeX,plumeY,plumeZ);
        plumeW = interp3(W_G,plumeX,plumeY,plumeZ);
        plumevelocity = [plumeU plumeV plumeW]';
        
      % Calculate magnitudes and directional components of plume unit
      % normal velocities (PUNV)
        PUNV_mag = dot(plumevelocity,unitnorm); 
        PUNV_X = PUNV_mag.*unitnorm(1,:);
        PUNV_Y = PUNV_mag.*unitnorm(2,:);
        PUNV_Z = PUNV_mag.*unitnorm(3,:); 
        
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
        
        dlmwrite(fullfile(sprintf('%s',dir),'coeff_avg-std.txt'),[avg_coeff(t) std_coeff(t)],'-append','delimiter','\t','precision','%0.6f');
        dlmwrite(fullfile(sprintf('%s',dir),'entr_avg-std.txt'),[avg_entr(t) std_entr(t)],'-append','delimiter','\t','precision','%0.6f');
        dlmwrite(fullfile(sprintf('%s',dir),'expn_avg-std.txt'),[avg_expn(t) std_expn(t)],'-append','delimiter','\t','precision','%0.6f');
        dlmwrite(fullfile(sprintf('%s',dir),'plot_time.txt'),time(t)','-append','delimiter','\t','precision','%0.6f');
        
      % Plot plume surface color-coded by entrainment coefficient
        e_color = zeros(length(e_coeff),3);

    nmap = 256;
%      colormap(jet(nmap));
    colormap([winter(round(nmap*-cmin/(cmax-cmin)));flipud(autumn(nmap-round(nmap*-cmin/(cmax-cmin))))]);
    caxis([cmin cmax])
    cmap = colormap;
    emap = linspace(cmin,cmax,nmap);
    e_color = zeros(length(e_coeff),3);   
    
    for k = 1:length(e_coeff)
        if e_coeff(k) > cmax
            e_coeff(k) = cmax;
        elseif e_coeff(k) < cmin
            e_coeff(k) = cmin;
        end
    end

    e_round = interp1(emap,emap,e_coeff,'nearest');
    emap = [emap; 1:nmap];
    for j = 1:length(e_coeff)
        e_color(j,:) = cmap(emap(2,emap(1,:) == e_round(j)),:);
    end

    patch('Vertices',[plumeX plumeY plumeZ],'Faces',1:length(plumeX),...
        'FaceVertexCData',e_color,'FaceColor','none','EdgeColor',...
        'none','Marker','o','MarkerFaceColor','flat')
    colorbar

% patch('Vertices',[plumeX plumeY plumeZ],'Faces',1:length(plumeX),...
%             'FaceVertexCData',ecolor,'FaceColor','none','EdgeColor',...
%             'none','Marker','o','MarkerFaceColor','flat')
        
      % Calculate plume volume
        numcells = sum(sum(sum(EPG <= plumeedge)));
        plumevolume(t) = numcells*XRES*YRES*ZRES;
        
      % Calculate entrainment coefficient using Morton linear assumption
        e_Morton = PUNV_mag./plumeV';
        
        
        
        tL = pulsetitle(varname,PULSE,time,t,run,FREQ);
        title(tL,'FontSize',12,'FontWeight','bold');
        
        cd(dir)
          vidfig = 'EntrCurrent.jpg';
          saveas(fig,vidfig);
          img = imread(vidfig);
          writeVideo(vidEntr,img);
            
          % Save each timestep as an individual figure in either a
          % multipage tif file or other image filetype (user-specified).
            if strcmp(imtype,'tif') == 1 || strcmp(imtype,'tiff') == 1
                imwrite(img,sprintf('Entr_tsteps_%d.tif',run),'WriteMode','append')
            else saveas(fig,sprintf('Entr_%03ds_%d.%s',time(t),run,imtype));
            end
            
    end
    
    cd(dir)
    close(vidEntr);

    
    % Plot total plume volume and change in plume volume over time
      if strcmp(PULSE,'T') == 1
        str = sprintf('ID#%d: Unsteady flow %.1f Hz',run,FREQ);
      elseif strcmp(PULSE,'F') == 1
        str = sprintf('ID#%d: Steady flow',run);
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
      saveas(fig_plumevol,sprintf('PlumeVolume_%d.jpg',run));
      
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
      saveas(fig_coeff,sprintf('Coefficients_%d.jpg',run));
            
    cd(postpath)
    
    sprintf('Entrainment processing complete. \nvidEntr_%d has been saved to %s',run,dir)

% end


function [ vidPartConc ] = particleConc3D( run,dir,vis,IMAX,JMAX,KMAX,...
    ghostcells,tickx,labelx,labelxunit,ticky,labely,labelyunit,tickz,...
    labelz,labelzunit,XRES,YRES,ZRES,postpath,sdistX,sdistY,sdistZ,...
    particleConc_cmin,particleConc_cmax,titlerun,timesteps,PULSE,FREQ,...
    time,imtype,plumeedge )
%particleConc3D plots a volume slice of the concentration of each particle
%size over time.
%   Detailed explanation goes here
%   
%   Special functions called: varchunk3D; pulsetitle
%   Last edit: Taryn Black, 3 December 2015

  % Clear directory of appending files from previous processing attempts
    cd(dir)
    delete('PartConc_*');
    
    
  % ----------------------- FIGURE INITIALIZATION ----------------------- %
  % Define variable names for figures
    varS1 = 'S1 concentration';
    varS2 = 'S2 concentration';
    varS3 = 'S3 concentration';
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
    fig = figure('Name','Particle Concentrations','units','normalized',...
        'outerposition',[0 0 1 1],'visible',vis);
    hold on
    subfigS1 = subplot(1,3,1);
      hold on
      view(saz,sel)
      axis equal
      axis([ghostcells-1,IMAX-(ghostcells/2),ghostcells-1,...
          KMAX-(ghostcells/2),ghostcells-1,JMAX-(ghostcells/2)]);
      set(gca,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,'FontSize',12)
        xlabel(sprintf('\\bf Distance_x (%s)',labelxunit),'FontSize',12)
      set(gca,'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,'FontSize',12)
        ylabel(sprintf('\\bf Distance_z (%s)',labelzunit),'FontSize',12)
      set(gca,'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely,'FontSize',12)
        zlabel(sprintf('\\bf Altitude (%s)',labelyunit),'FontSize',12)
    subfigS2 = subplot(1,3,2);
      hold on
      view(saz,sel)
      axis equal
      axis([ghostcells-1,IMAX-(ghostcells/2),ghostcells-1,...
          KMAX-(ghostcells/2),ghostcells-1,JMAX-(ghostcells/2)]);
      set(gca,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,'FontSize',12)
        xlabel(sprintf('\\bf Distance_x (%s)',labelxunit),'FontSize',12)
      set(gca,'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,'FontSize',12)
        ylabel(sprintf('\\bf Distance_z (%s)',labelzunit),'FontSize',12)
      set(gca,'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely,'FontSize',12)
        zlabel(sprintf('\\bf Altitude (%s)',labelyunit),'FontSize',12)
     subfigS3 = subplot(1,3,3);
       hold on
       view(saz,sel)
       axis equal
       axis([ghostcells-1,IMAX-(ghostcells/2),ghostcells-1,...
           KMAX-(ghostcells/2),ghostcells-1,JMAX-(ghostcells/2)]);
       set(gca,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,'FontSize',12)
         xlabel(sprintf('\\bf Distance_x (%s)',labelxunit),'FontSize',12)
       set(gca,'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,'FontSize',12)
         ylabel(sprintf('\\bf Distance_z (%s)',labelzunit),'FontSize',12)
       set(gca,'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely,'FontSize',12)
         zlabel(sprintf('\\bf Altitude (%s)',labelyunit),'FontSize',12)
            
  % Initialize video
    vidPartConc = VideoWriter(sprintf('vidPartConc_%s.avi',run));
    vidPartConc.Quality = 100;
    vidPartConc.FrameRate = 10;
    open(vidPartConc);
    set(gcf,'Visible',vis);
  % ===================================================================== %
  
    
  % File import specifications: columns to read or skip for each variable
    EPS1import = '%*f%f%*f%*f%*f%*f%*f';
    EPS2import = '%*f%*f%f%*f%*f%*f%*f';
    EPS3import = '%*f%*f%*f%f%*f%*f%*f';
    EGimport   = '%f%*f%*f%*f%*f%*f%*f';
    
    
  % =================== B E G I N   T I M E   L O O P =================== %
    t = 0;
    while t <= timesteps 
        
        t = t+1;
        
      % Queue up current timestep files
        cd(dir)
        fclose('all');
        clear fID*;
        fID_EPS1 = fopen(sprintf('EP_t%02d.txt',t));
        fID_EPS2 = fopen(sprintf('EP_t%02d.txt',t));
        fID_EPS3 = fopen(sprintf('EP_t%02d.txt',t));
        fID_EPG  = fopen(sprintf('EP_t%02d.txt',t));
        cd(postpath)
      
      % Prepare particle vol. fractions for full domain at current timestep
        try
            EPS1 = varchunk3D(fID_EPS1,EPS1import,IMAX,JMAX,KMAX,ghostcells);
        catch ME
            warning('Error in varchunk3D at t=%d s:\n%s\nContinuing to next simulation.',...
                time(t),ME.identifier)
            break
        end
        EPS2 = varchunk3D(fID_EPS2,EPS2import,IMAX,JMAX,KMAX,ghostcells);
        EPS3 = varchunk3D(fID_EPS3,EPS3import,IMAX,JMAX,KMAX,ghostcells);
        EPG  = varchunk3D(fID_EPG,EPimport,IMAX,JMAX,KMAX,ghostcells);
        
      % Skip processing for first timestep when there is no plume.
        if t==1;
            continue
        end
        cla;
        
      % Calculate log of volume fractions
        logS1 = log10(EPS1 + 1E-10);
        logS2 = log10(EPS2 + 1E-10);
        logS3 = log10(EPS3 + 1E-10);
        
        
      % ------------- PARTICLE VOLUME FRACTION SLICE FIGURES ------------ %
                    % ----- WITH PLUME OUTLINE OVERLAY ----- %
        subplot(1,3,1)
          view(saz,sel)
          hS1 = slice(logS1,sdistX*IMAX,sdistY*KMAX,sdistZ*JMAX);
            hS1.FaceColor = 'interp';
            hS1.EdgeColor = 'none';
          tLS1 = pulsetitle(varS1,PULSE,time,t,titlerun,FREQ);
          title(tLS1,'FontSize',12,'FontWeight','bold'); 
          hEP1 = contourslice(EPG,sdistX*IMAX,sdistY*KMAX,0,[plumeedge plumeedge]);
          set(hEP1,'EdgeColor',[1 1 1],'LineWidth',0.5);
        subplot(1,3,2)
          view(saz,sel)
          hS2 = slice(logS2,sdistX*IMAX,sdistY*KMAX,sdistZ*JMAX);
            hS2.FaceColor = 'interp';
            hS2.EdgeColor = 'none';
          tLS2 = pulsetitle(varS2,PULSE,time,t,titlerun,FREQ);
          title(tLS2,'FontSize',12,'FontWeight','bold');
          hEP2 = contourslice(EPG,sdistX*IMAX,sdistY*KMAX,0,[plumeedge plumeedge]);
          set(hEP2,'EdgeColor',[1 1 1],'LineWidth',0.5);
        subplot(1,3,3)
          view(saz,sel)
          hS3 = slice(logS3,sdistX*IMAX,sdistY*KMAX,sdistZ*JMAX);
            hS3.FaceColor = 'interp';
            hS3.EdgeColor = 'none';
          tLS3 = pulsetitle(varS3,PULSE,time,t,titlerun,FREQ);
          title(tLS3,'FontSize',12,'FontWeight','bold');
          hEP3 = contourslice(EPG,sdistX*IMAX,sdistY*KMAX,0,[plumeedge plumeedge]);
          set(hEP3,'EdgeColor',[1 1 1],'LineWidth',0.5);
          hc = colorbar('location','eastoutside');
            caxis([particleConc_cmin particleConc_cmax]);
            ylabel(hc,'\bf log_1_0(Particle volume fraction)','FontSize',12)  
        PosS1 = get(subfigS1,'position');
        PosS2 = get(subfigS2,'position');
        PosS3 = get(subfigS3,'position');
        PosS2(3:4) = PosS1(3:4);
        PosS3(3:4) = PosS1(3:4);
        set(subfigS2,'position',PosS2);
        set(subfigS3,'position',PosS3);
      % ================================================================= %

      
      % --------- SAVE CURRENT FRAMES TO VIDEOS AND IMAGE FILES --------- %
        cd(dir)
        
      % Append current particle concentration frame to vidPartConc
        vidfig = 'PartConcCurrent.jpg';
        saveas(fig,vidfig);
        img = imread(vidfig);
        writeVideo(vidPartConc,img);
          
      % If user-specified image filetype is tif, append current timestep
      % frame to multipage tif file. Otherwise, save frame as independent
      % image named by timestep.
        if strcmp(imtype,'tif') == 1 || strcmp(imtype,'tiff') == 1
            imwrite(img,sprintf('PartConc_tsteps_%s.tif',run),'WriteMode','append')
        else
            saveas(fig,sprintf('PartConc_%03ds_%s.%s',time(t),run,imtype));
        end
      % ================================================================= %
          
    end
  % ===================== E N D   T I M E   L O O P ===================== %
    
  
  % End video write and finish video files
    cd(dir)
    close(vidPartConc);
    
    cd(postpath)    
    disp('Particle concentration processing complete.')
    fprintf('vidPartConc_%s has been saved to %s.\n',run,dir)

end

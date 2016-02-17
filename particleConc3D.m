function [ vidPartConc ] = particleConc3D( run,dir,vis,IMAX,JMAX,KMAX,...
    ghostcells,tickx,labelx,labelXunit,ticky,labely,labelYunit,tickz,...
    labelz,labelZunit,XRES,YRES,ZRES,postpath,sdistX,sdistY,sdistZ,...
    particleConc_cmin,particleConc_cmax,titlerun,timesteps,PULSE,FREQ,...
    time,imtype,plumeedge,D_S1,D_S2,D_S3 )
%particleConc3D plots a volume slice of the concentration of each particle
%size over time.
%   Detailed explanation goes here
%   
%   Special functions called: varchunk3D; pulsetitle
%   Last edit: Taryn Black, 12 February 2016

  % Clear directory of appending files from previous processing attempts
    cd(dir)
    delete('PartConc_*');
    
    
  % ----------------------- FIGURE INITIALIZATION ----------------------- %
  % Define variable names for figures
    varS1 = sprintf('Solid phase 1 (d = %g mm)',D_S1*1E3);
    varS2 = sprintf('Solid phase 2 (d = %g mm)',D_S2*1E3);
    varS3 = sprintf('Solid phase 3 (d = %g mm)',D_S3*1E3);
    varPC = 'Particle concentration';
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
    
  % Subtightplot properties
    gap = [0.01 0.01];
    ht = 0.15;
    wd = 0.15;
    
  % Figure and axes properties
    figPC = figure('Name','Particle Concentrations','units','pixels',...
        'outerposition',[1 44 1110 780],'visible',vis,'color','w',...
        'PaperPositionMode','auto');
    axS1 = subtightplot(1,3,1,gap,ht,wd);
      hold on
      axis(axS1,'equal',[0,IMAX-ghostcells,0,KMAX-ghostcells,0,...
            JMAX-ghostcells]);
      set(axS1,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,...
          'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,...
          'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely);
      xlabel(axS1,sprintf('\\bf Distance_x (%s)',labelXunit))
      ylabel(axS1,sprintf('\\bf Distance_z (%s)',labelZunit))
      zlabel(axS1,sprintf('\\bf Altitude (%s)',labelYunit))                
    axS2 = subtightplot(1,3,2,gap,ht,wd);
      hold on
      axis(axS2,'equal',[0,IMAX-ghostcells,0,KMAX-ghostcells,0,...
            JMAX-ghostcells]);
      set(axS2,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,...
          'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,...
          'ZTick',ticky(2:end)/YRES,'ZTickLabel','')
      xlabel(axS2,sprintf('\\bf Distance_x (%s)',labelXunit))
      ylabel(axS2,sprintf('\\bf Distance_z (%s)',labelZunit))
     axS3 = subtightplot(1,3,3,gap,ht,wd);
       hold on
       axis(axS3,'equal',[0,IMAX-ghostcells,0,KMAX-ghostcells,0,...
            JMAX-ghostcells]);
       set(axS3,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,...
          'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,...
          'ZTick',ticky(2:end)/YRES,'ZTickLabel','')
       xlabel(axS3,sprintf('\\bf Distance_x (%s)',labelXunit))
       ylabel(axS3,sprintf('\\bf Distance_z (%s)',labelZunit))
       cbPC = colorbar(axS3,'Location','eastoutside','AxisLocation',...
           'out','FontSize',12);
       cbPC.Label.String = '\bflog_{10}(Particle volume fraction)';
    set([axS1 axS2 axS3],'box','on','TickDir','in','FontSize',12)
    grid(axS1,'on'); grid(axS2,'on'); grid(axS3,'on');
    axS1.Layer = 'top'; axS2.Layer = 'top'; axS3.Layer = 'top';
    view(axS1,saz,sel);view(axS2,saz,sel);view(axS3,saz,sel);
            
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
        EPG  = varchunk3D(fID_EPG,EGimport,IMAX,JMAX,KMAX,ghostcells);
        
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
        colormap jet
        cla(axS1);cla(axS2);cla(axS3);
        hS1 = slice(axS1,0.5:(IMAX-ghostcells-0.5),0.5:(KMAX-ghostcells-0.5),...
            0.5:(JMAX-ghostcells-0.5),logS1,sdistX*(IMAX-ghostcells),...
            sdistY*(KMAX-ghostcells),sdistZ*(JMAX-ghostcells));
          hS1.FaceColor = 'interp';
          hS1.EdgeColor = 'none';
          tLS1 = pulsetitle(varPC,PULSE,time,t,titlerun,FREQ);
          title(axS1,sprintf('%s\n%s',tLS1{1},varS1));
          hEP1 = contourslice(axS1,EPG,sdistX*(IMAX-ghostcells),...
              sdistY*(KMAX-ghostcells),0,[plumeedge plumeedge]);
          set(hEP1,'EdgeColor',[1 1 1],'LineWidth',0.5);
          caxis(axS1,[particleConc_cmin,particleConc_cmax]);
        hS2 = slice(axS2,0.5:(IMAX-ghostcells-0.5),0.5:(KMAX-ghostcells-0.5),...
            0.5:(JMAX-ghostcells-0.5),logS2,sdistX*(IMAX-ghostcells),...
            sdistY*(KMAX-ghostcells),sdistZ*(JMAX-ghostcells));
          hS2.FaceColor = 'interp';
          hS2.EdgeColor = 'none';
          title(axS2,sprintf('%s',varS2));
          hEP2 = contourslice(axS2,EPG,sdistX*(IMAX-ghostcells),...
              sdistY*(KMAX-ghostcells),0,[plumeedge plumeedge]);
          set(hEP2,'EdgeColor',[1 1 1],'LineWidth',0.5);
          caxis(axS2,[particleConc_cmin,particleConc_cmax]);
        hS3 = slice(axS3,0.5:(IMAX-ghostcells-0.5),0.5:(KMAX-ghostcells-0.5),...
            0.5:(JMAX-ghostcells-0.5),logS3,sdistX*(IMAX-ghostcells),...
            sdistY*(KMAX-ghostcells),sdistZ*(JMAX-ghostcells));
          hS3.FaceColor = 'interp';
          hS3.EdgeColor = 'none';
          title(axS3,sprintf('%s\n%s',tLS1{2},varS3));
          hEP3 = contourslice(axS3,EPG,sdistX*(IMAX-ghostcells),...
              sdistY*(KMAX-ghostcells),0,[plumeedge plumeedge]);
          set(hEP3,'EdgeColor',[1 1 1],'LineWidth',0.5);
          caxis(axS3,[particleConc_cmin particleConc_cmax]);
        PosS1 = get(axS1,'position');
        PosS2 = get(axS2,'position');
        PosS3 = get(axS3,'position');
        PosS2(3:4) = PosS1(3:4);
        PosS3(3:4) = PosS1(3:4);
        set(axS2,'position',PosS2);
        set(axS3,'position',PosS3);
      % ================================================================= %

      
      % --------- SAVE CURRENT FRAMES TO VIDEOS AND IMAGE FILES --------- %
        cd(dir)
        
      % Append current particle concentration frame to vidPartConc
        vidfig = 'PartConcCurrent.jpg';
        saveas(figPC,vidfig);
        img = imread(vidfig);
        writeVideo(vidPartConc,img);
          
      % If user-specified image filetype is tif, append current timestep
      % frame to multipage tif file. Otherwise, save frame as independent
      % image named by timestep.
        if strcmp(imtype,'tif') == 1 || strcmp(imtype,'tiff') == 1
            imwrite(img,sprintf('PartConc_tsteps_%s.tif',run),'WriteMode','append')
        else
            saveas(figPC,sprintf('PartConc_%03ds_%s.%s',time(t),run,imtype));
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

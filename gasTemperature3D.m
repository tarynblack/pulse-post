function [ vidGasTemp ] = gasTemperature3D( run,dir,vis,ghostcells,IMAX,...
    JMAX,KMAX,tickx,labelx,labelxunit,ticky,labely,labelyunit,tickz,...
    labelz,labelzunit,XRES,YRES,ZRES,sdistX,sdistY,sdistZ,postpath,...
    ATMOS,TROPO,Y,BC_TG,gasTemperature_cmin,PULSE,FREQ,time,titlerun,...
    timesteps,imtype,plumeedge )
%gasTemperature3D plots a volume slice of the gas temperature of the plume
%over time.
%   Detailed explanation goes here
%   
%   Special functions called: varchunk3D; pulsetitle
%   Last edit: Taryn Black, 3 December 2015

  % Clear directory of appending files from previous processing attempts
    cd(dir)
    delete('GasTemp_*');
    
    
  % ----------------------- FIGURE INITIALIZATION ----------------------- %
  % Define variable names for figures
    varTG = 'Gas temperature';
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
    fig = figure('Name','Gas Temperature','visible',vis);
    hold on
    view(saz,sel)
    axis equal
    axis([ghostcells-1,IMAX-(ghostcells/2),ghostcells-1,...
        KMAX-(ghostcells/2),ghostcells-1,JMAX-(ghostcells/2)]);
    set(gca,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,'FontSize',12)
        xlabel(sprintf('\\bf Distance (%s)',labelxunit),'FontSize',12)
    set(gca,'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,'FontSize',12)
        ylabel(sprintf('\\bf Distance (%s)',labelzunit),'FontSize',12)
    set(gca,'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely,'FontSize',12)
        zlabel(sprintf('\\bf Altitude (%s)',labelyunit),'FontSize',12)
    
  % Initialize video
    vidGasTemp = VideoWriter(sprintf('vidGasTemp_%s.avi',run));
    vidGasTemp.Quality = 100;
    vidGasTemp.FrameRate = 10;
    open(vidGasTemp);
    set(gcf,'Visible',vis);
  % ===================================================================== %

  
  % File import specifications: columns to read or skip for each variable
    TGimport = '%f%*f%*f%*f';
    EGimport = '%f%*f%*f%*f%*f%*f%*f';

    
  % =================== B E G I N   T I M E   L O O P =================== %
    t = 0;
    while t <= timesteps
        
        t = t+1;
        
      % Queue up current timestep files
        cd(dir)
        fclose('all');
        clear fID*;
        fID_TG = fopen(sprintf('T_G_t%02d.txt',t));
        fID_EPG = fopen(sprintf('EP_t%02d.txt',t));
        cd(postpath)
        
      % Prepare gas temperature for full domain at current timestep
        try
            TG = varchunk3D(fID_TG,TGimport,IMAX,JMAX,KMAX,ghostcells);
        catch ME
            warning('Error in varchunk3D at t=%d s:\n%s\nContinuing to next simulation.',...
                time(t),ME.identifier)
            break
        end
        EPG = varchunk3D(fID_EPG,EGimport,IMAX,JMAX,KMAX,ghostcells);
        
      % Skip processing for first timestep when there is no plume.
        if t==1;
            continue
        end
        cla;
        
      % Apply atmospheric correction (equation of state for gas)
        if strcmp(ATMOS,'T') == 1
            for i = 1:(JMAX-ghostcells)
                if Y(i) <= TROPO
                    TG(:,:,i) = TG(:,:,i) - 0.0098*Y(i);
                elseif Y(i) > TROPO
                    TG(:,:,i) = TG(:,:,i) - 0.0098*TROPO + 0.001*Y(i);
                end
            end
        end
        
        
      % ------------------- TEMPERATURE SLICE FIGURE -------------------- %
        hTG = slice(TG,sdistX*IMAX,sdistY*KMAX,sdistZ*JMAX);
          hTG.FaceColor = 'interp';
          hTG.EdgeColor = 'none';
        hc = colorbar;
          caxis([gasTemperature_cmin BC_TG]);
          ylabel(hc,'\bf Temperature [K]','FontSize',12)
        tL = pulsetitle(varTG,PULSE,time,t,titlerun,FREQ);
        title(tL,'FontSize',12,'FontWeight','bold');
      % ================================================================= %
      
      
      % --------------------- OVERLAY PLUME OUTLINE --------------------- %
        hEP = contourslice(EPG,sdistX*IMAX,sdistY*KMAX,0,[plumeedge plumeedge]);
        set(hEP,'EdgeColor',[1 1 1],'LineWidth',0.5);
      % ================================================================= %
        
      
      % --------- SAVE CURRENT FRAMES TO VIDEOS AND IMAGE FILES --------- %
        cd(dir)
        
      % Append current gas temperature frame to vidGasTemp
        vidfig = 'GasTempCurrent.jpg';
        saveas(fig,vidfig);
        img = imread(vidfig);
        writeVideo(vidGasTemp,img);
            
      % If user-specified image filetype is tif, append current timestep
      % frame to multipage tif file. Otherwise, save frame as independent
      % image named by timestep.
        if strcmp(imtype,'tif') == 1 || strcmp(imtype,'tiff') == 1
            imwrite(img,sprintf('GasTemp_tsteps_%s.tif',run),'WriteMode','append')
        else
            saveas(fig,sprintf('GasTemp_%03ds_%s.%s',time(t),run,imtype));
        end
      % ================================================================= %
                    
    end
  % ===================== E N D   T I M E   L O O P ===================== %
    
  
  % End video write and finish video files
    cd(dir)
    close(vidGasTemp)
    
    cd(postpath)
    disp('Gas temperature processing complete.')
    fprintf('vidGasTemp_%s has been saved to %s.\n',run,dir) 

end


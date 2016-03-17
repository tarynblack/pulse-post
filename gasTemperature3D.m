function [ vidGasTemp ] = gasTemperature3D( run,runpath,vis,ghostcells,IMAX,...
    JMAX,KMAX,tickx,labelx,labelXunit,ticky,labely,labelYunit,tickz,...
    labelz,labelZunit,XRES,YRES,ZRES,sdistX,sdistY,sdistZ,postpath,...
    ATMOS,TROPO,YGRID,BC_TG,gasTemperature_cmin,PULSE,FREQ,time,titlerun,...
    timesteps,imtype,plumeedge,savepath,readTG,fnameTG,readEPG,fnameEPG )
%gasTemperature3D plots a volume slice of the gas temperature of the plume
%over time.
%   Detailed explanation goes here
%   
%   Special functions called: loadTimestep3D; pulsetitle
%   Last edit: Taryn Black, 17 March 2016

  % Clear directory of appending files from previous processing attempts
    cd(savepath)
    delete('GasTemp_*');
    
    
  % ----------------------- FIGURE INITIALIZATION ----------------------- %
  % Define variable names for figures
    varTG = 'Gas temperature';
    cd(runpath)
    
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
        sel = 90;
    elseif isempty(sdistY) && isempty(sdistZ)
        saz = 90;
        sel = 0;
    elseif isempty(sdistX) && isempty(sdistZ)
        saz = 0;
        sel = 0;
    else [saz,sel] = view(3);
    end
    
  % Figure and axes properties
    figTemp = figure('Name','Gas Temperature','visible',vis,'units',...
        'normalized','outerposition',[0.5 0 0.45 1],'PaperPositionMode',...
        'auto','color','w');
    axTemp = axes('Parent',figTemp,'box','on','TickDir','in','FontSize',12);
    hold on
    grid(axTemp,'on');axTemp.Layer = 'top';
    view(axTemp,saz,sel)
    axis(axTemp,'equal',[0,IMAX-ghostcells,0,KMAX-ghostcells,0,...
        JMAX-ghostcells]);
    set(axTemp,'XTick',tickx(2:end)/XRES,'XTickLabel',labelx,...
        'YTick',tickz(2:end)/ZRES,'YTickLabel',labelz,...
        'ZTick',ticky(2:end)/YRES,'ZTickLabel',labely);
    xlabel(axTemp,sprintf('\\bf Distance_x (%s)',labelXunit))
    ylabel(axTemp,sprintf('\\bf Distance_z (%s)',labelZunit))
    zlabel(axTemp,sprintf('\\bf Altitude (%s)',labelYunit))
    cbTemp = colorbar(axTemp,'AxisLocation','in','FontSize',12);
    cbTemp.Label.String = '\bfTemperature (K)';
        
    cvalsPC = log10(-plumeedge+1):1:-2;
    cmapPC = colormap(flipud(bone(length(cvalsPC))));
    colormap(figTemp,'default');
    
  % Initialize video
    cd(savepath)
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
        cd(runpath)
        fclose('all');
        clear fID*;
        
        cd(postpath)
        fID_TG  = fileReadType(fnameTG,readTG,t,runpath,postpath);
        fID_EPG = fileReadType(fnameEPG,readEPG,t,runpath,postpath);
        
      % Prepare gas temperature for full domain at current timestep
        try
            TG = loadTimestep3D(fID_TG,TGimport,readTG,IMAX,JMAX,KMAX,ghostcells);
        catch ME
            warning('Error in loadTimestep3D at t=%d s:\n%s\nContinuing to next simulation.',...
                time(t),ME.identifier)
            break
        end
        EPG = loadTimestep3D(fID_EPG,EGimport,readEPG,IMAX,JMAX,KMAX,ghostcells);
        
      % Skip processing for first timestep when there is no plume.
        if t==1;
            continue
        end
        cla;
        
      % Apply atmospheric correction (equation of state for gas)
        if strcmp(ATMOS,'T') == 1
            for i = 1:(JMAX-ghostcells)
                if YGRID(i) <= TROPO
                    TG(:,:,i) = TG(:,:,i) - 0.0098*YGRID(i);
                elseif YGRID(i) > TROPO
                    TG(:,:,i) = TG(:,:,i) - 0.0098*TROPO + 0.001*YGRID(i);
                end
            end
        end
        
      % Calculate log of particle volume fraction
        logParticles = log10(1 - EPG + 1E-10);
        
        
      % ------------------- TEMPERATURE SLICE FIGURE -------------------- %
        figure(figTemp)
        cla(axTemp);
        hTG = slice(0.5:(IMAX-ghostcells-0.5),0.5:(KMAX-ghostcells-0.5),...
            0.5:(JMAX-ghostcells-0.5),TG,sdistX*(IMAX-ghostcells),...
            sdistY*(KMAX-ghostcells),sdistZ*(JMAX-ghostcells));
        hTG.FaceColor = 'interp';
        hTG.EdgeColor = 'none';
        caxis(axTemp,[gasTemperature_cmin BC_TG]);
        tL = pulsetitle(varTG,PULSE,time,t,titlerun,FREQ);
        title(tL,'FontSize',12,'FontWeight','bold');
        hold on
      % ================================================================= %
      
      
      % --------------------- OVERLAY PLUME OUTLINE --------------------- %
%         hEP = contourslice(EPG,sdistX*IMAX,sdistY*KMAX,0,[plumeedge plumeedge]);
%         set(hEP,'EdgeColor',[1 1 1],'LineWidth',0.5);
      % ================================================================= %
      
      
      % --------- OVERLAY LOG PARTICLE VOLUME FRACTION CONTOURS --------- %
        for j = 1:length(cvalsPC)
            hPS = contourslice(logParticles,sdistX*(IMAX-ghostcells),...
                sdistY*(KMAX-ghostcells),0,[cvalsPC(j) cvalsPC(j)]);
            set(hPS,'EdgeColor',cmapPC(j,:),'LineWidth',1.5);
            legnam{j} = sprintf('10^{%g}',cvalsPC(j));
            leg(j) = scatter(0,0,'s','filled','MarkerFaceColor',cmapPC(j,:));
            if leg(j).MarkerFaceColor == [1 1 1]
                set(leg(j),'MarkerEdgeColor',[0.8 0.8 0.8])
            end
            set(leg(j),'visible','off');
        end       
        hLeg = legend(leg,char(legnam));
        set(hLeg,'Box','off','Location','eastoutside');
        hLT = legendTitle(hLeg,sprintf('Particle\nConcentration\nContours'));
      % ================================================================= %
        
      
      % --------- SAVE CURRENT FRAMES TO VIDEOS AND IMAGE FILES --------- %
        cd(savepath)
        
      % Append current gas temperature frame to vidGasTemp
        vidfig = 'GasTempCurrent.jpg';
        saveas(figTemp,fullfile(savepath,vidfig));
        img = imread(vidfig);
        writeVideo(vidGasTemp,img);
            
      % If user-specified image filetype is tif, append current timestep
      % frame to multipage tif file. Otherwise, save frame as independent
      % image named by timestep.
        if strcmp(imtype,'tif') == 1 || strcmp(imtype,'tiff') == 1
            imwrite(img,fullfile(savepath,sprintf('GasTemp_tsteps_%s.tif',...
                run)),'tif','WriteMode','append')
        else
            saveas(figTemp,fullfile(savepath,...
                sprintf('GasTemp_%03ds_%s.%s',time(t),run,imtype)));
        end
      % ================================================================= %
                    
    end
  % ===================== E N D   T I M E   L O O P ===================== %
    
  
  % End video write and finish video files
    cd(savepath)
    close(vidGasTemp)
    
    cd(postpath)
    disp('Gas temperature processing complete.')
    fprintf('vidGasTemp_%s has been saved to %s.\n',run,savepath) 

end

function [ vidGasTemp ] = gasTemperature3D( run,dir,vis,ghostcells,IMAX,...
    JMAX,KMAX,tickx,labelx,labelxunit,ticky,labely,labelyunit,tickz,...
    labelz,labelzunit,XRES,YRES,ZRES,sdistX,sdistY,sdistZ,postpath,...
    ATMOS,TROPO,Y,BC_TG,gasTemperature_cmin,PULSE,FREQ,time,titlerun,...
    timesteps,imtype )
%gasTemperature3D plots a volume slice of the gas temperature of the plume
%over time.
%   Detailed explanation goes here
%   
%   Special functions called: varchunk3D; pulsetitle
%   Last edit: Taryn Black, 17 November 2015

    varname = 'Gas temperature';
    
%%% Ensure that 'no slice' directions are empty and determine figure
%%% viewing angle based on slice direction
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
    
%%% Initialize figure frames
    cd(dir)
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
    
%%% Initialize video
    vidGasTemp = VideoWriter(sprintf('vidGasTemp_%s.avi',run));
    vidGasTemp.Quality = 100;
    vidGasTemp.FrameRate = 10;
    open(vidGasTemp);
    set(gcf,'Visible',vis);
    
%%% Plot gas temperature slice at each timestep and save video.
%     fID_TG = fopen(sprintf('T_G_%s',run));
     TGimport = '%f%*f%*f%*f';
    
    t = 0;
    while t <= timesteps %~feof(fID_TG)
        
        t = t+1;
        
        fID_TG = fopen(sprintf('T_G_t%02d.txt',t));
        
        cd(postpath)
        try
            TG = varchunk3D(fID_TG,TGimport,IMAX,JMAX,KMAX,ghostcells);
        catch ME
            warning('Error in varchunk3D at t=%d s:\n%s\nContinuing to next simulation.',time(t),ME.identifier)
            break
        end
        
      % Skip processing for first timestep when there is no plume.
        if t==1;
            continue
        end
        cla;
        
      % Apply atmospheric correction if necessary ~ equation of state for
      % gas, above/below tropopause.
        if strcmp(ATMOS,'T') == 1
            for i = 1:(JMAX-ghostcells)
                if Y(i) <= TROPO
                    TG(:,:,i) = TG(:,:,i) - 0.0098*Y(i);
                elseif Y(i) > TROPO
                    TG(:,:,i) = TG(:,:,i) - 0.0098*TROPO + 0.001*Y(i);
                end
            end
        end
        
        hTG = slice(TG,sdistX*IMAX,sdistY*KMAX,sdistZ*JMAX);
            hTG.FaceColor = 'interp';
            hTG.EdgeColor = 'none';
        hc = colorbar;
            caxis([gasTemperature_cmin BC_TG]);
            ylabel(hc,'\bf Temperature [K]','FontSize',12)
        tL = pulsetitle(varname,PULSE,time,t,titlerun,FREQ);
        title(tL,'FontSize',12,'FontWeight','bold');
        
        cd(dir)
            vidfig = 'GasTempCurrent.jpg';
            saveas(fig,vidfig);
            img = imread(vidfig);
            writeVideo(vidGasTemp,img);
            
      % Save each timestep as an individual figure in either a
      % multipage tif file or other image filetype (user-specified).
        if strcmp(imtype,'tif') == 1 || strcmp(imtype,'tiff') == 1
            imwrite(img,sprintf('GasTemp_tsteps_%s.tif',run),'WriteMode','append')
        else saveas(fig,sprintf('GasTemp_%03ds_%s.%s',time(t),run,imtype));
        end
        
        fclose(fID_TG);
            
    end
    
    cd(dir)
    close(vidGasTemp)
    cd(postpath)
    sprintf('Gas temperature processing complete. \nvidGasTemp_%s has been saved to %s',run,dir) 

end


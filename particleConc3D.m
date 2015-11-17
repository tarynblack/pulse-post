function [ vidPartConc ] = particleConc3D( run,dir,vis,IMAX,JMAX,KMAX,...
    ghostcells,tickx,labelx,labelxunit,ticky,labely,labelyunit,tickz,...
    labelz,labelzunit,XRES,YRES,ZRES,postpath,sdistX,sdistY,sdistZ,...
    RO_S1,RO_S2,RO_S3,particleConc_cmin,particleConc_cmax )
%particleConc3D plots a volume slice of the concentration of each particle
%size over time.
%   Detailed explanation goes here
%   
%   Special functions called: varchunk3D; pulsetitle
%   Last edit: Taryn Black, 17 November 2015

    varS1 = 'S1 concentration';
    varS2 = 'S2 concentration';
    varS3 = 'S3 concentration';
    
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
    fig = figure('Name','Particle Concentrations','units','normalized','outerposition',[0 0 1 1],'visible',vis);
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
            
%%% Initialize video
    vidPartConc = VideoWriter(sprintf('vidPartConc_%s.avi',run));
    vidPartConc.Quality = 100;
    vidPartConc.FrameRate = 10;
    open(vidPartConc);
    set(gcf,'Visible',vis);
    
%%% Plot particle volume fraction concentrations at each timestep and save
%%% video.
    fID_S1 = fopen(sprintf('ROP_S1_%s',run));
    fID_S2 = fopen(sprintf('ROP_S2_%s',run));
    fID_S3 = fopen(sprintf('ROP_S3_%s',run));

    t = 0;
    while ~feof(fID_EPG)
        
        t = t+1;
        
        cd(postpath)
        try
            ROPS1 = varchunk3D(fID_S1,IMAX,JMAX,KMAX,ghostcells);
        catch ME
            warning('Error in varchunk3D at t=%d s:\n%s\nContinuing to next simulation.',time(t-1),ME.identifier)
            break
        end
        ROPS2 = varchunk3D(fID_S2,IMAX,JMAX,KMAX,ghostcells);
        ROPS3 = varchunk3D(fID_S3,IMAX,JMAX,KMAX,ghostcells);
        
      % Skip processing for first timestep - no isosurface at step one
        if t==1;
            continue
        end
        
        cla;
        
%%% #TODO# ROP*/RO* = EPS*. After switching to using Mary's output, change
%%% this calculation - she writes out EPS* rather than ROP*.
        logS1 = log10((ROPS1/RO_S1) + 1E-10);
        logS2 = log10((ROPS2/RO_S2) + 1E-10);
        logS3 = log10((ROPS3/RO_S3) + 1E-10);
        
        subplot(1,3,1)
            view(saz,sel)
            hS1 = slice(logS1,sdistX*IMAX,sdistY*KMAX,sdistZ*JMAX);
            hS1.FaceColor = 'interp';
            hS1.EdgeColor = 'none';
            tLS1 = pulsetitle(varS1,PULSE,time,t,run,FREQ);
            title(tLS1,'FontSize',12,'FontWeight','bold'); 
        subplot(1,3,2)
            view(saz,sel)
            hS2 = slice(logS2,sdistX*IMAX,sdistY*KMAX,sdistZ*JMAX);
            hS2.FaceColor = 'interp';
            hS2.EdgeColor = 'none';
            tLS2 = pulsetitle(varS2,PULSE,time,t,run,FREQ);
            title(tLS2,'FontSize',12,'FontWeight','bold');
        subplot(1,3,3)
            view(saz,sel)
            hS3 = slice(logS3,sdistX*IMAX,sdistY*KMAX,sdistZ*JMAX);
            hS3.FaceColor = 'interp';
            hS3.EdgeColor = 'none';
            tLS3 = pulsetitle(varS3,PULSE,time,t,run,FREQ);
            title(tLS3,'FontSize',12,'FontWeight','bold');
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

        cd(dir)
          vidfig = 'PartConcCurrent.jpg';
          saveas(fig,vidfig);
          img = imread(vidfig);
          writeVideo(vidPartConc,img);
          
    end
    
    cd(dir)
    close(vidPartConc)
    
    cd(postpath)
    
    sprintf('Particle concentration processing complete. \nvidPartConc_%s has been saved to %s',run,dir)

end

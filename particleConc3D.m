% function [ vidPartConc ] = particleConc3D( run,dir,vis,IMAX,JMAX,KMAX,...
%     ghostcells,tickx,labelx,labelxunit,ticky,labely,labelyunit,tickz,...
%     labelz,labelzunit,XRES,YRES,ZRES,postpath,X,Y,Z )
%particleConc3D plots a volume slice of the concentration of each particle
%size over time.
%   Detailed explanation goes here
%   
%   Special functions called: varchunk3D; pulsetitle
%   Last edit: Taryn Black, 14 November 2015

    varS1 = 'S1 concentration';
    varS2 = 'S2 concentration';
    varS3 = 'S3 concentration';
    
%%% Initialize figure frames
    cd(dir)
    fig = figure('Name','Particle Concentrations','units','normalized','outerposition',[0 0 1 1],'visible',vis);
    hold on
        subfigS1 = subplot(1,3,1);
            hold on
            view(90,0)
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
        subfigS2 = subplot(1,3,2);
            hold on
            view(90,0)
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
         subfigS3 = subplot(1,3,3);
            hold on
            view(90,0)
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
    vidPartConc = VideoWriter(sprintf('vidPartConc_%d.avi',run));
    vidPartConc.Quality = 100;
    vidPartConc.FrameRate = 10;
    open(vidPartConc);
    set(gcf,'Visible',vis);
    
%%% Plot particle volume fraction concentrations at each timestep and save
%%% video.
%%% #TODO# change from gas to particles when on cluster!
    fID_EPG = fopen(sprintf('EP_G_%d',run));
%     fID_S1 = fopen(sprintf('ROP_S1_%d',run));
%     fID_S2 = fopen(sprintf('ROP_S2_%d',run));
%     fID_S3 = fopen(sprintf('ROP_S3_%d',run));
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
        
        logEPG = log10(1-EPG+1E-10);
        
        subplot(1,3,1)
            hS1 = slice(logEPG,50,[],[]);
            hS1.FaceColor = 'interp';
            hS1.EdgeColor = 'none';
            tLS1 = pulsetitle(varS1,PULSE,time,t,run,FREQ);
            title(tLS1,'FontSize',12,'FontWeight','bold'); 
        subplot(1,3,2)
            hS2 = slice(logEPG,50,[],[]);
            hS2.FaceColor = 'interp';
            hS2.EdgeColor = 'none';
            tLS2 = pulsetitle(varS2,PULSE,time,t,run,FREQ);
            title(tLS2,'FontSize',12,'FontWeight','bold');
        subplot(1,3,3)
            hS3 = slice(logEPG,50,[],[]);
            hS3.FaceColor = 'interp';
            hS3.EdgeColor = 'none';
            tLS3 = pulsetitle(varS3,PULSE,time,t,run,FREQ);
            title(tLS3,'FontSize',12,'FontWeight','bold');
            hc = colorbar('location','eastoutside');
            
        PosS1 = get(subfigS1,'position');
        PosS2 = get(subfigS2,'position');
        PosS3 = get(subfigS3,'position');
        PosS2(3:4) = [PosS1(3:4)];
        PosS3(3:4) = [PosS1(3:4)];
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

% end

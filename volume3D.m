function [ vidEPG ] = volume3D( run,dir,vis,ghostcells,tickx,labelx,labelxunit,...
                ticky,labely,labelyunit,tickz,labelz,labelzunit,...
                plumeedge,XRES,YRES,ZRES,timesteps,IMAX,JMAX,KMAX,isoEPG,colEPG,trnEPG,time,...
                PULSE,FREQ,postpath )
%volume3D makes frames of the time evolution of gas volume fraction in a
%volcanic plume during a simulated volcanic eruption.
%   volume3D processes gas volume fraction (EP_G) data from an MFiX run
%   <run> to create video frames of the evolution of gas volume fraction
%   during an eruption, which are passed as output to a file
%   vidEPG_<run>.avi. The frames display specified isosurfaces <isoEPG> of
%   gas volume fraction, at specified colors <colEPG> and transparencies
%   <trnEPG>.
%   
%   Special functions called: varchunk3D; pulsetitle
%   Last edit: Taryn Black, 12 November 2015

    
    varname = 'Gas volume fraction';
    
%%% Initialize figure frames
    cd(dir)
    fig = figure('Name','Gas Volume Fraction','visible',vis);
    hold on
    view(3)
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
    %%% Initialize legend entries based on number of isosurfaces
%         names = cell(1,length(isoEPG));
%         cd(postpath)
%         initEPG = timeslice3D(EP_G,1,IMAX,JMAX,KMAX,ghostcells);
%         cd(dir)
%         for i = 1:length(isoEPG)
%             names{i} = sprintf('EPG = %0.5f',isoEPG(i));
%             surf = patch(isosurface(initEPG,isoEPG(i)));
%             set(surf,'FaceColor',colEPG(i,:),'EdgeColor','none','FaceAlpha',1.0);
%             hold on
%         end
%         legend(names)
    
%%% Initialize video
    vidEPG = VideoWriter(sprintf('vidEPG_%d.avi',run));
    vidEPG.Quality = 100;
    vidEPG.FrameRate = 10;
    open(vidEPG);
    set(gcf,'Visible',vis);
%     cd(postpath)
    
%%% Plot gas volume fraction isosurfaces at each timestep and save video.
%     cd(dir)
    fileID = fopen(sprintf('EP_G_%d',run));
    t = 0;
    while ~feof(fileID)
        
        t = t+1;
        
        cd(postpath)
        try
            EPG = varchunk3D(fileID,IMAX,JMAX,KMAX,ghostcells);
        catch ME
            warning('Error in varchunk3D at t=%d s:\n%s\nContinuing to next simulation.',time(t-1),ME.identifier)
            break
        end
%         cd(dir)

%     for t = 1:timesteps
%         
%         EPG = timeslice3D(EP_G,t,IMAX,JMAX,KMAX,ghostcells);
%         
        cla;
        for j = 1:length(isoEPG)
            surf(j) = patch(isosurface(EPG,isoEPG(j)));
            set(surf(j),'FaceColor',colEPG(j,:),'EdgeColor','none','FaceAlpha',trnEPG(j));
            hold on
        end
        
        tL = pulsetitle(varname,PULSE,time,t,run,FREQ);
        title(tL,'FontSize',12,'FontWeight','bold');
  
	cd(dir)      
        vidfig = 'EPGcurrent.jpg';
        saveas(fig,vidfig);
        img = imread(vidfig);
        writeVideo(vidEPG,img);
        
    end
    
%     close(fileID)
    cd(dir)
    close(vidEPG);
    cd(postpath)
    
    sprintf('EPG processing complete.\nvidEPG_%d has been saved to %s',run,dir)

end


function [mGV] = volume2D(EP_G,timesteps,IMAX,JMAX,X,Y,distance,altitude,run,PULSE,FREQ,time)
% volume2D makes frames of the time evolution of gas volume fraction during
% a simulated volcanic eruption.
%   This script processes gas volume fraction (EP_G) data from a specified
%   MFIX run (variable <run>) to create movie frames of the evolution of
%   gas volume fraction during an eruption, which are passed as output for
%   a file gVol_<run>.avi. The frames display log(EP_G) for better data
%   resolution at high EP_G, however the colorbar is set to display actual
%   values of EP_G rather than log(EP_G). 
%   Taryn Black, last edit 24 April 2015

    cd(sprintf('%d',run))
    
    contours = [0 0.9 0.95 0.99 0.995 0.999 0.9995 0.9999 0.99995 0.99999 1];
    log_contours = log10(1 - contours + 1E-10);
    labels = fliplr(contours);
    
    mGV = moviein(timesteps);
        hold on
        axis equal
        set(gca,'XTick',distance(2:end)*1000,'XTickLabel',distance(2:end),'FontSize',12)
            xlabel('\bf Distance (km)','FontSize',12)
        set(gca,'YTick',altitude(2:end)*1000,'YTickLabel',altitude(2:end),'FontSize',12)
            ylabel('\bf Altitude (km)','FontSize',12)
        hc = colorbar;
            set(hc,'YTickLabel',labels,'FontSize',12)
            ylabel(hc,'\bf Gas volume fraction','FontSize',12)

    for t = 1:timesteps
        
        EPG = reshape(EP_G((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);
%         EPG = EPG(2:end-1,2:end-1);
        EPGlog = log10(1 - EPG + 1E-10);        
        contourf(X,Y,EPGlog,log_contours);
        axis tight

        if strcmp(PULSE,'T') == 1
            str = 'Unsteady flow';
            tL = ({sprintf('ID#%d: %s, %.1f Hz',run,str,FREQ);sprintf('Gas volume fraction, t=%d s',time(t))});
        elseif strcmp(PULSE,'F') == 1
            str = 'Steady flow';
            tL = ({sprintf('ID#%d: %s',run,str);sprintf('Gas volume fraction, t=%d s',time(t))});
        end
        title(tL,'FontSize',12,'FontWeight','bold');
            
        mGV(:,t) = getframe(gcf);
        
        if time(t) == 300;
            saveas(gcf,sprintf('gVol300s_%d.jpg',run))
        end
        if time(t) == 450;
            saveas(gcf,sprintf('gVol450s_%d.jpg',run))
        end
        if t == timesteps;
            saveas(gcf,sprintf('gVolEND_%d.jpg',run))
        end
    end

% cd ..
    
end


function [mGP] = pressure2D(P_G,timesteps,IMAX,JMAX,X,Y,distance,altitude,run,P_range,PULSE,FREQ,time)
% pressure2D makes frames of the time evolution of gas pressure during
% a simulated volcanic eruption.
%   This script processes gas pressure (P_G) data from a specified MFIX
%   run (variable <run>) to create movie frames of the evolution of gas
%   pressure during an eruption, which are passed as output for a file
%   gPres_<run>.avi.
%   Taryn Black, last edit 24 April 2015

    cd(sprintf('%d',run))

    mGP = moviein(timesteps);
        hold on
        axis equal
        set(gca,'XTick',distance(2:end)*1000,'XTickLabel',distance(2:end),'FontSize',12)
            xlabel('\bf Distance (km)','FontSize',12)
        set(gca,'YTick',altitude(2:end)*1000,'YTickLabel',altitude(2:end),'FontSize',12)
            ylabel('\bf Altitude (km)','FontSize',12)
        hc = colorbar;
            ylabel(hc,'\bf Gas pressure (Pa)','FontSize',12)
            set(hc,'FontSize',12)
            caxis(P_range)

    for t = 1:timesteps
        
        PG = reshape(P_G((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);
        PG = PG(2:end-1,2:end-1);
        contourf(X,Y,PG)
        axis tight
        
        if strcmp(PULSE,'T') == 1
            str = 'Unsteady flow';
            tL = ({sprintf('ID#%d: %s, %.1f Hz',run,str,FREQ);sprintf('Pressure, t=%d s',time(t))});
        elseif strcmp(PULSE,'F') == 1
            str = 'Steady flow';
            tL = ({sprintf('ID#%d: %s',run,str);sprintf('Pressure, t=%d s',time(t))});
        end
        title(tL,'FontSize',12,'FontWeight','bold');
            
        mGP(:,t)=getframe(gcf);
        
        if t == timesteps;
            saveas(gcf,sprintf('gPresEND_%d.jpg',run))
        end
    end

% cd ..
    
end


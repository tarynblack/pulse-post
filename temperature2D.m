 function [mGT] = temperature2D(T_G,timesteps,IMAX,JMAX,X,Y,distance,altitude,run,T_range,PULSE,FREQ,time)
% temperature2D makes frames of the time evolution of gas temperature during
% a simulated volcanic eruption.
%   This script processes gas temperature (T_G) data from a specified MFIX
%   run (variable <run>) to create movie frames of the evolution of gas
%   temperature during an eruption, which are passed as output for a file
%   gTemp_<run>.avi.
%   Taryn Black, last edit 24 April 2015 

    cd(sprintf('%d',run))

    mGT = moviein(timesteps);
        hold on
        axis equal
        set(gca,'XTick',distance(2:end)*1000,'XTickLabel',distance(2:end),'FontSize',12)
            xlabel('\bf Distance (km)','FontSize',12)
        set(gca,'YTick',altitude(2:end)*1000,'YTickLabel',altitude(2:end),'FontSize',12)
            ylabel('\bf Altitude (km)','FontSize',12)
        hc = colorbar;
            ylabel(hc,'\bf Gas temperature (K)','FontSize',12)
            set(hc,'FontSize',12)
            caxis(T_range)

    for t = 1:timesteps
        
        TG = reshape(T_G((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);
        TG = TG(2:end-1,2:end-1);
        contourf(X,Y,TG)
        axis tight

        if strcmp(PULSE,'T') == 1
            str = 'Unsteady flow';
            tL = ({sprintf('ID#%d: %s, %.1f Hz',run,str,FREQ);sprintf('Temperature, t=%d s',time(t))});
        elseif strcmp(PULSE,'F') == 1
            str = 'Steady flow';
            tL = ({sprintf('ID#%d: %s',run,str);sprintf('Temperature, t=%d s',time(t))});
        end
        title(tL,'FontSize',12,'FontWeight','bold');
                    
        mGT(:,t) = getframe(gcf);
        
        if t == timesteps;
            saveas(gcf,sprintf('gTempEND_%d.jpg',run))
        end
    end
    
% cd ..

end


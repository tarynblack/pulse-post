function [mBD,ROB] = bulkdens2D(EP_G,P_G,T_G,RO_S1,RO_S2,timesteps,IMAX,JMAX,X,Y,distance,altitude,run,BD_range,PULSE,FREQ,time)
% bulkdens makes frames of the time evolution of bulk density during
% a simulated volcanic eruption.
%   This script processes volume fraction, pressure, temperature, and
%   density data from a specified MFIX run (variable <run>) to create movie
%   frames of the evolution of bulk density during an eruption, which are
%   passed as output for a file bulkDens_<run>.avi. 
%   Taryn Black, last edit 24 April 2015

cd(sprintf('%d',run))

    mBD = moviein(timesteps);
        hold on
        axis equal
        set(gca,'XTick',distance(2:end)*1000,'XTickLabel',distance(2:end),'FontSize',12)
            xlabel('\bf Distance (km)','FontSize',12)
        set(gca,'YTick',altitude(2:end)*1000,'YTickLabel',altitude(2:end),'FontSize',12)
            ylabel('\bf Altitude (km)','FontSize',12)
        hc = colorbar;
            ylabel(hc,'\bf Bulk density (kg/m^3)','FontSize',12)
            set(hc,'FontSize',12)
            caxis(BD_range)

    for t = 1:timesteps

        EPG = reshape(EP_G((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);
        PG  = reshape(P_G((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);
        TG  = reshape(T_G((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);
        
        % calculate bulk density and delete edges (they blow up contour)
        ROB = (EPG.*(PG./(461.5*TG))) + ((1 - EPG).*((RO_S1/2) + (RO_S2/2)));
        ROB = ROB(2:end-1,2:end-1);
        contourf(X,Y,ROB)
        axis tight
        
        if strcmp(PULSE,'T') == 1
            str = 'Unsteady flow';
            tL = ({sprintf('ID#%d: %s, %.1f Hz',run,str,FREQ);sprintf('Bulk density, t=%d s',time(t))});
        elseif strcmp(PULSE,'F') == 1
            str = 'Steady flow';
            tL = ({sprintf('ID#%d: %s',run,str);sprintf('Bulk density, t=%d s',time(t))});
        end
        title(tL,'FontSize',12,'FontWeight','bold');
        
        mBD(:,t)=getframe(gcf);
        
        if t == timesteps;
            saveas(gcf,sprintf('bulkDensEND_%d.jpg',run))
        end   
    end

% cd ..
    
end

function [mGVel] = velocity2D(V_G,U_G,timesteps,IMAX,JMAX,X,Y,distance,altitude,run,V_range,PULSE,FREQ,time)
% velocity2D makes frames of the time evolution of gas velocity during
% a simulated volcanic eruption.
%   This script processes gas velocity (V_G) data from a specified MFIX
%   run (variable <run>) to create movie frames of the evolution of gas
%   velocity during an eruption, which are passed as output for a file
%   gPres_<run>.avi.
%   Taryn Black, last edit 24 April 2015

%     cd(sprintf('%d',run))
    
    XX = repmat(X,length(X),1);
    YY = repmat(Y',1,length(Y));
    
    velmag = zeros(JMAX,IMAX);

    mGVel = moviein(timesteps);
        hold on
        axis equal
        hvel = quiver(0,0,'b');
%         hmag = contour([0 0],[0 0],[0 0;0 0]);
        set(gca,'XTick',distance(2:end)*1000,'XTickLabel',distance(2:end),'FontSize',12)
            xlabel('\bf Distance (km)','FontSize',12)
        set(gca,'YTick',altitude(2:end)*1000,'YTickLabel',altitude(2:end),'FontSize',12)
            ylabel('\bf Altitude (km)','FontSize',12)
        hc = colorbar;
            ylabel(hc,'\bf Velocity (m/s)','FontSize',12)
            set(hc,'FontSize',12)
            caxis(V_range)

    for t = 1:timesteps
        
        delete(hvel)
        
        VG = reshape(V_G((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);
        UG = reshape(U_G((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);
        
%         for i = 1:JMAX
%             for j = 1:IMAX
%                 velmag(i,j) = norm([VG(i,j) UG(i,j)]);
%             end
%         end
        
        hvel = quiver(XX,YY,UG,VG,2,'b');
%         hmag = contour(X,Y,velmag);
        axis tight

        if strcmp(PULSE,'T') == 1
            str = 'Unsteady flow';
            tL = ({sprintf('ID#%d: %s, %.1f Hz',run,str,FREQ);sprintf('Gas velocity, t=%d s',time(t))});
        elseif strcmp(PULSE,'F') == 1
            str = 'Steady flow';
            tL = ({sprintf('ID#%d: %s',run,str);sprintf('Gas velocity, t=%d s',time(t))});
        end
        title(tL,'FontSize',12,'FontWeight','bold');
        
        drawnow
        
        mGVel(:,t)=getframe(gcf);
        
        if time(t) == 300;
            saveas(gcf,sprintf('gVelo300s_%d.jpg',run))
        end
        if time(t) == 450;
            saveas(gcf,sprintf('gVelo450s_%d.jpg',run))
        end
        if t == timesteps;
            saveas(gcf,sprintf('gVeloEND_%d.jpg',run))
        end
        
    end

% cd ..
    
end


function [mEntr,hFig2,hFig3,avg_entrain,plumevolume,vol_inlet] = entrainment2D(EP_G,V_G,U_G,run,timesteps,distance,altitude,XRES,YRES,IMAX,JMAX,PULSE,FREQ,time)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    cd(sprintf('%d',run))

    EPG_inlet   = zeros(1,length(timesteps));
    avg_entrain = zeros(1,length(timesteps));
    plumevolume = zeros(1,length(timesteps));
    
    if strcmp(PULSE,'T') == 1
        str = sprintf('ID#%d: Unsteady flow %.1f Hz',run,FREQ);
    elseif strcmp(PULSE,'F') == 1
        str = sprintf('ID#%d: Steady flow',run);
    end

    edge = 0.99999; % gas volume fraction at edge of plume (## ARBITRARY ##)
    cellvol = XRES*YRES*1;  % volume of a grid cell [m3]

    %%% Figure initialization: set constant properties (labels, colorbar,
    %%% legend, etc.) before loop to reduce runtime
    figEntr = figure('Name','Entrainment animation','units','inches','outerposition',[1 1 9 9]);
        mEntr = moviein(timesteps,figEntr);    
        hold on
        axis equal
        xlim([0 distance(end)*1000])
        ylim([0 altitude(end)*1000])
        set(gca,'XTick',distance(2:end)*1000,'XTickLabel',distance(2:end),'FontSize',12)
            xlabel('\bf Distance (km)','FontSize',12)
        set(gca,'YTick',altitude(2:end)*1000,'YTickLabel',altitude(2:end),'FontSize',12)
            ylabel('\bf Altitude (km)','FontSize',12)
        set(gca,'FontSize',12)
%         grid on
%         grid(gca,'minor')
%         hedge = plot(0,0,'k');
%         hnorm = quiver(0,0,'b');
%         hvelo = quiver(0,0,'r');
%         hentr = surface([0 0;0 0],[0 0;0 0],[0 0;0 0],[0 0;0 0]);
        hentr = scatter(0,0,[],0);
        hc = colorbar;
            set(hc,'FontSize',12)
            ylabel(hc,'\bf Entrainment coefficient','FontSize',12)
            caxis([-0.1 1.1])    %%% ##### FIXME ##### %%%
            
%         legend([hedge hnorm hvelo],{sprintf('Gas volume fraction = %f',edge),...
%             'Surface unit normals','Plume-normal velocities'},'Location',...
%             'NorthEast','FontSize',12)
%         legend([hedge hvelo],{sprintf('Gas volume fraction = %f',edge),...
%             'Plume-normal velocities'},'Box','off','Location','North','Orientation','Horizontal','FontSize',12)
%             set(gca,'LegendColorbarListeners',[]);
%             setappdata(gca,'LegendColorbarManualSpace',1);
%             setappdata(gca,'LegendColorbarReclaimSpace',1);

    for t = 1:timesteps

        %%% Find gas volume fraction field at current timestep and filter out
        %%% values along domain edges (nonphysical)
            EPG = reshape(EP_G((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);
            EPG_inlet(t) = EPG(1,IMAX/2);
            EPG = EPG(2:end-1,2:end-1);
            [rowsEPG colsEPG] = size(EPG);

        %%% Extract the coordinates of the plume edge isoline. contourc returns
        %%% a 2-row matrix of coordinates, row1=xcoord (column!), row2=ycoord
        %%% (row!). Coordinates are preceded by a column containing the contour
        %%% value and number of points within that value - this column is
        %%% removed in the loop below, to isolate contour coordinates.
            iso = contourc(EPG,[edge edge]);
            coords = zeros(2,length(iso));
            for i = 1:length(iso)
                if iso(1,i) == edge
                    coords(:,i) = 0;
                elseif iso(1,i) ~= edge
                    coords(:,i) = iso(:,i);
                end
            end
            coords(:,all(coords==0,1))=[];
            x = coords(1,:);
            y = coords(2,:);

        %%% Calculate normals to plume surface. Each column in 'normal' is a
        %%% single normal vector.
            dx = gradient(x);
            dy = gradient(y);
            normal = [-dy; dx];
    %         figure(1)
    %             plot(x,y,'k')
    %             hold on
    %             quiver(x,y,-dy,dx,'b')

        %%% Convert plot (decimal) coordinates to nearest (integer) grid coordinates
            cells = round(coords);
            row = cells(2,:);
            col = cells(1,:);

        %%% Calculate vertical and horizontal velocity components at isosurface
            VG = reshape(V_G((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);
            UG = reshape(U_G((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);

            xvel = zeros(1,length(cells));          % x component of plume velocity
            yvel = zeros(1,length(cells));          % y component of plume velocity
            xvelpos = zeros(1,length(cells));       % x position of velocity
            yvelpos = zeros(1,length(cells));       % y position of velocity

            % This loop determines whether the isosurface grid coordinates
            % (which have been rounded from the exact isosurface coordinates)
            % are inside or outside the actual plume edge. If the grid
            % coordinates are outside the plume, the loop checks the
            % coordinates the the left and right {NOTE: only checking along
            % x-axis!}. At the first gridpoint *actually* inside the plume, the
            % x/y velocities and associated gridpoints are recorded.
            for j = 1:length(row)
                if EPG(row(j),col(j)) <= edge
                    yvel(j) = VG(row(j),col(j));
                    xvel(j) = UG(row(j),col(j));
                    yvelpos(j) = row(j);
                    xvelpos(j) = col(j);
                elseif col(j)+1<=colsEPG && EPG(row(j),col(j)+1) <= edge
                    yvel(j) = VG(row(j),col(j)+1);
                    xvel(j) = UG(row(j),col(j)+1);
                    yvelpos(j) = row(j);
                    xvelpos(j) = col(j)+1;
                elseif col(j)-1>=1 && EPG(row(j),col(j)-1) <= edge
                    yvel(j) = VG(row(j),col(j)-1);
                    xvel(j) = UG(row(j),col(j)-1);
                    yvelpos(j) = row(j);
                    xvelpos(j) = col(j)-1;
                end
            end
            
        %%% Convert grid coordinates back to physical coordinates for
        %%% plotting
            xvelpos = xvelpos*XRES;
            yvelpos = yvelpos*YRES;

        %%% Convert normals to unit normals and calculate the velocity
        %%% component in the surface-normal direction
            unitnorm = zeros(2,length(normal));
            for k = 1:length(normal)
                unitnorm(:,k) = normal(:,k)./norm(normal(:,k));
            end
            normvel_mag = dot([xvel;yvel],unitnorm);    % magnitude of velocity normal to plume isosurface (scalar)
            normvel_x = normvel_mag.*unitnorm(1,:);     % x-component of plume-normal velocity
            normvel_y = normvel_mag.*unitnorm(2,:);     % y-component of plume-normal velocity

        %%% Calculate entrainment coefficient (Morton linear assumption:
        %%% u_norm=k*u_vert) where k is calculated in each cell and averaged
        %%% for the entire plume. 
        %%% ##### TODO ##### entrain_coeff has Inf, NaN, numbers outside of 0-1 range, problem???
            entrain_coeff = normvel_mag./yvel;
            entrain_coeff_cut = entrain_coeff(~isinf(entrain_coeff));   %%% ##### FIXME #####
            avg_entrain(t) = nanmean(entrain_coeff_cut);                %%% ##### FIXME #####

        %%% Plot results, updating for each timestep: plume edge, surface unit
        %%% normals, plume-normal velocities, entrainment coefficient along
        %%% (grid) plume edge
            color = entrain_coeff;  % base colorbar on entrainment coefficient

%             delete(hedge,hnorm,hvelo,hentr)
%             delete(hedge,hvelo,hentr)
            delete(hentr)

%             hedge = plot(x,y,'k');
%             hnorm = quiver(x,y,unitnorm(1,:),unitnorm(2,:),'b');
%             hvelo = quiver(xvelpos,yvelpos,normvel_x,normvel_y,3,'r');
            hentr = scatter(xvelpos,yvelpos,20,entrain_coeff,'filled');
%             hentr = surface([xvelpos;xvelpos],[yvelpos;yvelpos],[entrain_coeff;entrain_coeff],[color;color]);
%                     set(hentr,'FaceColor','none','EdgeColor','interp','Marker','.','MarkerSize',10)

            tL = ({sprintf('%s',str);sprintf('Plume entrainment, t=%d s',time(t))});
            title(tL,'FontSize',12,'FontWeight','bold');
                    
            drawnow

            mEntr(:,t) = getframe(gcf);

        %%% Calculate total plume volume at each timestep   
            numcells = sum(sum(EPG<=edge));
            plumevolume(t) = numcells*cellvol;

    end
    
    hold off
    
    %%% Plot plume-averaged entrainment coefficient over time
        hFig2 = figure('Name','Entrainment variation','units','normalized','outerposition',[0 0 1 0.5],'visible','on');
            plot(time,avg_entrain);
            title(sprintf('%s: Variations in entrainment coefficient',str),'FontWeight','bold','FontSize',12)
            xlabel('Timestep','FontWeight','bold','FontSize',12)
            ylabel('Plume-averaged entrainment coefficient','FontWeight','bold','FontSize',12)
            ylim([-10 10])
                set(gca,'FontSize',12)

    %%% Calculate inlet volume flow rate
        vol_inlet = EPG_inlet*cellvol;
        sum_vol_inlet = cumsum(vol_inlet);
        
    %%% Calculate expected plume volume based on literature entrainment
        litk = 0.13;     % literature typical entrainment coeff (check)
        lit_volume = zeros(1,length(vol_inlet));
        lit_volume(1) = vol_inlet(1);
        for k = 2:length(lit_volume)
            lit_volume(k) = (1+litk)*(lit_volume(k-1) + vol_inlet(k));
        end
            
    %%% Plot total plume volumes (input, lit expected, modeled) and
    %%% entrained volumes (lit expected, modeled)
        hFig3 = figure('Name','Total and entrained volumes','units','normalized','outerposition',[0 0 1 1],'visible','on');
            h3vol = subplot(2,1,1);
                semilogy(time,plumevolume,'b',time,lit_volume,'r',time,sum_vol_inlet,'-k')
                title(h3vol,{sprintf('%s: Total plume volumes',str)},'FontWeight','bold','FontSize',12)
                xlabel(h3vol,{'Time (s)'},'FontWeight','bold','FontSize',12)
                ylabel(h3vol,{'Volume (m^3)'},'FontWeight','bold','FontSize',12)
                ylim([0 prod(size(EPG))*cellvol])
                    set(gca,'FontSize',12)
                legend(h3vol,{'Modeled growth','Linear growth','Cumulative input'},'Box','on','Location','NorthWest','FontWeight','bold','FontSize',12)
            h3diff = subplot(2,1,2);
                semilogy(time,plumevolume-sum_vol_inlet,'b',time,lit_volume-sum_vol_inlet,'r')
                title(h3diff,{sprintf('%s: Entrained volumes',str)},'FontWeight','bold','FontSize',12)
                xlabel(h3diff,{'Time (s)'},'FontWeight','bold','FontSize',12)
                ylabel(h3diff,{'Volume (m^3)'},'FontWeight','bold','FontSize',12)
                ylim([0 prod(size(EPG))*cellvol])
                    set(gca,'FontSize',12)
                legend(h3diff,{'Modeled entrainment','Linear entrainment'},'Box','on','Location','NorthWest','FontWeight','bold','FontSize',12)

% figure(4)
% plot(time,lit_volume,time,sum_vol_inlet,time,vol_inlet,time,lit_volume-sum_vol_inlet)
% legend('lit volume','sum inlet','inlet','entrainment')

end


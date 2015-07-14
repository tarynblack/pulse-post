function [h1,h2,h3] = momentum2D(h,EP_G,RO_G,ROP_S1,ROP_S2,V_G,V_S1,V_S2,IMAX,JMAX,LENGTH,XRES,YRES,ventwidth,timesteps,run)
% momentum plots the time evolution of the momentum density of gas and
% solid phases during a simulated volcanic eruption.
%   This script calculates and plots the momentum density from a specified
%   MFIX run (var <run>) at selected heights above the volcanic vent (var
%   <h>), both directly over the vent summed over the entire width of the
%   domain. 
%   Taryn Black, last edit 10 November 2014

cd(sprintf('%d',run))

    time = 1:timesteps;   

    %%% Set vent boundaries in domain units
    ventcells = ventwidth/XRES;
    L_vent = 0.5*((LENGTH/XRES)-ventcells); % LHS cell boundary of vent
    R_vent = 0.5*((LENGTH/XRES)+ventcells); % RHS cell boundary of vent

    %%% Preallocate height-data and momentum flux vectors
    epg_h = zeros(length(h),IMAX,timesteps);
    rog_h = zeros(length(h),IMAX,timesteps);
    rops1_h = zeros(length(h),IMAX,timesteps);
    rops2_h = zeros(length(h),IMAX,timesteps);
    vg_h = zeros(length(h),IMAX,timesteps);
    vs1_h = zeros(length(h),IMAX,timesteps);
    vs2_h = zeros(length(h),IMAX,timesteps);

    J_G = zeros(length(h),timesteps);
    J_Gvent = zeros(length(h),timesteps);
    J_S1 = zeros(length(h),timesteps);
    J_S2 = zeros(length(h),timesteps);

    %%% Reshape data into domain matrices over time and calculate fluxes
    for t = 1:timesteps
        EPG = reshape(EP_G((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);
        ROG = reshape(RO_G((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);
        ROPS1 = reshape(ROP_S1((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);
        ROPS2 = reshape(ROP_S2((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);
        VG = reshape(V_G((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);
        VS1 = reshape(V_S1((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);
        VS2 = reshape(V_S2((t-1)*IMAX*JMAX+1:t*IMAX*JMAX),[JMAX IMAX]);

        % Select data at a given height and calculate fluxes
        for i = 1:length(h)
            epg_h(i,:,t) = EPG(1+(h(i)/YRES),:);
            rog_h(i,:,t) = ROG(1+(h(i)/YRES),:);
            rops1_h(i,:,t) = ROPS1(1+(h(i)/YRES),:);
            rops2_h(i,:,t) = ROPS2(1+(h(i)/YRES),:);
            vg_h(i,:,t) = VG(1+(h(i)/YRES),:);
            vs1_h(i,:,t) = VS1(1+(h(i)/YRES),:);
            vs2_h(i,:,t) = VS2(1+(h(i)/YRES),:);
            
            J_G(i,t) = sum(epg_h(i,:,t).*rog_h(i,:,t).*vg_h(i,:,t));
            J_Gvent(i,t) = sum(epg_h(i,L_vent:R_vent,t).*rog_h(i,L_vent:R_vent,t).*vg_h(i,L_vent:R_vent,t));
            J_S1(i,t) = sum(rops1_h(i,:,t).*vs1_h(i,:,t));
            J_S2(i,t) = sum(rops2_h(i,:,t).*vs2_h(i,:,t));
        end
    end

    %%% Plot gas and solid momentum densities
    h1 = figure(1);
        subplot(2,2,[1 2])
            plot(time,J_G(1,:),'r',time,J_G(2,:),'b',time,J_G(3,:),'g')
            xlim([0 timesteps])
            title(sprintf('ID:%d: \\SigmaGas momentum density over domain width',run))
            xlabel('timestep')
            ylabel('momentum density, N/m^3')
            hleg1a = legend(sprintf('h=%d m',h(1)),sprintf('h=%d m',h(2)),sprintf('h=%d m',h(3)),'Location','EastOutside');
            set(get(hleg1a,'Title'),'String','Height above vent (m)')
        subplot(2,2,[3 4])
            plot(time,J_Gvent(1,:),'r',time,J_Gvent(2,:),'b',time,J_Gvent(3,:),'g')
            xlim([0 timesteps])
            title(sprintf('ID:%d: \\SigmaGas momentum density across vent',run))
            xlabel('timestep')
            ylabel('momentum density, N/m^3')
            hleg1b = legend(sprintf('h=%d m',h(1)),sprintf('h=%d m',h(2)),sprintf('h=%d m',h(3)),'Location','EastOutside');
            set(get(hleg1b,'Title'),'String','Height above vent (m)')
    h2 = figure(2);
        subplot(2,2,[1 2])
            plot(time,J_S1(1,:),'r',time,J_S1(2,:),'b',time,J_S1(3,:),'g')
            xlim([0 timesteps])
            title(sprintf('ID:%d: \\SigmaSolid1 momentum density over domain width',run))
            xlabel('timestep')
            ylabel('momentum density, N/m^3')
            hleg2a = legend(sprintf('h=%d m',h(1)),sprintf('h=%d m',h(2)),sprintf('h=%d m',h(3)),'Location','EastOutside');
            set(get(hleg2a,'Title'),'String','Height above vent (m)')
        subplot(2,2,[3 4])
            plot(time,J_S2(1,:),'r',time,J_S2(2,:),'b',time,J_S2(3,:),'g')
            xlim([0 timesteps])
            title(sprintf('ID:%d: \\SigmaSolid2 momentum density over domain width',run))
            xlabel('timestep')
            ylabel('momentum density, N/m^3')
            hleg2b = legend(sprintf('h=%d m',h(1)),sprintf('h=%d m',h(2)),sprintf('h=%d m',h(3)),'Location','EastOutside');
            set(get(hleg2b,'Title'),'String','Height above vent (m)')
    h3 = figure(3);
            plot(time,J_G(1,:),'r',time,J_G(2,:),'b',time,J_G(3,:),'g',time,J_Gvent(1,:),'r--',time,J_Gvent(2,:),'b--',time,J_Gvent(3,:),'g--')
            xlim([0 timesteps])
            title(sprintf('ID:%d: \\SigmaGas momentum density across domain (solid), vent (dashed)',run))
            xlabel('timestep')
            ylabel('momentum density, N/m^3')
            hleg3 = legend(sprintf('h=%d m',h(1)),sprintf('h=%d m',h(2)),sprintf('h=%d m',h(3)),sprintf('h=%d m',h(1)),sprintf('h=%d m',h(2)),sprintf('h=%d m',h(3)),'Location','EastOutside');
            set(get(hleg3,'Title'),'String','Height above vent (m)')

% cd ..

end


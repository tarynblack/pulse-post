function [ EPG_vent,EPS_vent ] = evolveVolumeFraction( MING,MAXG,FREQ,...
    PULSE,DT,TSTOP,vis,titlerun,imtype,time,savepath,postpath,run )
%evolveVolumeFraction Summary of this function goes here
%   Detailed explanation goes here

    cd(savepath)
    delete('GasSolFrac*')

    % Create time vector
      densetime = 0:(FREQ/20):TSTOP;    

    % Initialize figure
      varVFE = 'Volume fraction of gas and solids';
      figVolFrac = figure('Name','Volume fraction evolution','units',...
          'centimeters','outerposition',[0 0 33.33 18.75],'visible',vis,...
          'PaperPositionMode','auto','color','w');
      axVolFrac = axes('Parent',figVolFrac,'box','on','TickDir','in','FontSize',12);
      hold on
      grid(axVolFrac,'on');
      [axVolFrac,gasFrac,solFrac] = plotyy(0,0,0,0);
      hold on
      set(axVolFrac(1),'Parent',figVolFrac,'box','on','TickDir','in',...
          'FontSize',12,'YLim',[0.995 1],'YTick',0.995:0.001:1,...
          'XLim',[0 TSTOP]);
      set(axVolFrac(2),'Parent',figVolFrac,'TickDir','in',...
          'FontSize',12,'YLim',[0 0.005],'YTick',0:0.001:0.005,...
          'YTickLabel',{'0';'0.001';'0.002';'0.003';'0.004';'0.005'},...
          'XLim',[0 TSTOP]);
      xlabel(axVolFrac(1),'\bfTime (s)');
      ylabel(axVolFrac(1),'\bfGas volume fraction');
      ylabel(axVolFrac(2),'\bfSolid volume fraction');
      
    % Calculate solid and gas volume fractions (clipped) with time
      EPG_vent = MING*(1-MING)*abs(sin(pi*FREQ*densetime))+MING;
      EPG_vent(EPG_vent>MAXG) = 1;
      EPS_vent = 1 - EPG_vent;
      
    % Step through time to create evolving figure  
      for t = 1:length(time)
          
          % Create figure for current timestep
            figure(figVolFrac)
            set(gasFrac,'XData',densetime(1:(t-1)*DT*20/FREQ),...
                'YData',EPG_vent(1:(t-1)*DT*20/FREQ),'LineWidth',3);
            set(solFrac,'XData',densetime(1:(t-1)*DT*20/FREQ),...
                'YData',EPS_vent(1:(t-1)*DT*20/FREQ),'LineWidth',3);
            cd(postpath)
            tVFE = pulsetitle(varVFE,PULSE,time,t,titlerun,FREQ);
            title(axVolFrac(1),tVFE,'FontWeight','bold');
            set(figVolFrac,'Visible',vis);
           
          % Save current timestep
            cd(savepath)
            figVF = 'GasSolFracCurrent.jpg';
            saveas(figVolFrac,fullfile(savepath,figVF));
            imgVF = imread(figVF);
            if strcmp(imtype,'tif') == 1 || strcmp(imtype,'tiff') == 1
                imwrite(imgVF,fullfile(savepath,sprintf('GasSolFrac_tsteps_%s.tif',...
                    run)),'tif','WriteMode','append');
            else
                saveas(figVF,fullfile(savepath,sprintf('GasSolFrac_t%03d_%s.%s',...
                    densetime(t),run,imtype)));
            end
          
      end
      
      cd(postpath)
      
end


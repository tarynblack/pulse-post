function [ STavg_entr,STavg_jentr,massratio] = emergency( run,time,...
    PULSE,FREQ,titlerun,vis,DT,postpath,savepath,massflux_crange,...
    massflux_legend,massflux_alts,JMAX,ghostcells,ticky,YRES,labely,...
    labelYunit,jetheight )

cd(savepath)

plumevolume = load(sprintf('plumevolume_%s.txt',run));
ecoeff = load(sprintf('coeff_avg-std_%s.txt',run));
avg_coeff = ecoeff(:,1);
std_coeff = ecoeff(:,2);
expn = load(sprintf('expn_avg-std_%s.txt',run));
avg_expn = expn(:,1);
std_expn = expn(:,2);
entr = load(sprintf('entr_avg-std_%s.txt',run));
avg_entr = entr(:,1);
std_entr = entr(:,2);
jcoeff = load(sprintf('jcoeff_avg-std_%s.txt',run));
avg_jcoeff = jcoeff(:,1);
std_jcoeff = jcoeff(:,2);
jexpn = load(sprintf('jexpn_avg-std_%s.txt',run));
avg_jexpn = jexpn(:,1);
std_jexpn = jexpn(:,2);
jentr = load(sprintf('jentr_avg-std_%s.txt',run));
avg_jentr = jentr(:,1);
std_jentr = jentr(:,2);
e_Mconic = load(sprintf('e_MortonConic_%s.txt',run));
STavg_entr  = mean(avg_entr(2:end));
STavg_jentr = mean(avg_jentr(2:end));
netmassflux = load(sprintf('netmassflux_%s.txt',run));
netMF_alts = netmassflux(:,massflux_alts);
avgNMF = mean(netmassflux(:,2:end),1);
massratio = load(sprintf('finalMassRatio_%s.txt',run));
massratio = massratio(3);


% ------------------- ENTRAINMENT TIME SERIES PLOTS ------------------- %
    if strcmp(PULSE,'T') == 1
      str = sprintf('%s: Unsteady flow %g Hz',titlerun,FREQ);
    elseif strcmp(PULSE,'F') == 1
      str = sprintf('%s: Steady flow',titlerun);
    end
    
  % Total plume volume and change in plume volume   
    figPlumeVol = figure('Name','Plume Volume','visible',vis,'units',...
        'centimeters','outerposition',[0 0 33.33 18.75],'PaperPositionMode',...
        'auto','color','w');
    axVol1 = subplot(2,1,1);
      plot(time(2:end),plumevolume)
      title(axVol1,{sprintf('Total plume volume\n%s',str)},...
          'FontWeight','bold')
      xlabel(axVol1,{'Time (s)'},'FontWeight','bold')
      ylabel(axVol1,{'Volume (m^3)'},'FontWeight','bold')
    axVol2 = subplot(2,1,2);
      plot(time(2:end),[0 diff(plumevolume)'])
      title(axVol2,{sprintf('Change in plume volume\n%s',str)},...
          'FontWeight','bold')
      xlabel(axVol2,{'Time (s)'},'FontWeight','bold')
      ylabel(axVol2,{'\DeltaVolume (m^3)'},'FontWeight','bold')
    set([axVol1 axVol2],'box','on','FontSize',12)
    grid(axVol1,'on'); grid(axVol2,'on')
    saveas(figPlumeVol,fullfile(savepath,sprintf('PlumeVolume_%s.jpg',run)));
    
  % Plume-averaged entrainment/expansion
    figCoeff = figure('Name','Entrainment Coefficients','visible',vis,...
        'units','centimeters','outerposition',[0 0 33.33 18.75],...
        'PaperPositionMode','auto','color','w');
    axCoeff = axes('Parent',figCoeff,'box','on','FontSize',12);
    grid(axCoeff,'on');
    axCoeff.YMinorGrid = 'on';
    hold on
    for t = 2:length(time)
        coeff_all = load(sprintf('ecoeff_all_t%03d.txt',time(t)));
        hs = scatter(time(t)*ones(1,length(coeff_all)),coeff_all,...
            'MarkerEdgeColor',[0.65 0.65 0.65],'SizeData',16);
        hold on
    end
    hs.DisplayName = 'Full plume';
    he1 = errorbar(time(2:end),avg_expn,std_expn,'r',...
        'LineWidth',4,'LineStyle','none','Marker','+','MarkerSize',10,...
        'DisplayName','Expansion');
    he2 = errorbar(time(2:end),avg_entr,std_entr,'c',...
        'LineWidth',4,'LineStyle','none','Marker','+','MarkerSize',10,...
        'DisplayName','Entrainment');
    he3 = errorbar(time(2:end),avg_coeff,std_coeff,'k',...
        'LineWidth',1.5,'Marker','+','MarkerSize',10,...
        'DisplayName','Combined');
    tlCoeff1 = 'EEE averaged over full plume surface';
    tlCoeff2 = sprintf('%s',str); 
    tlCoeff3 = sprintf('Spatiotemporally averaged entrainment = %.4f',STavg_entr);
    title(axCoeff,{tlCoeff1;tlCoeff2;tlCoeff3},'FontWeight','bold')
    xlabel(axCoeff,'\bfTime (s)')
    ylabel(axCoeff,'\bfEEE')
    xlim(axCoeff,[0 time(end)+DT])
    ylim(axCoeff,[-1 1])
    line(xlim,[0 0],'color',[0.1 0.1 0.1],'LineWidth',0.5);
    hl = legend([hs he1 he2 he3],'Location','EastOutside');
    hl.Box = 'on';
    hl.FontWeight = 'bold';
    hl.FontSize = 10;
    saveas(figCoeff,fullfile(savepath,sprintf('Coefficients_%s.jpg',run)));
    
  % Plume-averaged entrainment/expansion, zoomed to show detail
    ylim(axCoeff,[-0.2 0.2])
    saveas(figCoeff,fullfile(savepath,sprintf('CoefficientsDetail_%s.jpg',run)));
    
  % Jet region entrainment/expansion
    figJCoeff = figure('Name','Jet Region Entrainment Coefficients',...
        'visible',vis,'units','centimeters','outerposition',...
        [0 0 33.33 18.75],'PaperPositionMode','auto','color','w');
    axJCoeff = axes('Parent',figJCoeff,'box','on','FontSize',12);
    grid(axJCoeff,'on');
    axJCoeff.YMinorGrid = 'on';
    hold on
    for t = 2:length(time)
        jcoeff_all = load(sprintf('jcoeff_all_t%03d.txt',time(t)));
        hjs = scatter(time(t)*ones(1,length(jcoeff_all)),jcoeff_all,...
            'MarkerEdgeColor',[0.65 0.65 0.65],'SizeData',16);
        hold on
    end
    hjs.DisplayName = 'Jet region';
    hje1 = errorbar(time(2:end),avg_jexpn,std_jexpn,'r',...
        'LineWidth',4,'LineStyle','none','Marker','+','MarkerSize',10,...
        'DisplayName','Expansion');
    hje2 = errorbar(time(2:end),avg_jentr,std_jentr,'c',...
        'LineWidth',4,'LineStyle','none','Marker','+','MarkerSize',10,...
        'DisplayName','Entrainment');
    hje3 = errorbar(time(2:end),avg_jcoeff,std_jcoeff,'k',...
        'LineWidth',1.5,'Marker','+','MarkerSize',10,...
        'DisplayName','Combined');
    tlJCoeff1 = 'EEE averaged over jet region';
    tlJCoeff2 = sprintf('%s',str); 
    tlJCoeff3 = sprintf('Spatiotemporally averaged entrainment = %.4f',STavg_jentr);
    title(axJCoeff,{tlJCoeff1;tlJCoeff2;tlJCoeff3},'FontWeight','bold')
    xlabel(axJCoeff,'\bfTime (s)')
    ylabel(axJCoeff,'\bfEEE')
    xlim(axJCoeff,[0 time(end)+DT])
    ylim(axJCoeff,[-1 1])
    line(xlim,[0 0],'color',[0.1 0.1 0.1],'LineWidth',0.5);
    hjl = legend([hjs hje1 hje2 hje3],'Location','EastOutside');
    hjl.Box = 'on';
    hjl.FontWeight = 'bold';
    hjl.FontSize = 10;
    saveas(figJCoeff,fullfile(savepath,sprintf('JetCoefficients_%s.jpg',run)));
    
  % Jet region entrainment/expansion, zoomed to show detail
    ylim(axJCoeff,[-0.2 0.2])
    saveas(figJCoeff,fullfile(savepath,sprintf('JetCoefficientsDetail_%s.jpg',run)));
    
  % Comparison conic coefficient
    figMorton = figure('Name','Conic entrainment coefficient','units',...
        'centimeters','outerposition',[0 6 33.33 11.25],'visible',vis,...
        'PaperPositionMode','auto','color','w');
    axMor = axes('Parent',figMorton,'box','on','FontSize',12); 
    plot(axMor,time(2:end),e_Mconic)
    title(axMor,sprintf('Morton entrainment coefficient\n%s',str),...
        'FontWeight','bold')
    xlabel(axMor,'\bfTime (s)')
    ylabel(axMor,'\bfCoefficient')
    ylim([0 1])
    grid(axMor,'on');
    saveas(figMorton,fullfile(savepath,sprintf('MortonConic_%s.jpg',run)));
  % =================================================================== %

  % ------------------ NET MASS FLUX TIME SERIES PLOTS ------------------ %
  figNetMF = figure('Name','Net Mass Flux','visible',vis,'units',...
      'centimeters','outerposition',[0 0 33.33 18.75],'PaperPositionMode',...
      'auto','color','w');
  axNetMF = axes('Parent',figNetMF,'box','on','FontSize',12);
  grid(axNetMF,'on');
  axis(axNetMF,[0,time(end),2*massflux_crange(1),2*massflux_crange(2)]);
  hold on
  plot(time(2:end),netMF_alts);
  legend(axNetMF,massflux_legend)
  title(axNetMF,sprintf('Net solid mass flux\n%s',str))
  xlabel(axNetMF,'\bfTime (s)')
  ylabel(axNetMF,'\bfNet mass flux (kg/m^2s)')
  saveas(figNetMF,fullfile(savepath,sprintf('NetMassFlux_tseries_%s.jpg',run)))
  
  % ---------------- TIME-AVERAGED MASS FLUX WITH HEIGHT ---------------- %
  % Average mass flux time series: figure and axes properties
    figAvgMFZ = figure('Name','Spatially averaged mass flux with altitude',...
        'units','centimeters','outerposition',[0 0 33.33 18.75],'visible',...
        vis,'PaperPositionMode','auto','color','w');
    axAvgMFZ = axes('Parent',figAvgMFZ,'box','on','TickDir','in','FontSize',12);
    hold on
    grid(axAvgMFZ,'on');
    axis(axAvgMFZ,[2*massflux_crange(1),2*massflux_crange(2),0,JMAX-ghostcells]);
    set(axAvgMFZ,'YTick',ticky(2:end)/YRES,'YTickLabel',labely);
    xlabel(axAvgMFZ,'\bfNet mass flux (kg/m^2s)')
    ylabel(axAvgMFZ,sprintf('\\bfAltitude (%s)',labelYunit))
    hMFZ = plot(0,0,'DisplayName','Previous profiles');
    hBlk = plot(0,0,'k','LineWidth',2,'DisplayName',...
        'Current profile');
    hJet = plot(2*massflux_crange(1):1E3:2*massflux_crange(2),...
        jetheight*ones(1,length(2*massflux_crange(1):1E3:2*massflux_crange(2))),...
        '--','Color',[0.2 0.5 0.2],'LineWidth',1.5,'DisplayName',...
        sprintf('Jet height (%.3f km)',jetheight*YRES/1000));
    hMFZleg = legend(axAvgMFZ,[hBlk hMFZ hJet]);
    set(hMFZleg,'FontSize',12,'Location','Northwest')

    set(hMFZ,'Color',[0.55 0.55 0.55],'LineWidth',0.5);
%     allNMF = load(sprintf('netmassflux_%s.txt',run));
%     avgNMF = mean(allNMF(:,2:end),1);
    hNMF = plot(avgNMF,1:JMAX-ghostcells,'-.','Color',[0 0.4 0.7],'LineWidth',3);
    set(hNMF,'DisplayName','Time-averaged profile')
    set(hMFZ,'DisplayName','Individual timestep profiles')
    title(axAvgMFZ,sprintf('Time-averaged solid mass flux\n%s',str))
    hMFZleg = legend(axAvgMFZ,[hNMF hMFZ hJet]);
    set(hMFZleg,'FontSize',12,'Location','Northwest')
    saveas(figAvgMFZ,fullfile(savepath,sprintf('TimeAvgNetMF_%s.jpg',run)));
  
  % ===================================================================== %
  
  
  cd(postpath)
  

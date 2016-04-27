% summaryPlots scrapes all figure directories for keyParameter files and
% makes summary plots of different parameters.

postpath = '~/scratch/pulse-post';
% Path to directory containing summaryPlots.m

% Path to directory containing subdirectories with keyParameter files
  savepath = '~/data2/ProductionRuns_Storage/';
  numKeyParams = 16;      % number of elements in keyParameters_*.txt files
                          % that are written in pulse_post3D
                        
% Define a color cycle for plots
  myColorOrder = [0.25 0.35 0.80    % blue
                  0.00 0.75 0.30    % green
                  0.85 0.00 0.20    % red
                  0.65 0.00 0.75    % purple
                  0.90 0.50 0.00];  % orange
            
% Define mass ratio values for transitions between buoyant/partial/collapse
  collapse2partial = 0.04;
  partial2buoyant  = 0.02;
  
% Define scatter symbology for buoyant/partial/collapse
  plumeStatSymbols = 'osv';
  

%%% =================================================================== %%%

cd(savepath)

% Scrape main directory for subdirectories containing simulation figures
  allContents = dir('*_*');
  dirFlags = [allContents.isdir];
  dirList = allContents(dirFlags);      % limit list to directories only
  N = size(dirList,1);

% Initialize cell array to hold key parameters for all simulations
  allKeyParams = cell(N,numKeyParams+1);

% Loop through each subdirectory and load parameters into cell array
  for k = 1:N
      
      cd(savepath)
      dirName = dirList(k).name;
      cd(dirName)
      
    % Load keyParameters.txt, or move on to next directory if no keyParams
      if exist(sprintf('keyParameters_%s.txt',dirName),'file') ~= 2
          warning('Excluding %s: no keyParameters file exists.',dirName);
          continue
      else
          keyParams = importdata(sprintf('keyParameters_%s.txt',dirName));
      end
      
      allKeyParams{k,1} = k;
      allKeyParams{k,2} = keyParams.textdata{2,1};
      allKeyParams(k,3:end) = num2cell(keyParams.data);
            
  end
  
  cd(savepath)

% Find and delete empty rows from cell array (empty rows occur if a
% directory did not contain the keyParameters file)
  emptyRows = cellfun('isempty',allKeyParams);
  allKeyParams(all(emptyRows,2),:) = [];
 
% Create table of key parameters for all simulations
  KPT = cell2table(allKeyParams,'VariableNames',...
      {'RunID'          % Unique number for each simulation
       'Pulse'          % Is pulsing on (T) or off (F)
       'Frequency'      % Pulse frequency (Hz)
       'MinEPG'         % Minimum gas volume fraction of pulse
       'MaxEPG'         % Maximum gas volume fraction of pulse (pre-clip)
       'DiamS1'         % Diameter of solid phase 1
       'DiamS2'         % Diameter of solid phase 2
       'DiamS3'         % Diameter of solid phase 3
       'StokesS1'       % Pulsed Stokes number for solid phase 1
       'StokesS2'       % Pulsed Stokes number for solid phase 2
       'StokesS3'       % Pulsed Stokes number for solid phase 3
       'AvgEPG'         % Average gas volume fraction of flow
       'AvgMFR'         % Average solid mass flow rate
       'AvgMFlux'       % Average solid mass flux
       'AvgEntr'        % Spatiotemporally averaged plume entrainment
       'AvgJEntr'       % Spatiotemporally averaged jet entrainment
       'MassRatio'});   % Ratio of total collapsed to erupted mass
   
% Delete empty rows (from if directory didn't contain keyParameters file)
  blankKPT = ismissing(KPT);
  KPT = KPT(~any(blankKPT,2),:);

% Replace steady 'frequencies' with fake value for plotting
  steadyFreq = min(KPT.Frequency)/10;
  for i = 1:size(KPT,1)
      if strcmp(KPT.Pulse{i},'F') == 1
          KPT.Frequency(i) = steadyFreq;
      end
  end
  
% Define index values for pulse frequency (generate cmap and legend entries) 
  allFrequencies = unique(KPT.Frequency);
  FreqIndex = zeros(length(KPT.Frequency),1);
  hFreqLeg = repmat(scatter(0,0),1,length(allFrequencies));
  fakeFig = figure('Visible','off'); fakeAx = axes('Parent',fakeFig);
  for i = 1:length(allFrequencies)
      idxFreq = find(KPT.Frequency == allFrequencies(i));
      FreqIndex(idxFreq) = i;
      hold on
      colorRow = mod(i,size(myColorOrder,1));
        if colorRow == 0
            colorRow = size(myColorOrder,1);
        end
      hFreqLeg(i) = scatter(fakeAx,0,0,[],myColorOrder(colorRow,:),'o','filled',...
          'MarkerEdgeColor','none','DisplayName',sprintf('%g Hz',allFrequencies(i)));
  end
  hFreqLeg(1).DisplayName = 'Steady';
  IDX = array2table(FreqIndex,'VariableNames',{'FreqIndex'});
  
% Define plume status based on mass ratio relative to defined transitions
  plumeStatus = cell(size(KPT,1),1);
  for i = 1:size(KPT,1)
      if KPT.MassRatio(i) > collapse2partial
          plumeStatus{i} = 'Collapse';
      elseif KPT.MassRatio(i) < partial2buoyant
          plumeStatus{i} = 'Buoyant';
      else plumeStatus{i} = 'Partial';
      end
  end
  PST = cell2table(plumeStatus,'VariableNames',{'PlumeStatus'});
          
% Combine key parameters, frequency index, and plume status and save table
  KPT = [KPT IDX PST];
  writetable(KPT,fullfile(savepath,'SummaryTable.txt'),'WriteVariableNames',1);
  
  
%%% ---------------- CREATE PARAMETRIZED SUMMARY PLOTS ---------------- %%%
            % Symbolically distinguish different column states %
            
% Set figure default values
  set(0,'DefaultFigurePaperPositionMode','auto',...
        'DefaultFigureColor','w',...
        'DefaultFigureVisible','off',...
        'DefaultFigureColorMap',myColorOrder,...
        'DefaultFigureUnits','normalized',...
        'DefaultFigurePosition',[0 0 1 1],...
        'DefaultAxesBox','on',...
        'DefaultAxesTickDir','in',...
        'DefaultAxesFontSize',12,...
        'DefaultAxesFontWeight','bold',...
        'DefaultAxesXGrid','on',...
        'DefaultAxesYGrid','on',...
        'DefaultAxesZGrid','on',...
        'DefaultAxesColorOrder',myColorOrder);

    
% Particle pulsed Stokes numbers vs. average plume entrainment
  fSTOKESvAVGENTR = figure;
  axSvAE = axes('Parent',fSTOKESvAVGENTR);
  hold on
  ylim([0 0.2])
  axSvAE.XScale = 'log';
  hSvAE1 = scatter(KPT.StokesS1,-KPT.AvgEntr,[],'o','filled',...
      'DisplayName',sprintf('Solid1 (d=%g mm)',KPT.DiamS1(1)*1E3));
  hSvAE2 = scatter(KPT.StokesS2,-KPT.AvgEntr,[],'o','filled',...
      'DisplayName',sprintf('Solid2 (d=%g mm)',KPT.DiamS2(1)*1E3));
  hSvAE3 = scatter(KPT.StokesS3,-KPT.AvgEntr,[],'o','filled',...
      'DisplayName',sprintf('Solid3 (d=%g mm)',KPT.DiamS3(1)*1E3));
  xlabel('Pulsed Stokes number');
  ylabel('Average plume entrainment');
  legend(axSvAE,[hSvAE1 hSvAE2 hSvAE3]);
  saveas(fSTOKESvAVGENTR,fullfile(savepath,'Stokes_AvgPEntr.jpg'));

% Particle pulsed Stokes numbers vs. average jet entrainment
  set(hSvAE1,'YData',-KPT.AvgJEntr);
  set(hSvAE2,'YData',-KPT.AvgJEntr);
  set(hSvAE3,'YData',-KPT.AvgJEntr);
  ylabel('Average jet entrainment');
  ylim([0 0.1])
  saveas(fSTOKESvAVGENTR,fullfile(savepath,'Stokes_AvgJEntr.jpg'));

  
% Particle pulsed Stokes numbers vs. mass ratio
  fSTOKESvMRATIO = figure;
  axSvMR = axes('Parent',fSTOKESvMRATIO);
  hold on
  ylim([0 0.5])
  axSvMR.XScale = 'log';
  hSvMR1 = scatter(KPT.StokesS1,KPT.MassRatio,[],'o','filled',...
      'DisplayName',sprintf('Solid1 (d=%g mm)',KPT.DiamS1(1)*1E3));
  hSvMR2 = scatter(KPT.StokesS2,KPT.MassRatio,[],'o','filled',...
      'DisplayName',sprintf('Solid2 (d=%g mm)',KPT.DiamS2(1)*1E3));
  hSvMR3 = scatter(KPT.StokesS3,KPT.MassRatio,[],'o','filled',...
      'DisplayName',sprintf('Solid3 (d=%g mm)',KPT.DiamS3(1)*1E3));
  xlabel('Pulsed Stokes number');
  ylabel('Ratio of collapsed to erupted mass');
  legend(axSvMR,[hSvMR1 hSvMR2 hSvMR3]);
  saveas(fSTOKESvMRATIO,fullfile(savepath,'Stokes_MassRatio.jpg'));

  
% Pulse frequency vs. average gas volume fraction
  fFREQvEPG = figure;
  axFvEPG = axes('Parent',fFREQvEPG);
  hold on
  axFvEPG.XTick = allFrequencies;
  axFvEPG.XScale = 'log';
  axFvEPG.XLim = [min(allFrequencies) max(allFrequencies)];
  axFvEPG.XTickLabel{1} = 'Steady';
  ylim([min(KPT.MinEPG) 1])
  hFvEPG = gscatter(KPT.Frequency,KPT.AvgEPG,KPT.PlumeStatus,...
      myColorOrder,plumeStatSymbols);
  for i = 1:length(hFvEPG)
      if strcmp(hFvEPG(i).DisplayName,'Collapse') == 1
          set(hFvEPG(i),'MarkerFaceColor',myColorOrder(1,:));
      elseif strcmp(hFvEPG(i).DisplayName,'Partial') == 1
          set(hFvEPG(i),'MarkerFaceColor',myColorOrder(2,:));
      elseif strcmp(hFvEPG(i).DisplayName,'Buoyant') == 1
          set(hFvEPG(i),'MarkerFaceColor',myColorOrder(3,:));
      end
  end
  set(hFvEPG,'MarkerEdgeColor','none');
  xlabel('Pulse frequency (Hz)');
  ylabel('Average gas volume fraction');
  saveas(fFREQvEPG,fullfile(savepath,'Freq_AvgEPG.jpg'));

  
% Pulse frequency vs. average mass flux
  fFREQvMFLUX = figure;
  axFvMF = axes('Parent',fFREQvMFLUX);
  hold on
  axFvMF.XTick = allFrequencies;
  axFvMF.XScale = 'log';
  axFvMF.XLim = [min(allFrequencies) max(allFrequencies)];
  axFvMF.XTickLabel{1} = 'Steady';
  hFvMF = gscatter(KPT.Frequency,KPT.AvgMFlux,KPT.PlumeStatus,...
      myColorOrder,plumeStatSymbols);
  for i = 1:length(hFvMF)
      if strcmp(hFvMF(i).DisplayName,'Collapse') == 1
          set(hFvMF(i),'MarkerFaceColor',myColorOrder(1,:));
      elseif strcmp(hFvMF(i).DisplayName,'Partial') == 1
          set(hFvMF(i),'MarkerFaceColor',myColorOrder(2,:));
      elseif strcmp(hFvMF(i).DisplayName,'Buoyant') == 1
          set(hFvMF(i),'MarkerFaceColor',myColorOrder(3,:));
      end
  end
  set(hFvMF,'MarkerEdgeColor','none');
  xlabel('Pulse frequency (Hz)');
  ylabel('Average mass flux (kg/m^2s)');
  saveas(fFREQvMFLUX,fullfile(savepath,'Freq_AvgMFlux.jpg'));
  
  
% Pulse frequency vs. mass ratio
  fFREQvMRATIO = figure;
  axFvMR = axes('Parent',fFREQvMRATIO);
  hold on
  ylim([0 0.5])
  axFvMR.XTick = allFrequencies;
  axFvMR.XScale = 'log';
  axFvMR.XLim = [min(allFrequencies) max(allFrequencies)];
  axFvMR.XTickLabel{1} = 'Steady';
  hFvMR = gscatter(KPT.Frequency,KPT.AvgMFlux,KPT.PlumeStatus,...
      myColorOrder,plumeStatSymbols);
  for i = 1:length(hFvMR)
      if strcmp(hFvMR(i).DisplayName,'Collapse') == 1
          set(hFvMR(i),'MarkerFaceColor',myColorOrder(1,:));
      elseif strcmp(hFvMR(i).DisplayName,'Partial') == 1
          set(hFvMR(i),'MarkerFaceColor',myColorOrder(2,:));
      elseif strcmp(hFvMR(i).DisplayName,'Buoyant') == 1
          set(hFvMR(i),'MarkerFaceColor',myColorOrder(3,:));
      end
  end
  set(hFvMR,'MarkerEdgeColor','none');
  xlabel('Pulse frequency (Hz)');
  ylabel('Ratio of collapsed to erupted mass');
  saveas(fFREQvMRATIO,fullfile(savepath,'Freq_MassRatio.jpg'));
  
  
% Average mass flux vs. mass ratio
  fMFLUXvMRATIO = figure;
  axMFvMR = axes('Parent',fMFLUXvMRATIO);
  hold on
  ylim([0 0.5])
  hMFvMR = scatter(KPT.AvgMFlux,KPT.MassRatio,[],KPT.FreqIndex,'o','filled');
  xlabel('Average mass flux (kg/m^2s)')
  ylabel('Ratio of collapsed to erupted mass')
  legMFvMR = legend(axMFvMR,hFreqLeg);
  saveas(fMFLUXvMRATIO,fullfile(savepath,'AvgMFlux_MassRatio.jpg'));


% Average mass flux vs. average plume entrainment
  fMFLUXvAVGENTR = figure;
  axMFvAE = axes('Parent',fMFLUXvAVGENTR);
  hold on
  ylim([0 0.2])
  hMFvAE = scatter(KPT.AvgMFlux,-KPT.AvgEntr,[],KPT.FreqIndex,'o','filled');
  xlabel('Average mass flux (kg/m^2s)');
  ylabel('Average plume entrainment');
  legMFvAE = legend(axMFvAE,hFreqLeg);
  saveas(fMFLUXvAVGENTR,fullfile(savepath,'AvgMFlux_AvgPEntr.jpg'));

% Average jet entrainment vs. average mass flux
  set(hMFvAE,'YData',-KPT.AvgJEntr);
  ylabel('Average jet entrainment');
  ylim([0 0.1])
  saveas(fMFLUXvAVGENTR,fullfile(savepath,'AvgMFlux_AvgJEntr.jpg'));

%%% =================================================================== %%%

close all
cd(postpath)
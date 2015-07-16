% test loop for variable number of isosurfaces

EPG = zeros([100,100,100]);
    EPG(1:10,:,:) = 0.1;
    EPG(11:20,:,:) = 0.2;
    EPG(21:30,:,:) = 0.3;
    EPG(31:40,:,:) = 0.4;
    EPG(41:50,:,:) = 0.5;
    EPG(51:60,:,:) = 0.6;
    EPG(61:70,:,:) = 0.7;
    EPG(71:80,:,:) = 0.8;
    EPG(81:90,:,:) = 0.9;
    EPG(91:100,:,:) = 1.0;
    
isoEPG = [0.1,0.3,0.5,0.7 0.9];
colEPG = [0.1 0.5 0.1;
          1.0 0.1 0.4;
          0.9 0.1 0.5;
          0.0 0.8 0.2
          0.1 0.6 0.3];
trnEPG = [0.1,0.2,0.3,0.5 0.7];

fig = figure('Name','testing3D_newgraphics','visible','on');
    hold on
    view(3)
    axis equal
    axis([3,102,3,102,3,102])
    grid on
    box on
    
% test initialize legend
names = cell(1,length(isoEPG));
for j = 1:length(isoEPG)
    names{j} = sprintf('EPG = %d',isoEPG(j));
end

% add variable number of isosurfaces
for i = 1:length(isoEPG)
    surf(i) = patch(isosurface(EPG,isoEPG(i)));
    set(surf(i),'FaceColor',colEPG(i,:),'EdgeColor','none','FaceAlpha',trnEPG(i))
    hold on
end
legend(names)
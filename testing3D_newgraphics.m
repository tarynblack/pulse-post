%%% Testing 3D data processing with new MATLAB graphics

run = 448706;
cd(sprintf('/Users/Taryn/Documents/MATLAB/MFIX_temp/%d',run))

EP_G = importdata(sprintf('EP_G_%d',run));

IMAX = 104;
JMAX = 104;
KMAX = 104;

timesteps = length(EP_G)/(IMAX*JMAX*KMAX);

fig = figure('Name','testing3D_newgraphics','visible','on');
hold on
view(3)
axis equal
axis([3,102,3,102,3,102])
grid on
box on


%%% testing all in one time loop
vidObj = VideoWriter('EPGtest.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 10;
open(vidObj);
set(gcf,'Visible','off');

for t = 50;%1:timesteps
    
    EPG = reshape(EP_G((t-1)*IMAX*JMAX*KMAX+1:t*IMAX*JMAX*KMAX),[JMAX IMAX KMAX]);
    EPG = EPG(3:end-2,3:end-2,3:end-2);
    EPG = permute(EPG,[3 2 1]);
        
    cla;
    
    surf1 = patch(isosurface(EPG,0.99999));
        set(surf1,'FaceColor','red','EdgeColor','none','FaceAlpha',0.1)  
    surf2 = patch(isosurface(EPG,0.99995));
        set(surf2,'FaceColor','yellow','EdgeColor','none','FaceAlpha',0.2)
    surf3 = patch(isosurface(EPG,0.9999));
        set(surf3,'FaceColor','cyan','EdgeColor','none','FaceAlpha',0.3)
    surf4 = patch(isosurface(EPG,0.9995));
        set(surf4,'FaceColor','black','EdgeColor','none','FaceAlpha',0.5)
    
    legend([surf1,surf2,surf3,surf4],'EPG=0.99999','EPG=0.99995','EPG=0.9999','EPG=0.9995')
    title(sprintf('Gas volume fraction, timestep=%d',t));
    
    vidfig = 'current.jpg';
    saveas(fig,vidfig);
%     newfig = sprintf('EPG_%dt',t);
%     savefig(fig,newfig,'compact');

    img = imread(vidfig);
    writeVideo(vidObj,img);
    
end

close(vidObj);
    
    
    

% for t = 1:timesteps
% 
%     EPG = reshape(EP_G((t-1)*IMAX*JMAX*KMAX+1:t*IMAX*JMAX*KMAX),[JMAX IMAX KMAX]);
%     EPG = EPG(3:end-2,3:end-2,3:end-2);
%     EPG = permute(EPG,[3 2 1]);
%     
%     filename = sprintf('EPGtest%d.jpg',t);
%     
%     cla;
%     
%     surf1 = patch(isosurface(EPG,0.99999));
%         set(surf1,'FaceColor','red','EdgeColor','none','FaceAlpha',0.1)
% %     surf2 = patch(isosurface(EPG,0.99995));
% %         set(surf2,'FaceColor','yellow','EdgeColor','none','FaceAlpha',0.2)
% %     surf3 = patch(isosurface(EPG,0.9999));
% %         set(surf3,'FaceColor','cyan','EdgeColor','none','FaceAlpha',0.3)
% %     surf4 = patch(isosurface(EPG,0.9995));
% %         set(surf4,'FaceColor','black','EdgeColor','none','FaceAlpha',0.5)
% %     legend([surf1,surf2,surf3,surf4],'EPG=0.99999','EPG=0.99995','EPG=0.9999','EPG=0.9995')
% 
%     title(sprintf('Gas volume fraction, timestep=%d',t));
% 
% %     saveas(fig,filename,'jpg');
% 
% end
% 
% vidObj = VideoWriter('EPGtest.avi');
% vidObj.Quality = 100;
% vidObj.FrameRate = 10;
% open(vidObj);
% set(gcf,'Visible','off');
% for t = 1:timesteps
%     cla;
%     filename = ['EPGtest',num2str(t),'.jpg'];
%     img = imread(filename);
%     writeVideo(vidObj,img);
% end
% close(vidObj);
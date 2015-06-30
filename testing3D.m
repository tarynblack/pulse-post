%%% Testing 3D data - loading, manipulating, displaying

run = 366892;

cd(sprintf('%d',run))

EP_G = importdata(sprintf('EP_G_%d',run));

IMAX = 104;
JMAX = 104;
KMAX = 104;

timesteps = length(EP_G)/(IMAX*JMAX*KMAX);

fig = figure('Name','testing_3D','units','normalized','outerposition',[0 0 1 1]);
m = moviein(timesteps,fig);
hold on
view(3)
axis equal
axis([3,102,3,102,3,102])
grid on
box on

surf1 = patch(isosurface(zeros(2,2,2),0));
    set(surf1,'FaceColor','red','EdgeColor','none','FaceAlpha',0.1)
surf2 = patch(isosurface(zeros(2,2,2),0));
    set(surf2,'FaceColor','yellow','EdgeColor','none','FaceAlpha',0.2)
surf3 = patch(isosurface(zeros(2,2,2),0));
    set(surf3,'FaceColor','cyan','EdgeColor','none','FaceAlpha',0.3)
surf4 = patch(isosurface(zeros(2,2,2),0));
    set(surf4,'FaceColor','black','EdgeColor','none','FaceAlpha',0.5)
legend([surf1,surf2,surf3,surf4],'EPG=0.99999','EPG=0.99995','EPG=0.9999','EPG=0.9995')


for t = 1:timesteps

    EPG = reshape(EP_G((t-1)*IMAX*JMAX*KMAX+1:t*IMAX*JMAX*KMAX),[JMAX IMAX KMAX]);
    EPG = EPG(3:end-2,3:end-2,3:end-2);
    EPG = permute(EPG,[3 2 1]);

%     surf1 = patch(isosurface(EPG,0.99999));
%     set(surf1,'FaceColor',[0.6 0.6 0.6],'EdgeColor','none','FaceAlpha',0.1)
%     surf2 = patch(isosurface(EPG,0.99995));
%     set(surf2,'FaceColor',[0.6 0.6 0.6],'EdgeColor','none','FaceAlpha',0.3)
%     surf3 = patch(isosurface(EPG,0.9999));
%     set(surf3,'FaceColor',[0.6 0.6 0.6],'EdgeColor','none','FaceAlpha',0.5)
%     surf4 = patch(isosurface(EPG,0.9995));
%     set(surf4,'FaceColor',[0.6 0.6 0.6],'EdgeColor','none','FaceAlpha',0.7)

    surf1 = patch(isosurface(EPG,0.99999));
        set(surf1,'FaceColor','red','EdgeColor','none','FaceAlpha',0.1)
    surf2 = patch(isosurface(EPG,0.99995));
        set(surf2,'FaceColor','yellow','EdgeColor','none','FaceAlpha',0.2)
    surf3 = patch(isosurface(EPG,0.9999));
        set(surf3,'FaceColor','cyan','EdgeColor','none','FaceAlpha',0.3)
    surf4 = patch(isosurface(EPG,0.9995));
        set(surf4,'FaceColor','black','EdgeColor','none','FaceAlpha',0.5)
    
%     refreshdata(surf1)
%     refreshdata(surf2)
%     refreshdata(surf3)
%     refreshdata(surf4)
        
    title(sprintf('Gas volume fraction, timestep=%d',t));

    m(:,t) = getframe(gcf);
    
    cla
    
%     clear('surf1','surf2','surf3','surf4')
%     clear surf2
%     clear surf3 
%     clear surf4

%     delete(surf1,surf2,surf3,surf4)

end

movie2avi(m,sprintf('testing_3D_%dsteady.avi',run));

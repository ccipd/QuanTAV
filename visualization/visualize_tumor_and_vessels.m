
function visualize_tumor_and_vessels(seg_v,seg_t)

figure;
hold on;

view(3)
camlight
hold on;
FV = isosurface(seg_t,0.5);
patch(FV,'facecolor',[246/255 210/255 140/255],'facealpha',.6,'edgecolor','none');
set(gcf,'units','normalized','outerposition',[0 0 1 1])

FV = isosurface(seg_v,.5);
patch(FV,'facevertexcdata',[215/255 71/255 31/255],'facecolor','flat','facealpha',.65,'edgecolor','none');

ax = gcf;
axis equal;
axis tight

lighting gouraud     
axis on
grid on
set(gca,'Yticklabel',[])
set(gca,'Xticklabel',[])
set(gca,'Zticklabel',[])
box on
set(gcf,'InvertHardCopy','off')


set(gca,'GridAlpha',0)
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
set(gca,'Color','k')
set(gcf,'color','k')
axis off
grid off



function visualize_quantav_morphology(seg_v,seg_t, branches,curv_range,mode)

if ~exist('mode','var')
    mode = 'all'; 
end
disp(mode)

%% CURVATURE
if strcmp(mode,'curvature') || strcmp(mode,'all')
    lineStyles = linspecer(200,'sequential');

    figure;
    hold on;

    for n = 1:length(branches)
        branch_info = branches{n}; 
        L = branch_info.points; 
        curv = branch_info.curvature; 
        for i=2:size(L,1)-1
            curv_norm = max([min([(curv(i)-curv_range(1))/(curv_range(2)-curv_range(1)),1]),0]);
            color = lineStyles(1+round(199*curv_norm),:);
            linewidth = 2+1.5*curv_norm;
            plot3(L(i-1:i+1,2),L(i-1:i+1,1),L(i-1:i+1,3),'-','Color', color,'LineWidth',linewidth);
        end
    end

    view(3)
    camlight
    hold on;
    FV = isosurface(seg_t,0.5);
    patch(FV,'facecolor',[246/255 210/255 140/255],'facealpha',.6,'edgecolor','none');
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
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
end

%% TORSION 
if strcmp(mode,'torsion') || strcmp(mode,'all')

    tors = [];
    pts = [];
    figure;
    hold on;
    for n = 1:length(branches)
        branch_info = branches{n}; 
        L = branch_info.points; 
        pts = [pts;L];
        tors = [tors;branch_info.torsion*ones(size(L,1),1)]; 
        plot3(L(:,2),L(:,1),L(:,3),'-','Color',[1,1,1],'LineWidth',2);

    end

    [X,Y,Z] = ind2sub(size(seg_v),find(seg_v(:)));
    vpts = [X,Y,Z];

    [D,I] = pdist2(pts,vpts,'squaredeuclidean','Smallest',1);
    ID = D.^-1;
    interpolated_torsion=tors(I);%sum(tors(I).*ID,2)./sum(ID,2);

    color_vol = zeros(size(seg_t)); 
    color_vol(seg_v(:)) = interpolated_torsion; 
    max_val = .75;
    color_vol(color_vol>max_val) = max_val;

    FV = isosurface(seg_v,.5,color_vol);

    patch(FV,'facevertexcdata',FV.facevertexcdata,'facecolor','flat','facealpha',.65,'edgecolor','none');
    %     isocolors(X,Y,Z,color_vol,p)
    view(3)
    camlight
    % Display the skeleton
    hold on;
    FV = isosurface(seg_t,0.5);
    patch(FV,'facecolor',[246/255 210/255 140/255],'facealpha',.8,'edgecolor','none');
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
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

    %      set(gcf,'Color','k')
    colormap(inferno)
    camlight
%     campos([1810.80660060858,-1523.81180211257,1146.07146697799]);
    colorbar
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

end


function masked_crops = my_masked_crops(structure_name, fit_again_later,starting_pt,stopping_pt)
masked_crops = {};
for i = starting_pt:stopping_pt
    close all
    j = fit_again_later(i); 
    
    
    figure(1);
    set(gcf,'Position',[958.3333,   36.6667,  316.3333,  608.6667]);
    title(structure_name(j).name);
    s1=subplot(1,3,1);imagesc(structure_name(j).crop),colormap(gray);
    title('kymograph');
    xlabel('x, px');
    ylabel('y, px');
    caxis([0 8]);
    s2=subplot(1,3,2);imagesc(structure_name(j).crop),colormap(gray);
    caxis([0 8]);
    hold on
    plot(s2,structure_name(j).pixel,structure_name(j).timepix,'-r', 'LineWidth', 1);
    title('overlay');
    xlabel('x, px');
    ylabel('y, px');
    hold off
    s3=subplot(1,3,3);
    plot(s3,structure_name(j).pixel,structure_name(j).timepix,'ro');
    s3=gca;
    sz=size(structure_name(j).crop);
    t=sz(:,1);
    w=sz(:,2);
    xlim(gca,[0 w]);
    ylim(gca,[0 t]);
    set(gca,'YDir','Reverse');
    title('trace');
    xlabel('x, px');
    ylabel('y, px');
    
    
    figure(2)
    crop = structure_name(j).crop;
    imagesc(crop);
    p = drawpolygon();
    xy = round(p.Position);
    x = xy(:,1);
    y = xy(:,2);
    clf
    bw = poly2mask(x,y,size(crop, 1),size(crop,2));
    figure
    masked_crops{j} = crop.*bw;
    imagesc(masked_crops{j})
end
end 
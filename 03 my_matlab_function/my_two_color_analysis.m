function instances_of_crossover = my_two_color_analysis(Green_MY_DATA, Red_MY_DATA)
for i = 1:length(Green_MY_DATA)
    Green_name_contents{i} = strsplit(Green_MY_DATA(i).name);
end 

for i = 1:length(Green_name_contents)
    Green_pull_these_kymos{i}{1} = Green_name_contents{1,i}{1};
    Green_pull_these_kymos{i}{2} = join([Green_name_contents{1,i}{2},Green_name_contents{1,i}{3}]);
end

for i = 1:length(Red_MY_DATA)
    Red_name_contents{i} = strsplit(Red_MY_DATA(i).name);
end 

for i = 1:length(Red_name_contents)
    Red_pull_these_kymos{i}{1} = Red_name_contents{1,i}{1};
    Red_pull_these_kymos{i}{2} = join([Red_name_contents{1,i}{2},Red_name_contents{1,i}{3}]);
end
for i = start_var:end_var
    filename = pre_MSD(i).name;
    timestep = pre_MSD(i).line_time;
    intensity = transpose(pre_MSD(i).intensity);
    time= pre_MSD(i).timepix*timestep;
    time=transpose(time)/1000;
    position= pre_MSD(i).pixel.*0.1;
    position=transpose(position);
    figure(1);
    yyaxis left
    plot(time, position)
    ylabel('position, \mum');
    yyaxis right
    plot(time, intensity)
    set(gcf,'Position',[0,35,480,300]);% for laptop
    title(pre_MSD(i).name, 'Interpreter', 'none');
    xlabel('time, s');
    xlim([0 max(time)]);
    ylabel('intensity, a.u.');
    figure(2);
    set(gcf,'Position',[482.0000,   31.0000,  472.6667,  615.6667]);
    imagesc(rot90(pre_MSD(i).full_kymo,3));
    hold on
    rectangle('Position',[pre_MSD(i).crop_coordinates(3),pre_MSD(i).crop_coordinates(1),...
        pre_MSD(i).crop_coordinates(4)-pre_MSD(i).crop_coordinates(3),pre_MSD(i).crop_coordinates(2)-...
        pre_MSD(i).crop_coordinates(1)], 'EdgeColor','y'); 
    hold off
    figure(3);
    set(gcf,'Position',[958.3333,   36.6667,  316.3333,  608.6667]);
    title(pre_MSD(i).name);
    s1=subplot(1,3,1);imagesc(pre_MSD(i).crop),colormap(gray);
    title('kymograph');
    xlabel('x, px');
    ylabel('y, px');
    s2=subplot(1,3,2);imagesc(pre_MSD(i).crop),colormap(gray);
    hold on
    plot(s2,pre_MSD(i).pixel,pre_MSD(i).timepix,'-y');
    title('overlay');
    xlabel('x, px');
    ylabel('y, px');
    hold off
    s3=subplot(1,3,3);
    plot(s3,pre_MSD(i).pixel,pre_MSD(i).timepix,'ro');
    s3=gca;
    sz=size(pre_MSD(i).crop);
    t=sz(:,1);
    w=sz(:,2);
    xlim(gca,[0 w]);
    ylim(gca,[0 t]);
    set(gca,'YDir','Reverse');
    title('trace');
    xlabel('x, px');
    ylabel('y, px');
    dispthis = ['you are on trace ', num2str(i), ' out of ', num2str(end_var)];
    choice = menu(dispthis,'Good Fit','Bad Fit','truncate_top','truncate_bottom','truncate_between(top-bottom)', 'Break');
    disp(['you are on trace ', num2str(i), ' out of ', num2str(end_var)]);
    if choice == 1
    keep_these(i).name = pre_MSD(i).name;
    keep_these(i).particle_tracked = horzcat(pre_MSD(i).timepix',pre_MSD(i).pixel',pre_MSD(i).intensity');
    keep_these(i).crop_coordinates = pre_MSD(i).crop_coordinates;
    keep_these(i).crop = pre_MSD(i).crop;
    keep_these(i).line_time = pre_MSD(i).line_time;
    keep_these(i).kymograph = pre_MSD(i).full_kymo;
    keep_these(i).truncated = 0;
    elseif choice == 2
    keep_these(i).name = [];
    keep_these(i).particle_tracked = [];
    keep_these(i).crop_coordinates = [];
    keep_these(i).crop = [];
    keep_these(i).line_time = [];
    keep_these(i).kymograph = [];
    elseif choice == 3
    [x,y] = ginput(1);
    y = round(y);
    keep_these(i).name = pre_MSD(i).name;
    keep_these(i).particle_tracked = horzcat(pre_MSD(i).timepix(1:end-y+1)',pre_MSD(i).pixel(y:end)',pre_MSD(i).intensity(y:end)');
    keep_these(i).crop_coordinates = pre_MSD(i).crop_coordinates + [y, 0, 0, 0];
    keep_these(i).crop = pre_MSD(i).crop(y:end,:);
    keep_these(i).line_time = pre_MSD(i).line_time;
    keep_these(i).kymograph = pre_MSD(i).full_kymo;
    keep_these(i).truncated = 1;
    elseif choice == 4
    [x,y] = ginput(1);
    y = round(y);
    keep_these(i).name = pre_MSD(i).name;
    keep_these(i).particle_tracked = horzcat(pre_MSD(i).timepix(1:y)',pre_MSD(i).pixel(1:y)',pre_MSD(i).intensity(1:y)');
    keep_these(i).crop_coordinates = pre_MSD(i).crop_coordinates - [0, y, 0, 0];
    keep_these(i).crop = pre_MSD(i).crop(1:y,:);
    keep_these(i).line_time = pre_MSD(i).line_time;
    keep_these(i).kymograph = pre_MSD(i).full_kymo;
    keep_these(i).truncated = 1;
       elseif choice == 5
    [x,y] = ginput(2);
    y = round(y);
    keep_these(i).name = pre_MSD(i).name;
    keep_these(i).particle_tracked = horzcat(pre_MSD(i).timepix(1:y(2)-y(1)+1)',pre_MSD(i).pixel(y(1):y(2))',pre_MSD(i).intensity(y(1):y(2))');
    keep_these(i).crop_coordinates = pre_MSD(i).crop_coordinates + [y(1), -y(2), 0, 0];
    keep_these(i).crop = pre_MSD(i).crop(y(1):y(2),:);
    keep_these(i).line_time = pre_MSD(i).line_time;
    keep_these(i).kymograph = pre_MSD(i).full_kymo;
    keep_these(i).truncated = 1;
    elseif choice == 6
        break
    end
end
end
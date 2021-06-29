%% TWO COLOR CROSS-OVER ANALYSIS
%Inserts a column in the Green_My_Data with the corresponding location of
%dCas9 kymo information

for i = 1:length(Green_MY_DATA)
    Green_name_contents{i} = strsplit(Green_MY_DATA(i).name);
end 

for i = 1:length(Green_name_contents)
    Green_pull_these_kymos{i}{1} = [Green_name_contents{1,i}{1},' '];
    Green_kymonames{i} = [Green_name_contents{1,i}{1},' '];
    temp = strsplit(Green_name_contents{1,i}{3}, '_');
    Green_pull_these_kymos{i}{2} = [Green_name_contents{1,i}{2},' ',temp{1}];
end
Red_names = {Red_MY_DATA.name}.';
for i = 1:length(Green_pull_these_kymos)
    x = 1;
    for j = 1:length(Red_names)
        if contains(Red_names{j}, Green_pull_these_kymos{i}{1}) && contains(Red_names{j}, Green_pull_these_kymos{i}{2})
            Green_MY_DATA(i).Red_Info_Loc(x) = j;
            x=x+1;
        end
    end
end
for i = 1:length(Green_MY_DATA)
    if isempty(Green_MY_DATA(i).Red_Info_Loc)
        for k = 1:length(Green_pull_these_kymos)
            y= 1;
            for j = 1:length(Red_names)
                if contains(Red_names{j}, Green_pull_these_kymos{k}{2}) && not(contains(Red_names{j}, Green_pull_these_kymos{k}{1}))
                    Green_MY_DATA(k).Red_Info_Loc(y) = j;
                    Green_MY_DATA(k).not_original_red = 1;
                    y= y+1;
                end
            end
        end
    end
end

% summarize the kymo numbers
uniq_kymo_names = unique(Green_kymonames);

for i = 1:length(uniq_kymo_names)
    for j = 1:length(Green_MY_DATA)
        a =  Green_pull_these_kymos{j}{1};
        if strcmp(a, uniq_kymo_names(i))
            Green_MY_DATA(j).Kymo_Number = i;
        end
    end
end      


% %% Side by side red (initial) and green (full)
% for i = 2
%     fh1 = figure;
%     sfh1 = subplot(1,2,1,'Parent',fh1);
%     imagesc(flip(Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(1)).kymograph(:,1:80)))
% %     ax = gca;
%     colormap(sfh1,red_map)
%     sfh1.XTickLabel = [];
%     sfh1.YTickLabel = [];
%     % ax.DataAspectRatio = [1 1 1];
%     sfh1.Position = sfh1.Position-[0 0 0.25 0];
%     caxis([0 8]) 
%     sfh2 = subplot(1,2,2,'Parent',fh1);
% 
%     imagesc(flip(Green_MY_DATA(i).kymograph))
% %     ax = gca;
%     colormap(sfh2,green_map)
%     sfh2.XTickLabel = [];
%     sfh2.YTickLabel = [];
%     % ax.DataAspectRatio = [1 1 1];
%     sfh2.Position = sfh2.Position-[0.35 0 0 0];
%     sfh2.Position = sfh2.Position+[0 0 0.4 0];
%     caxis([0 8])
% end


% %% Extend the red line into the green figure 
% 
% for lmnop = 2%length(uniq_kymo_names)
%     if isempty(Green_MY_DATA(lmnop).Red_Info_Loc)
%         continue
%     end
%     ii = find([Green_MY_DATA.Kymo_Number] == lmnop);
%     i = ii(1);
%     threshold_for_red_trace_in_subplot = 80;
%     num_of_red_traces = length(Green_MY_DATA(i).Red_Info_Loc);
%     fh1 = figure;%('visible','off');% select to suppress figures when saving
%     sfh1 = subplot(1,2,1,'Parent',fh1);
%     imagesc(flip(Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(1)).kymograph(:,1:threshold_for_red_trace_in_subplot)))
%     colormap(sfh1,red_map)
%     sfh1.XTickLabel = [];
%     sfh1.YTickLabel = [];
%     sfh1.Position = sfh1.Position-[0 0 0 0];
%     caxis([0 8]) 
%     
%     hold on
%     for j = 1:num_of_red_traces
%         red_chords = Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(j)).crop_coordinates;
%         red_position = Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(j)).particle_tracked;
%         length_o_trace = length(red_position);
%         if length_o_trace < threshold_for_red_trace_in_subplot
%             plot(red_position(1:length_o_trace,1)+red_chords(1)-1,red_position(1:length_o_trace,2)+red_chords(3)-1,'Color','g','LineWidth',2)
%         else
%             plot(red_position(1:threshold_for_red_trace_in_subplot,1)+red_chords(1)-1,red_position(1:threshold_for_red_trace_in_subplot,2)+red_chords(3)-1,'Color','g','LineWidth',2)
%         end
%     end
%     hold off
%     
%     sfh2 = subplot(1,2,2,'Parent',fh1);
%     imagesc(flip(Green_MY_DATA(i).kymograph))
%     colormap(sfh2,green_map)
%     sfh2.XTickLabel = [];
%     sfh2.YTickLabel = [];
%     sfh2.Position = sfh2.Position-[0.35 0 0 0];
%     sfh2.Position = sfh2.Position+[0 0 0.4 0];
%     caxis([0 8])
%     
%     hold on
%     
% %     red_chords = Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(1)).crop_coordinates;
% % %     holdingplace1 = red_chords(4)-red_chords(3);
% % %     holdingplace2 = red_chords(2)- red_chords(1);
% % %     rectangle('Position',[red_chords(1),red_chords(3), holdingplace2 , holdingplace1], 'EdgeColor','r'); 
% %     
% %     red_position = Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(1)).particle_tracked;
% %     plot(red_position(:,1)+red_chords(1)-1,red_position(:,2)+red_chords(3)-1,'Color','r','LineWidth',2)
% for j = 1:num_of_red_traces
%     red_chords = Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(j)).crop_coordinates;
%     red_position = Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(j)).particle_tracked;
%     red_avg_position = mean(red_position(:,2));
%     red_avg_pos_for_full_kymo = red_avg_position+red_chords(3)-1;
%     length_green_trace = length(Green_MY_DATA(i).particle_tracked);
%     yline(red_avg_pos_for_full_kymo,'Color','r','LineWidth',2)
% end
%     hold off
% %     saveas(fh1,['two_color ',uniq_kymo_names{lmnop},'.png'],'png')
% 
% end
%% Plotting the SWR1 trajectories with Cas9 trajectories 

% green map
green_map = zeros(3);
for i = 1:11
green_map(i,2) = i*0.1-0.1;
end
% red map
red_map = zeros(3);
for i = 1:11
red_map(i,1) = i*0.1-0.1;
end

for lmnop = 1:length(uniq_kymo_names)
    if isempty(Green_MY_DATA(lmnop).Red_Info_Loc)
        continue
    end
    ii = find([Green_MY_DATA.Kymo_Number] == lmnop);
    i = ii(1);
    threshold_for_red_trace_in_subplot = 80;
    num_of_red_traces = length(Green_MY_DATA(i).Red_Info_Loc);
    fh1 = figure('visible','off');% select to suppress figures when saving
%     fh1 = figure;
    sfh1 = subplot(1,2,1,'Parent',fh1);
    imagesc(flip(Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(1)).kymograph(:,1:threshold_for_red_trace_in_subplot)))
    colormap(sfh1,red_map)
    sfh1.XTickLabel = [];
    sfh1.YTickLabel = [];
    sfh1.Position = sfh1.Position-[0 0 0 0];
    caxis([0 8]) 
    
    hold on
    for j = 1:num_of_red_traces
        red_chords = Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(j)).crop_coordinates;
        red_position = Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(j)).particle_tracked;
        length_o_trace = length(red_position);
        if length_o_trace < threshold_for_red_trace_in_subplot
            plot(red_position(1:length_o_trace,1)+red_chords(1)-1,red_position(1:length_o_trace,2)+red_chords(3)-1,'Color','g','LineWidth',2)
        else
            plot(red_position(1:threshold_for_red_trace_in_subplot,1)+red_chords(1)-1,red_position(1:threshold_for_red_trace_in_subplot,2)+red_chords(3)-1,'Color','g','LineWidth',2)
        end
    end
    hold off
    
    sfh2 = subplot(1,2,2,'Parent',fh1);
    imagesc(flip(Green_MY_DATA(i).kymograph))
    colormap(sfh2,green_map)
    sfh2.XTickLabel = [];
    sfh2.YTickLabel = [];
    sfh2.Position = sfh2.Position-[0.35 0 0 0];
    sfh2.Position = sfh2.Position+[0 0 0.4 0];
    caxis([0 8])
    
    hold on
    
    
 
    
    %     red_chords = Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(1)).crop_coordinates;
% %     holdingplace1 = red_chords(4)-red_chords(3);
% %     holdingplace2 = red_chords(2)- red_chords(1);
% %     rectangle('Position',[red_chords(1),red_chords(3), holdingplace2 , holdingplace1], 'EdgeColor','r'); 
%     
%     red_position = Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(1)).particle_tracked;
%     plot(red_position(:,1)+red_chords(1)-1,red_position(:,2)+red_chords(3)-1,'Color','r','LineWidth',2)

%Plots dCas9
for j = 1:num_of_red_traces
    red_chords = Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(j)).crop_coordinates;
    red_position = Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(j)).particle_tracked;
    length_red_line = length(red_position);
    red_avg_position = mean(red_position(:,2));
    red_avg_pos_for_full_kymo = red_avg_position+red_chords(3)-1;
    length_green_trace = length(Green_MY_DATA(i).particle_tracked);
%     yline(red_avg_pos_for_full_kymo,'Color','r','LineWidth',2)
    if isempty(Green_MY_DATA(i).not_original_red)
        x_is = [0:length_red_line];
        red_avg_pos_for_full_kymo_plot = red_avg_pos_for_full_kymo.* ones(length(x_is),1);
        plot(x_is, red_avg_pos_for_full_kymo_plot,'Color','r','LineWidth',2)

        x_is = [length_red_line+1:length(Green_MY_DATA(i).kymograph)];
        red_avg_pos_for_full_kymo_plot = red_avg_pos_for_full_kymo.* ones(length(x_is),1);
        plot(x_is, red_avg_pos_for_full_kymo_plot,'Color','m','LineStyle', ':','LineWidth',2)
    else
        x_is = [0:length(Green_MY_DATA(i).kymograph)];
        red_avg_pos_for_full_kymo_plot = red_avg_pos_for_full_kymo.* ones(length(x_is),1);
        plot(x_is, red_avg_pos_for_full_kymo_plot,'Color','m','LineStyle', ':','LineWidth',2)
    end
end
    hold off
    
%Plots SWR1   
    hold on
    for j = 1:length(ii)
        green_chords = Green_MY_DATA(ii(j)).crop_coordinates;
        green_position = Green_MY_DATA(ii(j)).particle_tracked;
        plot(green_position(:,1)+green_chords(1)-1,green_position(:,2)+green_chords(3)-1,'Color','y','LineWidth',2)
    end
    hold off
    saveas(fh1,['two_color with positions after second pass magenta ',uniq_kymo_names{lmnop},'.png'],'png')

end


% %% random test site in script
% 
%     sfh1 = figure;
%     imagesc(Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(1)).crop)
%     colormap(sfh1,red_map)
% %     sfh1.XTickLabel = [];
% %     sfh1.YTickLabel = [];
%     caxis([0 8]) 
% 
%     hold on
%     red_position = Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(1)).particle_tracked;
%     plot(red_position(1:80,2),red_position(1:80,1),'Color','g','LineWidth',2)
%     hold off
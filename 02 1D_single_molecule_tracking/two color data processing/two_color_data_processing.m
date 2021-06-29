%% Plots
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
for i = 6%:length(correlated_struct)
    disp(correlated_struct(i).DNA)
    if ~isempty(correlated_struct(i).green_info) && ~isempty(correlated_struct(i).red_info)
        for j = 1:length(correlated_struct(i).green_info)
            if j == 1
                color_is = "r";
            else
                color_is = "c";
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sfh2 = figure(j);
            set(sfh2, "Visible", "on");
            hold on
            axis off
            imagesc(flip(correlated_struct(i).green_info(j).kymograph))
            length_of_kymo = length(correlated_struct(i).green_info(j).kymograph);
            colormap(sfh2,green_map)
            caxis([0 8])
            for k = 1:length(correlated_struct(i).green_info(j).particle_tracked)
                green_chords = correlated_struct(i).green_info(j).crop_coordinates{k};
                green_position = correlated_struct(i).green_info(j).particle_tracked{k};
                plot(green_position(:,1)+green_chords(1)-1,green_position(:,2)+green_chords(3)-1,'Color','m','LineWidth',1)
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k = 1:length(correlated_struct(i).red_info.particle_tracked)
                red_chords = correlated_struct(i).red_info.crop_coordinates{k};
                red_position = correlated_struct(i).red_info.particle_tracked{k};
                plot(red_position(:,1)+red_chords(1)-1,red_position(:,2)+red_chords(3)-1,'Color',color_is,'LineWidth',1)
                length_of_trace = length(red_position);
                red_avg_position = mean(red_position(:,2));
                red_avg_pos_for_full_kymo = red_avg_position+red_chords(3)-1;
                x_is = [length_of_trace+red_chords(1)-1:length_of_kymo];
                red_avg_pos_for_full_kymo_plot = red_avg_pos_for_full_kymo.* ones(length(x_is),1);
                plot(x_is, red_avg_pos_for_full_kymo_plot,'Color',color_is,'LineStyle', ':','LineWidth',1)
            end
            hold off
            disp(correlated_struct(i).green_info(j).kymoname)
%             saveas(sfh2, [correlated_struct(i).DNA, correlated_struct(i).green_info(j).kymoname], "png")
%             close all
        end 
        sfh1 = figure(i);
        set(sfh1, "Visible", "on");
        colormap(sfh1,red_map)
        caxis([0 8])
        hold on
        axis off
        imagesc(flip(correlated_struct(i).red_info.kymograph))
        length_of_kymo = length(correlated_struct(i).red_info.kymograph);
        color_is = "c";
        for k = 1:length(correlated_struct(i).red_info.particle_tracked)
                red_chords = correlated_struct(i).red_info.crop_coordinates{k};
                red_position = correlated_struct(i).red_info.particle_tracked{k};
                plot(red_position(:,1)+red_chords(1)-1,red_position(:,2)+red_chords(3)-1,'Color',color_is,'LineWidth',1)
                length_of_trace = length(red_position);
                red_avg_position = mean(red_position(:,2));
                red_avg_pos_for_full_kymo = red_avg_position+red_chords(3)-1;
                x_is = [length_of_trace+red_chords(1)-1:length_of_kymo];
                red_avg_pos_for_full_kymo_plot = red_avg_pos_for_full_kymo.* ones(length(x_is),1);
                plot(x_is, red_avg_pos_for_full_kymo_plot,'Color',color_is,'LineStyle', ':','LineWidth',1)
        end
%         saveas(sfh1, strcat("red only ", correlated_struct(i).DNA), "png")
%         close all
    end
end
%% Colocalizations
for i = 1:length(correlated_struct)
    disp(correlated_struct(i).DNA)
    if ~isempty(correlated_struct(i).green_info) && ~isempty(correlated_struct(i).red_info)
        for j = 1:length(correlated_struct(i).green_info)
            correlated_struct(i).correlated_data(j).kymoname =correlated_struct(i).green_info(j).kymoname; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            length_of_kymo = length(correlated_struct(i).green_info(j).kymograph);
            for k = 1:length(correlated_struct(i).green_info(j).particle_tracked)
                green_chords = correlated_struct(i).green_info(j).crop_coordinates{k};
                green_position = correlated_struct(i).green_info(j).particle_tracked{k};
%                 plot(green_position(:,1)+green_chords(1)-1,green_position(:,2)+green_chords(3)-1,'Color','m','LineWidth',1)
                correlated_struct(i).correlated_data(j).green_data{1, k}= green_position(:,1)+green_chords(1)-1;
                correlated_struct(i).correlated_data(j).green_data{2, k}= green_position(:,2)+green_chords(3)-1;
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k = 1:length(correlated_struct(i).red_info.particle_tracked)
                red_chords = correlated_struct(i).red_info.crop_coordinates{k};
                red_position = correlated_struct(i).red_info.particle_tracked{k};
%                 plot(red_position(:,1)+red_chords(1)-1,red_position(:,2)+red_chords(3)-1,'Color',color_is,'LineWidth',1)
                correlated_struct(i).correlated_data(j).red_data{1, k}= red_position(:,1)+red_chords(1)-1;
                correlated_struct(i).correlated_data(j).red_data{2, k}= red_position(:,2)+red_chords(3)-1;
                correlated_struct(i).correlated_data(j).red_data{5, k}= (red_position(:,1)+red_chords(1)-1)';
                correlated_struct(i).correlated_data(j).red_data{6, k}= (red_position(:,2)+red_chords(3)-1)';
                length_of_trace = length(red_position);
                red_avg_position = mean(red_position(:,2));
                red_avg_pos_for_full_kymo = red_avg_position+red_chords(3)-1;
                x_is = [length_of_trace+red_chords(1)-1:length_of_kymo];
                x_is_2 = [red_chords(1):length_of_kymo];
                red_avg_pos_for_full_kymo_plot = red_avg_pos_for_full_kymo.* ones(length(x_is),1);
                red_avg_pos_for_full_kymo_plot_2 = red_avg_pos_for_full_kymo.* ones(length(x_is_2),1);
%                 plot(x_is, red_avg_pos_for_full_kymo_plot,'Color',color_is,'LineStyle', ':','LineWidth',1)
                correlated_struct(i).correlated_data(j).red_data{3, k}= x_is;
                correlated_struct(i).correlated_data(j).red_data{4, k}= red_avg_pos_for_full_kymo_plot;
                correlated_struct(i).correlated_data(j).red_data{5, k}= [correlated_struct(i).correlated_data(j).red_data{5, k} x_is];
                correlated_struct(i).correlated_data(j).red_data{6, k}= [correlated_struct(i).correlated_data(j).red_data{6, k} red_avg_pos_for_full_kymo_plot'];
                correlated_struct(i).correlated_data(j).red_data{7, k}= x_is_2;
                correlated_struct(i).correlated_data(j).red_data{8, k}= red_avg_pos_for_full_kymo_plot_2';
            end
            disp(correlated_struct(i).green_info(j).kymoname)
        end 
    end
end

%% Colocalization plots test

for i = 16%:length(correlated_struct)
    if ~isempty(correlated_struct(i).green_info) && ~isempty(correlated_struct(i).red_info)
        for j = 1:length(correlated_struct(i).green_info)
            hold on
            axis off
            imagesc(flip(correlated_struct(i).red_info(j).kymograph))
            colormap(green_map)
            caxis([0 8])
            for k = 1:size(correlated_struct(i).correlated_data(j).green_data, 2)
                plot(correlated_struct(i).correlated_data(j).green_data{1, k}, correlated_struct(i).correlated_data(j).green_data{2, k}, "Color", 'g')
            end
            for k = 1:size(correlated_struct(i).correlated_data(j).red_data, 2)
%                 plot(correlated_struct(i).correlated_data(j).red_data{1, k}, correlated_struct(i).correlated_data(j).red_data{2, k}, "Color", 'r')
%                 plot(correlated_struct(i).correlated_data(j).red_data{3, k}, correlated_struct(i).correlated_data(j).red_data{4, k}, "Color", 'r')
%                 plot(correlated_struct(i).correlated_data(j).red_data{5, k}, correlated_struct(i).correlated_data(j).red_data{6, k}, "Color", 'r')
                plot(correlated_struct(i).correlated_data(j).red_data{7, k}, correlated_struct(i).correlated_data(j).red_data{8, k}, "Color", 'r')
            end
            hold off
        end
    end
end
%%
% for i = 1:length(correlated_struct)
%     if ~isempty(correlated_struct(i).green_info) && ~isempty(correlated_struct(i).red_info)
%         for j = 1:length(correlated_struct(i).green_info)
%             for k = 1:size(correlated_struct(i).correlated_data(j).red_data, 2)
%                 for l = 1:size(correlated_struct(i).correlated_data(j).green_data, 2)
%                     [these_intersect, i_red, i_green] =  intersect(correlated_struct(i).correlated_data(j).red_data{1, k},correlated_struct(i).correlated_data(j).green_data{1, l}); 
%                     correlated_struct(i).correlated_data(j).intersect_red_real{1,k} = these_intersect;
%                     correlated_struct(i).correlated_data(j).intersect_red_real{2,k} = i_red;
%                     correlated_struct(i).correlated_data(j).intersect_red_real{3,k} = i_green;
%                 end
%             end
%         end
%     end
% end
%%
% for i = 1:length(correlated_struct)
%     if ~isempty(correlated_struct(i).green_info) && ~isempty(correlated_struct(i).red_info)
%         for j = 1:length(correlated_struct(i).green_info)
%             for k = 1:size(correlated_struct(i).correlated_data(1).red_data, 2)
%                 for l = 1:size(correlated_struct(i).correlated_data(j).green_data, 2)
%                     these_intersect = [];
%                     i_red = [];
%                     i_green = [];
%                     [these_intersect, i_red, i_green] =  intersect(correlated_struct(i).correlated_data(j).red_data{1, k},correlated_struct(i).correlated_data(j).green_data{1, l});
%                     if ~isempty(these_intersect)
%                         correlated_struct(i).correlated_data(j).xintersect_red_real{k}{1,l} = these_intersect;
%                         correlated_struct(i).correlated_data(j).xintersect_red_real{k}{2,l} = i_red;
%                         correlated_struct(i).correlated_data(j).xintersect_red_real{k}{3,l} = i_green;
%                     end
%                     these_intersect = [];
%                     i_red = [];
%                     i_green = [];
%                     LIA =  ismembertol(correlated_struct(i).correlated_data(j).red_data{2, k},correlated_struct(i).correlated_data(j).green_data{2, l}, 0.01);
%                     if any(LIA)
%                         correlated_struct(i).correlated_data(j).yintersect_red_real{k}{1,l} = correlated_struct(i).correlated_data(j).red_data{2, k}(LIA);
%                     end
%                 end
%             end
%         end
%     end
% end

%% did it find the right things


% for i = 1:length(correlated_struct)
%     if ~isempty(correlated_struct(i).green_info) && ~isempty(correlated_struct(i).red_info)
%         for j = 1:length(correlated_struct(i).green_info)
%             fig = figure(i)
% %             axis off
% %             imagesc(flip(correlated_struct(i).green_info(j).kymograph))
% %             colormap(green_map)
% %             caxis([0 8])
%             counter = 1;
%             these_to_plot = [];
%             for z = 1:size(correlated_struct(i).correlated_data(j).yintersect_red_real, 2)
%                 if ~isempty(correlated_struct(i).correlated_data(j).yintersect_red_real{z})
%                     for zz = 1:length(correlated_struct(i).correlated_data(j).yintersect_red_real{z})
%                         if ~isempty(correlated_struct(i).correlated_data(j).yintersect_red_real{z}{zz})
%                             these_to_plot(counter)= zz;
%                             counter = counter +1;
%                         end
%                     end
%                 end
%             end
%             hold on
%             for k = these_to_plot
%                 plot(correlated_struct(i).correlated_data(j).green_data{1, k}, correlated_struct(i).correlated_data(j).green_data{2, k}, "Color", 'g')
%             end
%             for k = 1:size(correlated_struct(i).correlated_data(j).red_data, 2)
% %                 plot(correlated_struct(i).correlated_data(j).red_data{1, k}, correlated_struct(i).correlated_data(j).red_data{2, k}, "Color", 'r')
% %                 plot(correlated_struct(i).correlated_data(j).red_data{3, k}, correlated_struct(i).correlated_data(j).red_data{4, k}, "Color", 'r')
%                 plot(correlated_struct(i).correlated_data(j).red_data{5, k}, correlated_struct(i).correlated_data(j).red_data{6, k}, "Color", 'r')
% 
%             end
%             hold off
%         end
%     end
% end
%%


for i = 1:length(correlated_struct)
    if ~isempty(correlated_struct(i).green_info) && ~isempty(correlated_struct(i).red_info)
        for j = 1:length(correlated_struct(i).green_info)
            for k = 1:size(correlated_struct(i).correlated_data(1).red_data, 2)
                for l = 1:size(correlated_struct(i).correlated_data(j).green_data, 2)
                    correlated_struct(i).correlated_data(j).difference{k, l} = correlated_struct(i).correlated_data(j).green_data{2,l} - correlated_struct(i).correlated_data(j).red_data{8,k}(1);
                    correlated_struct(i).correlated_data(j).difference_abs{k, l} = abs(correlated_struct(i).correlated_data(j).green_data{2,l} - correlated_struct(i).correlated_data(j).red_data{8,k}(1));
                    correlated_struct(i).correlated_data(j).small_diff{k, l} = any(correlated_struct(i).correlated_data(j).difference_abs{k, l}<1);
                end
            end
        end
    end
end
%%
Yes = 0;
N0 = 0;
for i = 1:length(correlated_struct)
    if ~isempty(correlated_struct(i).green_info) && ~isempty(correlated_struct(i).red_info)
        for j = 1:length(correlated_struct(i).green_info)
            if j == 1
                color_is = "r";
            else
                color_is = "m";
            end
            for k = 1:size(correlated_struct(i).correlated_data(1).red_data, 2)
                for l = 1:size(correlated_struct(i).correlated_data(j).green_data, 2)
                    if correlated_struct(i).correlated_data(j).small_diff{k, l}
%                         fig = figure;
%                         hold on
                        Yes = Yes +1;
                        pointsare = [];
                        a = [];
                        b = [];
                        [pointsare, a, b] = intersect(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).red_data{5, k});
%                         plot(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).green_data{2, l}, "Color", 'g')
%                         plot(correlated_struct(i).correlated_data(j).red_data{5, k}(b), correlated_struct(i).correlated_data(j).red_data{6, k}(b), "Color", color_is)
%                         hold off
%                         hold on
%                         plot(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).difference{k, l})
%                         plot(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).difference_abs{k, l})
%                         hold off
                        avg_difference = mean(correlated_struct(i).correlated_data(j).difference{k, l});
                        if avg_difference > 0
                            other_side_amount = correlated_struct(i).correlated_data(j).difference_abs{k, l}(correlated_struct(i).correlated_data(j).difference{k, l}<0);
                        else
                            other_side_amount = correlated_struct(i).correlated_data(j).difference_abs{k, l}(correlated_struct(i).correlated_data(j).difference{k, l}>0);
                        end
                        if ~isempty(other_side_amount)
                          correlated_struct(i).correlated_data(j).other_side_amount{k, l} = other_side_amount;
                        else 
                          correlated_struct(i).correlated_data(j).other_side_amount{k, l} = [];
                        end
                    else
                        No = No + 1;
                    end
                end
            end
        end
    end
end

%% Plot just the portions over close proximity
Yes = 0;
N0 = 0;
for i = 1:length(correlated_struct)
    if ~isempty(correlated_struct(i).green_info) && ~isempty(correlated_struct(i).red_info)
        for j = 1:length(correlated_struct(i).green_info)
            if j == 1
                color_is = "r";
            else
                color_is = "m";
            end
            for k = 1:size(correlated_struct(i).correlated_data(1).red_data, 2)
                for l = 1:size(correlated_struct(i).correlated_data(j).green_data, 2)
                    if correlated_struct(i).correlated_data(j).small_diff{k, l}
                        fig = figure;
                        hold on
                        Yes = Yes +1;
                        pointsare = [];
                        a = [];
                        b = [];
                        [pointsare, a, b] = intersect(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).red_data{5, k});
                        plot(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).green_data{2, l}, "Color", 'g')
                        plot(correlated_struct(i).correlated_data(j).red_data{5, k}(b), correlated_struct(i).correlated_data(j).red_data{6, k}(b), "Color", color_is)
                        hold off
                    else
                        No = No + 1;
                    end
                end
            end
        end
    end
end

%%



Yes = 0;
N0 = 0;
for i = 1:length(correlated_struct)
    if ~isempty(correlated_struct(i).green_info) && ~isempty(correlated_struct(i).red_info)
        for j = 1:length(correlated_struct(i).green_info)
            if j == 1
                color_is = "r";
            else
                color_is = "m";
            end
            for k = 1:size(correlated_struct(i).correlated_data(1).red_data, 2)
                for l = 1:size(correlated_struct(i).correlated_data(j).green_data, 2)
                    if correlated_struct(i).correlated_data(j).small_diff{k, l} && ~isempty(correlated_struct(i).correlated_data(j).other_side_amount{k, l}) && any(correlated_struct(i).correlated_data(j).other_side_amount{k, l}>2)
%                         fig = figure;
                        hold on
                        Yes = Yes +1;
                        pointsare = [];
                        a = [];
                        b = [];
                        [pointsare, a, b] = intersect(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).red_data{5, k});
                        plot(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).green_data{2, l}, "Color", 'g')
                        plot(correlated_struct(i).correlated_data(j).red_data{5, k}(b), correlated_struct(i).correlated_data(j).red_data{6, k}(b), "Color", color_is)
%                         hold off
                    elseif correlated_struct(i).correlated_data(j).small_diff{k, l} && ~isempty(correlated_struct(i).correlated_data(j).other_side_amount{k, l}) && ~any(correlated_struct(i).correlated_data(j).other_side_amount{k, l}>2)
                        pointsare = [];
                        a = [];
                        b = [];
                        [pointsare, a, b] = intersect(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).red_data{5, k});
                        plot(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).green_data{2, l}, "Color", 'b')
                        plot(correlated_struct(i).correlated_data(j).red_data{5, k}(b), correlated_struct(i).correlated_data(j).red_data{6, k}(b), "Color", color_is)
                    elseif correlated_struct(i).correlated_data(j).small_diff{k, l} 
                        pointsare = [];
                        a = [];
                        b = [];
                        [pointsare, a, b] = intersect(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).red_data{5, k});
                        plot(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).green_data{2, l}, "Color", 'c')
                        plot(correlated_struct(i).correlated_data(j).red_data{5, k}(b), correlated_struct(i).correlated_data(j).red_data{6, k}(b), "Color", color_is)
                    end
                end
            end
        end
    end
end
hold off
%% Colocalizations and Crossovers
Yes = 1;
for i = 1:length(correlated_struct)
    disp(correlated_struct(i).DNA)
    if ~isempty(correlated_struct(i).green_info) && ~isempty(correlated_struct(i).red_info)
        for j = 1:length(correlated_struct(i).green_info)
            correlated_struct(i).correlated_data(j).kymoname =correlated_struct(i).green_info(j).kymoname; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            length_of_kymo = [];
            length_of_kymo = size(correlated_struct(i).green_info(j).kymograph,2);
            for k = 1:length(correlated_struct(i).green_info(j).particle_tracked)
                green_chords = [];
                green_position =[];
                green_chords = correlated_struct(i).green_info(j).crop_coordinates{k};
                green_position = correlated_struct(i).green_info(j).particle_tracked{k};
                correlated_struct(i).correlated_data(j).green_data{1, k}= green_position(:,1)+green_chords(1)-1;
                correlated_struct(i).correlated_data(j).green_data{2, k}= green_position(:,2)+green_chords(3)-1;
                correlated_struct(i).correlated_data(j).green_uncertainty{k} = sqrt(abs(correlated_struct(i).green_info(j).yinter{k}));
%                 NOTE TO SELF: The sqrt(variance) gives you the standard
%                 deviation of the data, which is the root mean squared
%                 differences between the data points. Which means this is
%                 the radius of the localization uncertainty around the
%                 particle. 
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for m = 1:length(correlated_struct(i).red_info(1).particle_tracked)
                red_chords = correlated_struct(i).red_info(1).crop_coordinates{m};
                red_position = correlated_struct(i).red_info(1).particle_tracked{m};
                correlated_struct(i).correlated_data(j).red_data{1, m}= red_position(:,1)+red_chords(1)-1;
                correlated_struct(i).correlated_data(j).red_data{2, m}= red_position(:,2)+red_chords(3)-1;
                correlated_struct(i).correlated_data(j).red_data{5, m}= (red_position(:,1)+red_chords(1)-1)';
                correlated_struct(i).correlated_data(j).red_data{6, m}= (red_position(:,2)+red_chords(3)-1)';
                length_of_trace = length(red_position);
                red_avg_position = mean(red_position(:,2));
                red_avg_pos_for_full_kymo = red_avg_position+red_chords(3)-1;
                x_is = [length_of_trace+red_chords(1)-1:length_of_kymo];
                x_is_2 = [red_chords(1):length_of_kymo];
                red_avg_pos_for_full_kymo_plot = red_avg_pos_for_full_kymo.* ones(length(x_is),1);
                red_avg_pos_for_full_kymo_plot_2 = red_avg_pos_for_full_kymo.* ones(length(x_is_2),1);
                correlated_struct(i).correlated_data(j).red_data{3, m}= x_is;
                correlated_struct(i).correlated_data(j).red_data{4, m}= red_avg_pos_for_full_kymo_plot;
                correlated_struct(i).correlated_data(j).red_data{5, m}= [correlated_struct(i).correlated_data(j).red_data{5, m} x_is];
                correlated_struct(i).correlated_data(j).red_data{6, m}= [correlated_struct(i).correlated_data(j).red_data{6, m} red_avg_pos_for_full_kymo_plot'];
                correlated_struct(i).correlated_data(j).red_data{7, m}= x_is_2;
                correlated_struct(i).correlated_data(j).red_data{8, m}= red_avg_pos_for_full_kymo_plot_2';
                correlated_struct(i).correlated_data(j).red_uncertainty{m} = sqrt(abs(correlated_struct(i).red_info(1).yinter{m}));
            end
            disp(correlated_struct(i).green_info(j).kymoname)
        end 
    end
end

for i = 1:length(correlated_struct)
    if ~isempty(correlated_struct(i).green_info) && ~isempty(correlated_struct(i).red_info)
        for j = 1:length(correlated_struct(i).green_info)
            for k = 1:size(correlated_struct(i).correlated_data(1).red_data, 2)
                for l = 1:size(correlated_struct(i).correlated_data(j).green_data, 2)
                    correlated_struct(i).correlated_data(j).difference{k, l} = correlated_struct(i).correlated_data(j).green_data{2,l} - correlated_struct(i).correlated_data(j).red_data{8,k}(1);
                    correlated_struct(i).correlated_data(j).difference_abs{k, l} = abs(correlated_struct(i).correlated_data(j).green_data{2,l} - correlated_struct(i).correlated_data(j).red_data{8,k}(1));
                    correlated_struct(i).correlated_data(j).small_diff{k, l} = any(correlated_struct(i).correlated_data(j).difference_abs{k, l}<0.7);%%%%%%%%%%%%%%%%%%%%%%CRITERIA FOR CROSSOVER%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
        end
    end
end


for i = 1:length(correlated_struct)
    if ~isempty(correlated_struct(i).green_info) && ~isempty(correlated_struct(i).red_info)
        for j = 1:length(correlated_struct(i).green_info)
            for k = 1:size(correlated_struct(i).correlated_data(1).red_data, 2)
                for l = 1:size(correlated_struct(i).correlated_data(j).green_data, 2)
                    if correlated_struct(i).correlated_data(j).small_diff{k, l}
                        Yes = Yes +1;
                        pointsare = [];
                        a = [];
                        b = [];
                        [pointsare, a, b] = intersect(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).red_data{5, k});
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
                    end
                end
            end
        end
    end
end

%% Plot colocalizations and crossovers

for i = 1:length(correlated_struct)
    if ~isempty(correlated_struct(i).green_info) && ~isempty(correlated_struct(i).red_info)
        for j = 1:length(correlated_struct(i).green_info)
            if j == 1
                color_is = 'r';
            else
                color_is = 'm';
            end
            for k = 1:size(correlated_struct(i).correlated_data(1).red_data, 2)
                for l = 1:size(correlated_struct(i).correlated_data(j).green_data, 2)
                    if correlated_struct(i).correlated_data(j).small_diff{k, l} && ~isempty(correlated_struct(i).correlated_data(j).other_side_amount{k, l}) && any(correlated_struct(i).correlated_data(j).other_side_amount{k, l}>2)
% % %                         fig = figure;
% %                         hold on
% %                         Yes = Yes +1;
% %                         pointsare = [];
% %                         a = [];
% %                         b = [];
% %                         b_len = [];
% %                         green_errorbar = [];
% %                         red_errorbar = [];
% %                         [pointsare, a, b] = intersect(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).red_data{5, k});
% %                         b_len = length(b);
% %                         len_green = length(correlated_struct(i).correlated_data(j).green_data{1, l});
% %                         green_errorbar = ones(1, len_green)*10*(correlated_struct(i).correlated_data(j).green_uncertainty{l}^2);
% %                         red_errorbar = ones(1, b_len)*10*(correlated_struct(i).correlated_data(j).red_uncertainty{k}^2);
% %                         shadedErrorBar(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).green_data{2, l},green_errorbar, 'lineProps', 'g')
% %                         shadedErrorBar(correlated_struct(i).correlated_data(j).red_data{5, k}(b), correlated_struct(i).correlated_data(j).red_data{6, k}(b),red_errorbar, 'lineProps', 'r')
% %                         disp(["i = ", num2str(i)])
% %                         disp(["k = ", num2str(k)])
% %                         disp(["l = ", num2str(l)])
% % %                                                 hold off
                    elseif correlated_struct(i).correlated_data(j).small_diff{k, l} && ~isempty(correlated_struct(i).correlated_data(j).other_side_amount{k, l}) && ~any(correlated_struct(i).correlated_data(j).other_side_amount{k, l}>2)
%                         pointsare = [];
%                         a = [];
%                         b = [];
%                         [pointsare, a, b] = intersect(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).red_data{5, k});
%                         b_len = length(b);
%                         len_green = length(correlated_struct(i).correlated_data(j).green_data{1, l});
%                         green_errorbar = ones(1, len_green)*10*(correlated_struct(i).correlated_data(j).green_uncertainty{l}^2);                        
%                         red_errorbar = ones(1, b_len)*10*(correlated_struct(i).correlated_data(j).red_uncertainty{k})^2;
%                         shadedErrorBar(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).green_data{2, l},green_errorbar, 'lineProps', 'b')
%                         shadedErrorBar(correlated_struct(i).correlated_data(j).red_data{5, k}(b), correlated_struct(i).correlated_data(j).red_data{6, k}(b),red_errorbar, 'lineProps', color_is)
                    elseif correlated_struct(i).correlated_data(j).small_diff{k, l} 
                        pointsare = [];
                        a = [];
                        b = [];
                        [pointsare, a, b] = intersect(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).red_data{5, k});
                        b_len = length(b);
                        len_green = length(correlated_struct(i).correlated_data(j).green_data{1, l});
                        green_errorbar = ones(1, len_green)*10*(correlated_struct(i).correlated_data(j).green_uncertainty{l}^2);                        
                        red_errorbar = ones(1, b_len)*10*(correlated_struct(i).correlated_data(j).red_uncertainty{k})^2;
                        shadedErrorBar(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).green_data{2, l},green_errorbar, 'lineProps', 'c')
                        shadedErrorBar(correlated_struct(i).correlated_data(j).red_data{5, k}(b), correlated_struct(i).correlated_data(j).red_data{6, k}(b),red_errorbar, 'lineProps', color_is)
                    end
                end
            end
        end
    end
end
hold off

%% Plot colocalizations and crossovers with 0 error

for i = 1:length(correlated_struct)
    if ~isempty(correlated_struct(i).green_info) && ~isempty(correlated_struct(i).red_info)
        for j = 1:length(correlated_struct(i).green_info)
            if j == 1
                color_is = 'r';
            else
                color_is = 'm';
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
                        b_len = [];
                        green_errorbar = [];
                        red_errorbar = [];
                        [pointsare, a, b] = intersect(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).red_data{5, k});
                        b_len = length(b);
                        len_green = length(correlated_struct(i).correlated_data(j).green_data{1, l});
                        green_errorbar = ones(1, len_green)*0*(correlated_struct(i).correlated_data(j).green_uncertainty{l}^2);
                        red_errorbar = ones(1, b_len)*0*(correlated_struct(i).correlated_data(j).red_uncertainty{k}^2);
                        shadedErrorBar(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).green_data{2, l},green_errorbar, 'lineProps', 'g')
                        shadedErrorBar(correlated_struct(i).correlated_data(j).red_data{5, k}(b), correlated_struct(i).correlated_data(j).red_data{6, k}(b),red_errorbar, 'lineProps', 'r')
                        disp(["i = ", num2str(i)])
                        disp(["k = ", num2str(k)])
                        disp(["l = ", num2str(l)])
                        %                         hold off
                    elseif correlated_struct(i).correlated_data(j).small_diff{k, l} && ~isempty(correlated_struct(i).correlated_data(j).other_side_amount{k, l}) && ~any(correlated_struct(i).correlated_data(j).other_side_amount{k, l}>2)
                        pointsare = [];
                        a = [];
                        b = [];
                        [pointsare, a, b] = intersect(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).red_data{5, k});
                        b_len = length(b);
                        len_green = length(correlated_struct(i).correlated_data(j).green_data{1, l});
                        green_errorbar = ones(1, len_green)*0*(correlated_struct(i).correlated_data(j).green_uncertainty{l}^2);                        
                        red_errorbar = ones(1, b_len)*0*(correlated_struct(i).correlated_data(j).red_uncertainty{k})^2;
                        shadedErrorBar(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).green_data{2, l},green_errorbar, 'lineProps', 'b')
                        shadedErrorBar(correlated_struct(i).correlated_data(j).red_data{5, k}(b), correlated_struct(i).correlated_data(j).red_data{6, k}(b),red_errorbar, 'lineProps', color_is)
                    elseif correlated_struct(i).correlated_data(j).small_diff{k, l} 
                        pointsare = [];
                        a = [];
                        b = [];
                        [pointsare, a, b] = intersect(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).red_data{5, k});
                        b_len = length(b);
                        len_green = length(correlated_struct(i).correlated_data(j).green_data{1, l});
                        green_errorbar = ones(1, len_green)*0*(correlated_struct(i).correlated_data(j).green_uncertainty{l}^2);                        
                        red_errorbar = ones(1, b_len)*0*(correlated_struct(i).correlated_data(j).red_uncertainty{k})^2;
                        shadedErrorBar(correlated_struct(i).correlated_data(j).green_data{1, l}, correlated_struct(i).correlated_data(j).green_data{2, l},green_errorbar, 'lineProps', 'c')
                        shadedErrorBar(correlated_struct(i).correlated_data(j).red_data{5, k}(b), correlated_struct(i).correlated_data(j).red_data{6, k}(b),red_errorbar, 'lineProps', color_is)
                    end
                end
            end
        end
    end
end
hold off

%% Average Localization Uncertainty Determination
% Green Localization Uncertainty 
counter = 1;
for i = 1:length(correlated_struct)
    if ~isempty(correlated_struct(i).green_info) && ~isempty(correlated_struct(i).red_info)
        for j = 1:length(correlated_struct(i).green_info)
            for k = 1:length(correlated_struct(i).green_info(j).yinter)
                green_intercepts(counter) = correlated_struct(i).green_info(j).yinter{k};
                counter = counter + 1;
            end
        end
    end
end
mean(green_intercepts)
mean(sqrt(abs(green_intercepts))) % Localization Precision is 0.2um which is +/- 2pixels, pretty bad


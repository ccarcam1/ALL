%% Add Cas9 data from the dual color imaging 
b = 0;
c = 1;
a = 1;
for i = 1:length(combined)
    if combined(i).red_exist == 1 && combined(i).green_exist == 1;
        for k = 1%:length(combined(i).red_info)
            for j = 1:length(combined(i).red_info(k).particle_tracked)
%             Cas9(a+b).number = c;
            Cas9(a+b).name = strcat(combined(i).DNA," ", combined(i).red_info(k).kymoname, " ", combined(i).date);
            Cas9(a+b).DNA = strcat(combined(i).DNA, " ", combined(i).date);
            Cas9(a+b).experiment = "cas9 in dual color";
            Cas9(a+b).line_time = combined(i).red_info(k).line_time{j};
            Cas9(a+b).particle_tracked = combined(i).red_info(k).particle_tracked{j};   
            Cas9(a+b).MSD = combined(i).red_info(k).MSD{j};
            Cas9(a+b).crop_coords = combined(i).red_info(k).crop_coordinates{j};
            Cas9(a+b).crop = combined(i).red_info(k).crop{j};
            Cas9(a+b).kymograph = combined(i).red_info(k).kymograph;
            Cas9(a+b).pvalue = combined(i).red_info(k).pval{j};
            Cas9(a+b).rsquared = combined(i).red_info(k).r2{j};
            Cas9(a+b).yintercept = combined(i).red_info(k).yinter{j};
            Cas9(a+b).fitcutoff = NaN;
            Cas9(a+b).fitcutoff_NaN = NaN;
            Cas9(a+b).slope = combined(i).red_info(k).diffusion_coeff{j}*2;
            Cas9(a+b).length_k = size(combined(i).red_info(k).kymograph,2);
            a = a+1;
            end
%         unique_names_Cas9{c} = strcat(combined(i).DNA," ", combined(i).red_info(k).kymoname," ", combined(i).date);             
        unique_names_Cas9{c} = strcat(combined(i).DNA," ", combined(i).date);
        c = c+1;
            
        end
    end
end
Cas9_name = {Cas9.name}.';
Cas9_DNA = {Cas9.DNA}.';
Cas9_particle_tracked = {Cas9.particle_tracked}.';
Cas9_crop_coords = {Cas9.crop_coords}.';
Cas9_yintercept = {Cas9.yintercept}.';
Cas9_length_k = {Cas9.length_k}.';
%% Add the SWR1 data from the two color data
b=0;
a = 1;
for i = 1:length(combined)
    if combined(i).green_exist == 1 && combined(i).red_exist == 1;
        for k = 1:length(combined(i).green_info)
            for j = 1:length(combined(i).green_info(k).particle_tracked)
            SWR1(a+b).name = strcat(combined(i).DNA," ", combined(i).green_info(k).kymoname, " ", combined(i).date);
            SWR1(a+b).DNA = strcat(combined(i).DNA," ", combined(i).date);
            SWR1(a+b).kymo = k;%strcat("kymo num ", num2str(k)," ", combined(i).date);
            SWR1(a+b).experiment = "SWR1 70mM KCl 1mM ATP dual color";
            SWR1(a+b).line_time = combined(i).green_info(k).line_time{j};
            SWR1(a+b).particle_tracked = combined(i).green_info(k).particle_tracked{j};   
            SWR1(a+b).MSD = combined(i).green_info(k).MSD{j};
            SWR1(a+b).crop_coords = combined(i).green_info(k).crop_coordinates{j};
            SWR1(a+b).crop = combined(i).green_info(k).crop{j};
            SWR1(a+b).kymograph = combined(i).green_info(k).kymograph;
            SWR1(a+b).pvalue = combined(i).green_info(k).pval{j};
            SWR1(a+b).rsquared = combined(i).green_info(k).r2{j};
            SWR1(a+b).yintercept = combined(i).green_info(k).yinter{j};
            SWR1(a+b).fitcutoff = NaN;
            SWR1(a+b).fitcutoff_NaN = NaN;
            SWR1(a+b).slope = combined(i).green_info(k).diffusion_coeff{j}*2;
            SWR1(a+b).length_k = size(combined(i).green_info(k).kymograph,2);
            a = a+1;
            end
        end
    end
end
SWR1_name = {SWR1.name}.';
SWR1_DNA = {SWR1.DNA}.';
SWR1_kymo = [SWR1.kymo].';
% SWR1_kymo = {SWR1.kymo}.';
SWR1_particle_tracked = {SWR1.particle_tracked}.';
SWR1_crop_coords = {SWR1.crop_coords}.';
SWR1_yintercept= {SWR1.yintercept}.';
SWR1_length_k = {SWR1.length_k}.';
% %% Find common locations:
% for j = 1:length(unique_names_Cas9)
%     for i = 1:length(SWR1_name)
%         if matches(SWR1_name{i}, unique_names_Cas9{j})
%             SWR1_name1(i) = 1; 
%         else
%             SWR1_name1(i) = 0;
%         end
%     end
%     SWR1_name_bools{j} = logical(SWR1_name1);
% end
% for j = 1:length(unique_names_Cas9)
%     for i = 1:length(Cas9_name)
%         if matches(Cas9_name{i}, unique_names_Cas9{j})
%             Cas9_name1(i) = 1; 
%         else
%             Cas9_name1(i) = 0;
%         end
%     end
%     Cas9_name_bools{j} = logical(Cas9_name1);
% end
%% Find common locations:
for j = 1:length(unique_names_Cas9)
    for i = 1:length(SWR1_DNA)
        if matches(SWR1_DNA{i}, unique_names_Cas9{j})
            SWR1_DNA1(i) = 1; 
        else
            SWR1_DNA1(i) = 0;
        end
    end
    SWR1_name_bools{j} = logical(SWR1_DNA1);
end
for j = 1:length(unique_names_Cas9)
    for i = 1:length(Cas9_DNA)
        if matches(Cas9_DNA{i}, unique_names_Cas9{j})
            Cas9_DNA1(i) = 1; 
        else
            Cas9_DNA1(i) = 0;
        end
    end
    Cas9_name_bools{j} = logical(Cas9_DNA1);
end
%% Organize the trajectory data
for i = 1:length(unique_names_Cas9)
    if sum(SWR1_name_bools{i})>0
        SWR1_where = SWR1_name_bools{i};
        Cas9_where = Cas9_name_bools{i};
        SWR1_particle_info = {};
        Cas9_particle_info = {};
        SWR1_chords = {};
        Cas9_chords = {};
        SWR1_len_k = {};
        SWR1_kymonum = [];
        SWR1_kymograph_cell = {};
        Cas9_kymograph_cell = {};
        SWR1_particle_info = SWR1_particle_tracked(SWR1_where);
        Cas9_particle_info = Cas9_particle_tracked(Cas9_where);
        SWR1_chords = SWR1_crop_coords(SWR1_where);
        Cas9_chords = Cas9_crop_coords(Cas9_where);
        SWR1_len_k = SWR1_length_k(SWR1_where);
        SWR1_yinter = SWR1_yintercept(SWR1_where);
        Cas9_yinter = Cas9_yintercept(Cas9_where);
        SWR1_kymonum = SWR1_kymo(SWR1_where);
        SWR1_kymograph_cell = SWR1_kymograph(SWR1_where);
        Cas9_kymograph_cell = Cas9_kymograph(Cas9_where);
        holdingplace = {};
        for m = 1:length(Cas9_particle_info)
            red_chords = Cas9_chords{m};
            red_position = Cas9_particle_info{m};
            holdingplace{1, m}= red_position(:,1)+red_chords(1)-1;
            holdingplace{2, m}= red_position(:,2)+red_chords(3)-1;
            holdingplace{5, m}= (red_position(:,1)+red_chords(1)-1)';
            holdingplace{6, m}= (red_position(:,2)+red_chords(3)-1)';
            length_of_trace = length(red_position);
            red_avg_position = mean(red_position(:,2));
            red_avg_pos_for_full_kymo = red_avg_position+red_chords(3)-1;
            x_is = [length_of_trace+red_chords(1)-1:SWR1_len_k{1}];
            x_is_2 = [red_chords(1):SWR1_len_k{1}];
            red_avg_pos_for_full_kymo_plot = red_avg_pos_for_full_kymo.* ones(length(x_is),1);
            red_avg_pos_for_full_kymo_plot_2 = red_avg_pos_for_full_kymo.* ones(length(x_is_2),1);
            holdingplace{3, m}= x_is;
            holdingplace{4, m}= red_avg_pos_for_full_kymo_plot;
            holdingplace{5, m}= [holdingplace{5, m} x_is];
            holdingplace{6, m}= [holdingplace{6, m} red_avg_pos_for_full_kymo_plot'];
            holdingplace{7, m}= x_is_2;
            holdingplace{8, m}= red_avg_pos_for_full_kymo_plot_2';
            holdingplace{9, m}= sqrt(abs(Cas9_yinter{m}));
            holdingplace{10, m}= std(red_position(:,2));
            holdingplace{11, m} = std(red_position(:,2))/sqrt(length(red_position(:,2)));
            holdingplace{12, m} = Cas9_kymograph_cell{1};
        end
        Cas9_info{i} = holdingplace;
        holdingplace = {};
        for m = 1:length(SWR1_particle_info)
            green_chords = SWR1_chords{m};
            green_position = SWR1_particle_info{m};
            holdingplace{1, m}= green_position(:,1)+green_chords(1)-1;
            holdingplace{2, m}= green_position(:,2)+green_chords(3)-1;
            holdingplace{3, m}= sqrt(abs(SWR1_yinter{m}));
            holdingplace{4, m} = SWR1_kymonum(m);
            holdingplace{5, m} = SWR1_kymograph_cell{m};
        end
        SWR1_info{i} = holdingplace;
    else
        disp(i)
    end
end

%% Find the distance to nearest dCas9
for i = 1:length(unique_names_Cas9)
%     fig= figure(j);
%     set(fig, "Visible", "on");
%     hold on
    holdingplace2 = {};
    green_info = {};
    red_info = {};
    absoluteval= {};
    small_val= {};
    small_val_yes_no= {};
    time_overlap = {};
    SWR1_where = SWR1_name_bools{i};
    Cas9_where = Cas9_name_bools{i};
    for m = 1:sum(SWR1_where)
        for n = 1:sum(Cas9_where)
            C= [];
            ia = [];
            ib = [];
            [C,ia,ib] = intersect(SWR1_info{i}{1, m},Cas9_info{i}{5, n});
            time_overlap{m, n} = ia;
            holdingplace2{m, n} = SWR1_info{i}{2, m}(ia) - [Cas9_info{i}{6, n}(ib)]';
            green_info{m, n} = [SWR1_info{i}{1, m}(ia)'; SWR1_info{i}{2, m}(ia)'];
            red_info{m, n} = [Cas9_info{i}{5, n}(ib); Cas9_info{i}{6, n}(ib)];
            absoluteval{m, n} = abs(holdingplace2{m, n});
            small_val{m, n} = absoluteval{m, n}<1;
            small_val_yes_no{m, n} = any(small_val{m, n});
% %             [C,ia,ib] = intersect(SWR1_info{i}{1, m},Cas9_info{i}{7, n});
% %             time_overlap{m, n} = ia;
% %             holdingplace2{m, n} = SWR1_info{i}{2, m}(ia) - [Cas9_info{i}{8, n}(ia)]';
% %             green_info{m, n} = [SWR1_info{i}{1, m}(ia)'; SWR1_info{i}{2, m}(ia)'];
% %             red_info{m, n} = [Cas9_info{i}{7, n}(ib); Cas9_info{i}{8, n}(ib)];
% %             absoluteval{m, n} = abs(holdingplace2{m, n});
% %             small_val{m, n} = absoluteval{m, n}<1;
% %             small_val_yes_no{m, n} = any(small_val{m, n});
        end
    end
    dual_info{i, 1} = holdingplace2;
    dual_info{i, 2} = green_info;
    dual_info{i, 3} = red_info;
    dual_info{i, 4} = absoluteval;
    dual_info{i, 5} = small_val;
    dual_info{i, 6} = small_val_yes_no;
    time_overlap_all{i} = time_overlap;
%     hold off
end

%%

for i = 1:length(unique_names_Cas9)
    fig= figure(j);
    set(fig, "Visible", "off");
    hold on
    SWR1_where = SWR1_name_bools{i};
    Cas9_where = Cas9_name_bools{i};
    for m = 1:sum(SWR1_where)
        for n = 1:sum(Cas9_where)
            if SWR1_info{i}{4, m} == 1
                color_SWR1 = 'g';
            else
                color_SWR1 = 'b';
            end
            if dual_info{i, 6}{m, n} == 1
                plot(dual_info{i, 2}{m, n}(1,:),dual_info{i, 2}{m, n}(2,:), "Color",color_SWR1) 
                plot(dual_info{i, 3}{m, n}(1,:),dual_info{i, 3}{m, n}(2,:), "Color",'r')
            disp("Yes")
            end
        end
    end
    hold off
    ylim([0, 180])
    if any([dual_info{i,6}{:}] == 1)
        saveas(fig, unique_names_Cas9{i}, "png")
    end
    close all
end

%%
for i = 1%:length(unique_names_Cas9)
    fig= figure(i);
    set(fig, "Visible", "on");
    hold on
    SWR1_where = SWR1_name_bools{i};
    Cas9_where = Cas9_name_bools{i};
    for m = 1:sum(SWR1_where)
        if SWR1_info{i}{4, m} == 1
            color_SWR1 = 'g';
        else
            color_SWR1 = 'b';
        end
            plot(SWR1_info{i}{1, m},SWR1_info{i}{2, m},"Color",color_SWR1)  
    end
    for n = 1:sum(Cas9_where)
%         plot(Cas9_info{i}{7, n},Cas9_info{i}{8, n},"Color",'r')
%         plot(Cas9_info{i}{1, n},Cas9_info{i}{2, n},"Color",'m')
%         plot(Cas9_info{i}{3, n},Cas9_info{i}{4, n},"Color",'y')
        plot(Cas9_info{i}{5, n},Cas9_info{i}{6, n},"Color",'k')
    end      
    hold off
%     saveas(fig, unique_names_Cas9{i}, "png")
%     close all
end



% % % %%
% % % for i = 1:length(unique_names_Cas9)
% % %     fig= figure(i);
% % %     set(fig, "Visible", "off");
% % %     SWR1_where = SWR1_name_bools{i};
% % %     Cas9_where = Cas9_name_bools{i};
% % % %     for m = 1:sum(SWR1_where)
% % % %         if SWR1_info{i}{4, m} == 1
% % % %             color_SWR1 = 'g';
% % % %         else
% % % %             color_SWR1 = 'b';
% % % %         end
% % %     RGB = rand(size(SWR1_info{i}{5, 1},1), size(SWR1_info{i}{5, 1}, 2), 3);
% % %     RGB(:,:,2) = SWR1_info{i}{5, 1};
% % %     RGB(:,:,1) = Cas9_info{i}{12, 1};
% % %     RGB(:,:,3) = 0;
% % %     imagesc(RGB)
% % % %     if any([dual_info{i,6}{:}] == 1)
% % %         saveas(fig, unique_names_Cas9{i}, "png")
% % % %     end
% % % %     figure(2);
% % % %     imagesc(SWR1_info{i}{5, 1})
% % % %     figure(3);
% % % %     imagesc(Cas9_info{i}{12, 1})
% % % %     end
% % % end
%% Tally by eye based on tracked colocalization data 

% % % 
% % % for i = 1:length(unique_names_Cas9)
% % %     fig= figure(j);
% % %     set(fig, "Visible", "on");
% % %     hold on
% % %     SWR1_where = SWR1_name_bools{i};
% % %     Cas9_where = Cas9_name_bools{i};
% % %     for m = 1:sum(SWR1_where)
% % %         for n = 1:sum(Cas9_where)
% % %             if SWR1_info{i}{4, m} == 1
% % %                 color_SWR1 = 'g';
% % %             else
% % %                 color_SWR1 = 'b';
% % %             end
% % %             if dual_info{i, 6}{m, n} == 1
% % %                 plot(dual_info{i, 2}{m, n}(1,:),dual_info{i, 2}{m, n}(2,:), "Color",color_SWR1) 
% % %                 plot(dual_info{i, 3}{m, n}(1,:),dual_info{i, 3}{m, n}(2,:), "Color",'r')
% % %             disp("Yes")
% % %             end
% % %         end
% % %     end
% % %     hold off
% % %     ylim([0, 180])
% % % %     if any([dual_info{i,6}{:}] == 1)
% % % %         saveas(fig, unique_names_Cas9{i}, "png")
% % % %     end
% % % counteris = 1;
% % % choice = menu('what kind of event is it?', 'Stuck','Bounce','Stick + Bounce','unlikely meaningful colocalization','Break');
% % % if choice == 1
% % %     colocs(counteris).type = 'stuck';
% % %     colocs(counteris).num_of_mols = input('number of molecules?');
% % %     choice2 = menu('go to the next one?', 'YES', 'NO');
% % %     if choice2 == 2
% % %         
% % % 
% % %     close all
% % % end

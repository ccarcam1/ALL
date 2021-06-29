function structure_name = my_second_pass_two_color_particle_picking(Green_MY_DATA, Red_MY_DATA, start_here, end_here)

green_map = zeros(3);
for i = 1:11
green_map(i,2) = i*0.1-0.1;
end

counter = 1;
for lmnop = start_here:end_here%length(uniq_kymo_names)
    if isempty(Green_MY_DATA(lmnop).Red_Info_Loc)
        continue
    end
    ii = find([Green_MY_DATA.Kymo_Number] == lmnop);
    i = ii(1);
%     threshold_for_red_trace_in_subplot = 80;
    num_of_red_traces = length(Green_MY_DATA(i).Red_Info_Loc);
    fh1 = figure;%('visible','off');% select to suppress figures when saving
    kymograph = rot90(Green_MY_DATA(i).kymograph, 3);
    imagesc(my_threshold_function(kymograph,30));
%     imagesc(flip(Green_MY_DATA(i).kymograph))
%     colormap(fh1,green_map)
    caxis([0 8])
    hold on

%Plots dCas9
for j = 1:num_of_red_traces
    red_chords = Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(j)).crop_coordinates;
    red_position = Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(j)).particle_tracked;
    red_avg_position = mean(red_position(:,2));
    red_avg_pos_for_full_kymo = red_avg_position+red_chords(3)-1;
    length_green_trace = length(Green_MY_DATA(i).particle_tracked);
    xline(red_avg_pos_for_full_kymo,'Color','r','LineWidth',2)
end
    hold off
    
%Plots SWR1   
    hold on
    for j = 1:length(ii)
        green_chords = Green_MY_DATA(ii(j)).crop_coordinates;
        green_position = Green_MY_DATA(ii(j)).particle_tracked;
        plot(green_position(:,2)+green_chords(3)-1,green_position(:,1)+green_chords(1)-1,'Color','g','LineWidth',2)
    end
    hold off


        counterthis = counter;
        j = 1;
        ksize = size(kymograph);
        kw=ksize(:,2);
        kh=ksize(:,1);
        filename = [Green_MY_DATA(i).name, '_tracenum_', num2str(j)];
        title(filename);
        prompt = 'Should skip this kymo. Type d to skip. else press enter key with blank:  ';
        skipthisone = input(prompt);
        close all
    while isempty(skipthisone)
%        num_of_red_traces = length(Green_MY_DATA(i).Red_Info_Loc);
    






fh1 = figure;%('visible','off');% select to suppress figures when saving
    kymograph = rot90(Green_MY_DATA(i).kymograph, 3);
    imagesc(my_threshold_function(kymograph,30));
%     imagesc(flip(Green_MY_DATA(i).kymograph))
%     colormap(fh1,green_map)
    caxis([0 8])
    hold on

%Plots dCas9
for j = 1:num_of_red_traces
    red_chords = Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(j)).crop_coordinates;
    red_position = Red_MY_DATA(Green_MY_DATA(i).Red_Info_Loc(j)).particle_tracked;
    red_avg_position = mean(red_position(:,2));
    red_avg_pos_for_full_kymo = red_avg_position+red_chords(3)-1;
    length_green_trace = length(Green_MY_DATA(i).particle_tracked);
    xline(red_avg_pos_for_full_kymo,'Color','r','LineWidth',2)
end
    hold off
    
%Plots SWR1   
    hold on
    for j = 1:length(ii)
        green_chords = Green_MY_DATA(ii(j)).crop_coordinates;
        green_position = Green_MY_DATA(ii(j)).particle_tracked;
        plot(green_position(:,2)+green_chords(3)-1,green_position(:,1)+green_chords(1)-1,'Color','g','LineWidth',2)
    end
    hold off
    
    
    
    
    
        ksize = size(kymograph);
        kw=ksize(:,2);
        kh=ksize(:,1);
        structure_name(counter).name = filename;
        if j ~= 1
            for asdf = counterthis:counter-1
                rectangle('Position',[structure_name(asdf).crop_coordinates(3),structure_name(asdf).crop_coordinates(1),...
                structure_name(asdf).crop_coordinates(4)-structure_name(asdf).crop_coordinates(3),structure_name(asdf).crop_coordinates(2)-...
                structure_name(asdf).crop_coordinates(1)], 'EdgeColor','y'); 
            end
        end
        zoom on;
        pause() % you can zoom with your mouse and when your image is okay, you press any key
        zoom off; % to escape the zoom mode
        p = ginput(2); 
        sp(1) = min(floor(p(1)), floor(p(2))); %xmin
        sp(2) = min(floor(p(3)), floor(p(4))); %ymin
        sp(3) = max(ceil(p(1)), ceil(p(2)));   %xmax
        sp(4) = max(ceil(p(3)), ceil(p(4)));   %ymax
        if sp(1)<1
            sp(1)=1;
        end
        if sp(2)<1
            sp(2)=1;
        end
        if sp(3)>(kw-1)
            sp(3)=kw;
        end
        if sp(4)>(kh-1)
            sp(4)=kh;
        end
        crop = kymograph(sp(2):sp(4), sp(1): sp(3),:);
        structure_name(counter).crop = crop;
        structure_name(counter).crop_coordinates = [sp(2),sp(4),sp(1),sp(3)];
        crop_fig = figure;
        if max(kymograph,[], 'all') > 30
            imagesc(crop,[0,max(crop, [], 'all')/2])  
          else
            imagesc(crop);      
          end
        ksize = size(crop);
        kw=ksize(:,2);
        kh=ksize(:,1);
        p = ginput(2); 
        sp(1) = min(floor(p(1)), floor(p(2))); %xmin
        sp(2) = min(floor(p(3)), floor(p(4))); %ymin
        sp(3) = max(ceil(p(1)), ceil(p(2)));   %xmax
        sp(4) = max(ceil(p(3)), ceil(p(4)));   %ymax
            if sp(1)<1
                sp(1)=1;
            end
            if sp(2)<1
                sp(2)=1;
            end
            if sp(3)>(kw-1)
                sp(3)=kw;
            end
            if sp(4)>(kh-1)
                sp(4)=kh;
            end
        BG = crop(sp(2):sp(4), sp(1): sp(3),:);
        BG = reshape(BG.',1,[]);
        d=mean(BG);
        dsd=sqrt(2)*std(BG);%% fix a mistake in the or
        dmin=d-(2*dsd);
        dmax=d+(2*dsd);
        structure_name(counter).dmin = dmin;
        structure_name(counter).dmax = dmax;        
        structure_name(counter).d = d;
        crop_fig = figure;
        imagesc(crop,[0,max(crop, [], 'all')/2]);
%         set(crop_fig,'Position',[594.0000, 33.3333, 334.0000, 612.6667]);
        p = ginput(2); 
        sp(11) = min(floor(p(1)), floor(p(2))); %xmin
        sp(12) = min(floor(p(3)), floor(p(4))); %ymin
        sp(13) = max(ceil(p(1)), ceil(p(2)));   %xmax
        sp(14) = max(ceil(p(3)), ceil(p(4)));   %ymax
            if sp(11)<1
                sp(11)=1;
            end
            if sp(12)<1
                sp(12)=1;
            end
            if sp(13)>(kw-1)
                sp(13)=kw;
            end
            if sp(14)>(kh-1)
                sp(14)=kh;
            end
        peak = crop(sp(12):sp(14), sp(11): sp(13),:);
        peak = reshape(peak.',1,[]);
        a=mean(peak);
        amin=a-(3*dsd);
        amax=a+(3*dsd);
        structure_name(counter).amin = amin;
        structure_name(counter).amax = amax;
        structure_name(counter).a = a;
        sz=size(crop);
        t=sz(:,1);
        w=sz(:,2);
        b=(sp(11)+sp(13))/2;
        bmin=b-5;
        bmax=b+5;
        structure_name(counter).bmin = bmin;
        structure_name(counter).bmax = bmax;
        structure_name(counter).b = b;
        wrange=[1:1:w];
        transpose(wrange);
        structure_name(counter).wrange = wrange;
        clf
        counter = counter + 1;
%         cd(segmentation_dir);
%         if colornumber == 1
%             save(['red_segmentation ',num2str(start_var), ' to ', num2str(end_var), '.mat'], 'structure_name');
%         elseif colornumber == 2
%             save(['green_segmentation ',num2str(start_var), ' to ', num2str(end_var), '.mat'], 'structure_name');
%         end
        close all
%         kymofig = figure;
%         if j ~= num_of_traces
%             imagesc(my_threshold_function(kymograph,30));
%             set(kymofig,'Position',[2.0000, 26.6667,589.6667,619.6667]); %laptop
%             ksize = size(kymograph);
%             kw=ksize(:,2);
%             kh=ksize(:,1);
%         else
%         end
    prompt = 'Should we break? To break enter number 1.';
    breakinner = input(prompt);
    if isempty(breakinner)
    else
        break
    end
%         zoom on;
%         pause() % you can zoom with your mouse and when your image is okay, you press any key
%         zoom off; % to escape the zoom mode
    j = j+1;
    
    end
end
% disp(['you are on trace ', num2str(i), ' out of ', num2str(length(data3))])
% end
    close
save(['two-color second pass segmentation ',num2str(start_here), ' to ', num2str(end_here), '.mat'], 'structure_name')
end


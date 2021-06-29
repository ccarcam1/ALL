function fitting_MSD_struct = my_condense_relevant_info2(color, data, gaussfitting, structure_name)
for i = 1:length(gaussfitting)
    fitting_MSD_struct(i).name = gaussfitting(i).name;
    fitting_MSD_struct(i).pixel = gaussfitting(i).data.pixel;
    fitting_MSD_struct(i).timepix = gaussfitting(i).data.timepix;
    fitting_MSD_struct(i).intensity = gaussfitting(i).data.intensity;
    fitting_MSD_struct(i).original_position = gaussfitting(i).original_position;
    for j = 1:length(structure_name)
        if strcmp(fitting_MSD_struct(i).name, structure_name(j).name)
            fitting_MSD_struct(i).crop = structure_name(j).crop;
            fitting_MSD_struct(i).crop_coordinates = structure_name(j).crop_coordinates;
        else
        end
    end
    for j = 1:length(data)
        if contains(fitting_MSD_struct(i).name, strcat(data(j).name, '_'))
            if strcmp(color, 'green')
                fitting_MSD_struct(i).full_kymo = data(j).green_kymo;
                fitting_MSD_struct(i).line_time = data(j).line_time;
            elseif strcmp(color, 'red')
                fitting_MSD_struct(i).full_kymo = data(j).red_kymo;
                fitting_MSD_struct(i).line_time = data(j).line_time;
            end
        else
        end
    end    
end
end

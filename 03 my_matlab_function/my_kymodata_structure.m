function my_kymodata_structure(color)
% Initialize
% color= 'green';
coloris = color;
kymoname = {}; 
kymofilename = {};
file_names = {};
%Load matlab structures containing path information
mat = dir('*.mat'); 
for q = 1:length(mat)
    load(mat(q).name);
end
if contains(coloris, 'red')
    redkymo = dir(kymo_mat_red);
elseif contains (coloris, 'green')
    greenkymo = dir(kymo_mat_green);
end
linetime = dir(linetime_matstruct);
pattern = [".", ".."];
counter = 1;
for i = 1:length(linetime) % Get name of kymos
    if not(startsWith(linetime(i).name, pattern))
    str = linetime(i).name;
    match = [".mat"];
    name = erase(str,match);
%     name = strip(name,'left','_');
%     name = name(1:(end-14));
%     name = name(1:(end-11));
    kymoname{counter,1} = name;
    counter = counter +1;
    end
end
disp('.....');
listing = dir();
stringis = strcat(listing.name);
if contains(stringis, 'green')
  cd(kymo_mat_green)
  counter = 1;
    for i = 1:length(greenkymo) %Get file names of green kymos
        if not(startsWith(greenkymo(i).name, pattern))
        name = greenkymo(i).name;
        kymofilename_green{counter,1} = name;
        counter = counter +1;
        end
    end
    for i = 1: length(kymoname) % Get kymo object array
        file_to_open = find(contains(kymofilename_green,kymoname{i}));
        load (kymofilename_green{file_to_open})
        kymoname{i,2} = double(obj_arr);
    end  
else
end

disp('....');
if contains(stringis, 'red')
    cd(kymo_mat_red)
    counter = 1;
    for i = 1:length(redkymo) %Get file names of red kymos
        if not(startsWith(redkymo(i).name, pattern))
        name = redkymo(i).name;
        kymofilename_red{counter,1} = name;
        counter = counter +1;
        end
    end
    for i = 1: length(kymoname) % Get kymo object array
        file_to_open = find(contains(kymofilename_red,kymoname{i}));
        load (kymofilename_red{file_to_open})
        kymoname{i,3} = double(obj_arr);
    end
end
disp('...');
cd(linetime_matstruct)
counter = 1;
for i = 1:length(linetime) %Get file names of linetimes
    if not(startsWith(linetime(i).name, pattern))
    name = linetime(i).name;
    linetime_filename{counter,1} = name;
    counter = counter +1;
    end
end
for i = 1: length(kymoname) % Get linetimes
    file_to_open = find(contains(linetime_filename,kymoname{i}));
    load (linetime_filename{file_to_open})
    kymoname{i,4} = double(linescan_time_);
end
disp('..');
for i = 1:length(kymoname) % Save into a structured array
    data(i).name = kymoname{i,1};
    data(i).green_kymo = kymoname{i,2};
    data(i).red_kymo = kymoname{i,3};
    data(i).line_time = kymoname{i,4};
end
cd(container_path)
save('data.mat', 'data');
disp('Done.');
end
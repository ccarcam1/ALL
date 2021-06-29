%% The purpose of this code is to quickly rename the results of linking 
%  data created using kymograph direct. This way the data can be compiled
%  and analyzed using msdanalyzer matlab package. Claudia Carcamo 8-12-19
% cd '2019-08-27 SWR1 ATP'; % Two directories out 
% A = dir('New folder'); % One directory out 
A = dir; % One directory out 

name_of_trace = {};

for i = 3: size(A,1)
    name_of_trace{i-2} = A(i).name;
end

%% Rename the results of the linked traces with the DNA and tracenumber
% cd '2019-08-27 SWR1 ATP';
A = dir('all_traces');
name_of_trace = {};

for i = 3: size(A,1)
    name_of_trace{i-2} = A(i).name;
end

cd 'all_traces';
for i = 1: size(name_of_trace,2)
    cd (name_of_trace{i});
    cd 'Linked lines results';
    movefile('Intensity_vs_position.txt', ['Intensity_vs_position ' name_of_trace{i} '.txt']);
    movefile('lines_coordinates_in_pixel.txt', ['lines_coordinates_in_pixel ' name_of_trace{i} '.txt']);
    movefile('Particles_positions_vs_time.txt', ['Particles_positions_vs_time ' name_of_trace{i} '.txt']);
    movefile('Velocity_vs_position.txt', ['Velocity_vs_position ' name_of_trace{i} '.txt']);
    cd  ..\..\
%     disp(name_of_trace{i})
end

%% Move all results to the correct folder
for i = 1: size(name_of_trace,2)
    cd (name_of_trace{i});
    cd 'Linked lines results';
    copyfile ** ..\..\'All Linked Traces';
    cd  ..\..\;
end

%% Rename files in folder (If the above three sections were not performed
%  previously, this section can rename files quickly
%  Part 1

A = dir('All Linked Traces');
name_of_trace = {};

for i = 3: size(A,1)
    name_of_trace{i-2} = A(i).name;
end

%% Next part for renaming files in folder 
%  Part 2
for i = 1: size(name_of_trace,2)
     movefile(['' name_of_trace{i} ''], ['3-' name_of_trace{i} '']);
end


%% Next part for renaming files in folder 
%  Part 2
A = dir('All Linked Traces 2019-08-20');
for i = 3: size(A,1)
    name_of_trace{i-2} = A(i).name;
end
%%

for i = 1: size(name_of_trace,2)
     movefile(['' name_of_trace{i} ''], ['cas9 day1 ' name_of_trace{i} '']);
end
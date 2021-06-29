function my_directory_function(varargin)

path_start = strcat(pwd, '\');
listing = dir();
stringis = strcat(listing.name);
if not(contains(stringis, 'container'))
    mkdir container;
else
    rmdir container s;
    mkdir container;
end
container_path = [path_start,'container\'];
cd(container_path)
save('container_path', 'container_path');
save('path_start', 'path_start');
sarray = [varargin{:}];
if contains(sarray,'green')
    kymo_mat_green = [path_start,'kymo_mat_green\'];
    save('kymo_mat_green', 'kymo_mat_green');
end
if contains(sarray,'red')
    kymo_mat_red = [path_start,'kymo_mat_red\'];
    save('kymo_mat_red', 'kymo_mat_red');
end
% change here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linetime_matstruct = [path_start,'linescan_time_mat\'];
save('linetime_matstruct', 'linetime_matstruct');
segmentation_dir = [path_start,'segmentation\'];
save('segmentation_dir', 'segmentation_dir');
MSD_path = [path_start,'MSD analysis\'];
save('MSD_path', 'MSD_path');
RESULTS = [path_start,'final\'];
save('RESULTS', 'RESULTS');
container_path = [path_start,'container\'];
save('container_path', 'container_path');
disp('Done');
end
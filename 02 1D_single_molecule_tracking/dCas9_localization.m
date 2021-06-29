% My goal is to perform 2D guassian fitting to locate position of dCas9 on
% the lambda DNA. This is using scan data. Eventually, I would like to
% correlate the FWHM for different z-positions, once I find where this
% information is stored in the metadata. 

my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\Projects\SWR1 Project\cas9\point spread functions for dCas9 data\20200310 red dCas9 cy5 specific binding\';
cd(my_directory)
path_start = strcat(pwd, '\');
listing = dir();
stringis = strcat(listing.name);
if not(contains(stringis, 'container'))
    mkdir container;
else
    rmdir container s;
    mkdir container;
end
container_path_dir = [path_start,'container\'];
cd(container_path_dir)
save('container_path_dir', 'container_path_dir');
save('path_start', 'path_start');
% scan_mat_green_dir = [path_start,'scan_mat_green\'];
% save('scan_mat_green_dir', 'scan_mat_green_dir');
    scan_mat_red_dir = [path_start,'scan_mat_red\'];
    save('scan_mat_red_dir', 'scan_mat_red_dir');
% change here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scantime_dir = [path_start,'scantime\'];
save('scantime_dir', 'scantime_dir');
pixeltime_dir = [path_start,'pixeltime\'];
save('pixeltime_dir', 'pixeltime_dir');
pixelsize_dir = [path_start,'pixelsize\'];
save('pixelsize_dir', 'pixelsize_dir');
disp('Done');
%%
% Initialize
scan_name = {}; 
scan_file_name = {};
file_names = {};
%Load matlab structures containing path information
mat = dir('*.mat'); 
for q = 1:length(mat)
    load(mat(q).name);
end
% red_scan = dir(scan_mat_red_dir);
% green_scan = dir(scan_mat_green_dir);
red_scan = dir(scan_mat_red_dir);
scantime = dir(scantime_dir);
pattern = [".", ".."];
counter = 1;
for i = 1:length(scantime) % Get name of kymos
    if not(startsWith(scantime(i).name, pattern))
    str = scantime(i).name;
    match = ["scantime_",".mat"];
    name = erase(str,match);
    scan_name{counter,1} = name;
    counter = counter +1;
    end
end

disp('....');
listing = dir();
stringis = strcat(listing.name);
if contains(stringis, 'green')
  cd(scan_mat_green_dir)
  counter = 1;
    for i = 1:length(green_scan) %Get file names of green kymos
        if not(startsWith(green_scan(i).name, pattern))
        name = green_scan(i).name;
        scan_file_name{counter,1} = name;
        counter = counter +1;
        end
    end
    for i = 1: length(scan_name) % Get kymo object array
        file_to_open = find(contains(scan_file_name,scan_name{i}));
        load (scan_file_name{file_to_open})
        scan_name{i,2} = double(green);
    end  
elseif contains(stringis, 'red')
    cd(scan_mat_red_dir)
    counter = 1;
    for i = 1:length(red_scan) %Get file names of green kymos
        if not(startsWith(red_scan(i).name, pattern))
        name = red_scan(i).name;
        scan_file_name{counter,1} = name;
        counter = counter +1;
        end
    end
    for i = 1: length(scan_name) % Get kymo object array
        file_to_open = find(contains(scan_file_name,scan_name{i}));
        load (scan_file_name{file_to_open})
        scan_name{i,2} = double(red);
    end  
end

disp('...');
cd(pixeltime_dir)
counter = 1;
pixeltime = dir(pixeltime_dir);
for i = 1:length(pixeltime) %Get file names of linetimes
    if not(startsWith(pixeltime(i).name, pattern))
    name = pixeltime(i).name;
    linetime_filename{counter,1} = name;
    counter = counter +1;
    end
end
for i = 1: length(scan_name) % Get linetimes
    file_to_open = find(contains(linetime_filename,scan_name{i}));
    load (linetime_filename{file_to_open})
    scan_name{i,3} = double(pixeltime);
end


cd(scantime_dir)
counter = 1;
for i = 1:length(scantime) %Get file names of linetimes
    if not(startsWith(scantime(i).name, pattern))
    name = scantime(i).name;
    linetime_filename{counter,1} = name;
    counter = counter +1;
    end
end
for i = 1: length(scan_name) % Get linetimes
    file_to_open = find(contains(linetime_filename,scan_name{i}));
    load (linetime_filename{file_to_open})
    scan_name{i,4} = double(scantime);
end
disp('..');


cd(pixelsize_dir)
pixelsize = dir(pixelsize_dir);

counter = 1;
for i = 1:length(pixelsize) %Get file names of linetimes
    if not(startsWith(pixelsize(i).name, pattern))
    name = pixelsize(i).name;
    linetime_filename{counter,1} = name;
    counter = counter +1;
    end
end
for i = 1: length(scan_name) % Get linetimes
    file_to_open = find(contains(linetime_filename,scan_name{i}));
    load (linetime_filename{file_to_open})
    scan_name{i,5} = double(pixelsize);
end
disp('.');

for i = 1:length(scan_name) % Save into a structured array
    scan_data(i).name = scan_name{i,1};
%     scan_data(i).green_scan = scan_name{i,2};
    scan_data(i).red_scan = scan_name{i,2};    
    scan_data(i).pixeltime = scan_name{i,3};
    scan_data(i).scantime = scan_name{i,4};
    scan_data(i).pixelsize = scan_name{i,5};
end
cd(container_path_dir)
save('scan_data.mat', 'scan_data');
disp('Done.');
%% Scan Data Segmenting Particles 
for i = 1:length(scan_data)
    shapeis = size(scan_data(i).red_scan);
    if length(shapeis) == 2
        scan = scan_data(i).red_scan;
            else
        scan = reshape(scan_data(i).red_scan(end-1,:,:),[shapeis(2),shapeis(3)]);
    end
        x = [];
        y = [];
        subscan = {};
        subcrop = {};
        % Finding centroids
        imagesc(scan);
        axis(gca, 'image')
        prompt = 'how many centroids?';
        num_particles = input(prompt);
        if num_particles>0
            [x,y] = ginput(num_particles);
            for j = 1:num_particles
                if scan_data(i).pixelsize == 100
                    [a,b] = imcrop(scan, [round(x(j))-5, round(y(j))-5, 10, 10]);
                else
                	[a,b] = imcrop(scan, [round(x(j))-15, round(y(j))-15, 30, 30]);
                end
                scan_data(i).subscan{j} = a;
                scan_data(i).subcrop{j} = b;
            end
            hold on
            for k = 1:num_particles
                rectangle('Position', scan_data(i).subcrop{k},'LineWidth',2);
                
            end
            hold off
             prompt = 'hold?';
             blank = input(prompt);
        else
        end
end
%% Fitting 2D gaussing to the cropped particles
for i = 1:length(scan_data)
    for j = 1:length(scan_data(i).subscan)
        if scan_data(i).pixelsize == 100
            MdataSize = 11; % Size of nxn data matrix
            [X,Y] = meshgrid(-5:5);
        else
            MdataSize = 31; % Size of nxn data matrix
            [X,Y] = meshgrid(-15:15);
        end
        x0(1) = max(max(scan_data(i).subscan{j})); %amplitude_guess
        x0(2) = 0; %startingx_guess
        x0(3) = 3; %sigmax_guess 
        x0(4) = 0; % startingy_guess
        x0(5) = 3; % sigmay_guess 
        x0(6) = 0; % rotation degrees
        InterpolationMethod = 'nearest'; % 'nearest','linear','spline','cubic'
        xdata = zeros(size(X,1),size(Y,2),2);
        xdata(:,:,1) = X;
        xdata(:,:,2) = Y;
        Z = scan_data(i).subscan{j}; %DATA HERE
        lb = [0,-MdataSize/2,0,-MdataSize/2,0,-pi/4];
        ub = [realmax('double'),MdataSize/2,(MdataSize/2)^2,MdataSize/2,(MdataSize/2)^2,pi/4];
        [x,resnorm,residual,exitflag] = lsqcurvefit(@D2GaussFunctionRot,x0,xdata,Z);
        scan_data(i).fit{j} = x;
    end
end

%%

count1 = 1;
count2 = 1;

for i = 4%1:length(scan_data)
    for j = 1:length(scan_data(i).subscan)
        if scan_data(i).pixelsize == 100
            xy_sigma_100{count1} = [round(scan_data(i).fit{j}(3)*100),round(scan_data(i).fit{j}(5)*100)];
            xy_FWHM_100(count1)= 2.355*((xy_sigma_100{count1}(1)+xy_sigma_100{count1}(2))/2);
            count1 = count1+1;
        else
            xy_sigma_30{count2} = [round(scan_data(i).fit{j}(3)*30),round(scan_data(i).fit{j}(5)*30)];
            xy_FWHM_30(count2) = 2.355*((xy_sigma_30{count2}(1)+xy_sigma_30{count2}(2))/2);
            count2 = count2+1;
        end
    end
end

%%
figure(1)
% histogram(xy_FWHM_100)
% figure(2)
% histogram(xy_FWHM_30)
% figure(3) 
histogram(all_FWHM)
title('FWHM Cy5-dCas9 Scans - panning through confocal plane')
xlabel('distance (nm)')
ylabel('counts')


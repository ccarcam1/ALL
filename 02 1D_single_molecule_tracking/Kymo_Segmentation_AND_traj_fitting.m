%% save the directory information here for the rest
path_start = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\Projects\SWR1 Project\sliding + Cas9 marker\2020-03-04-Cy5_dCas9_Cy3_SWR1\';
cd(path_start);
mkdir container;
cd('container');
save('path_start', 'path_start');
% change here
kymo_mat_green = [path_start,'kymo_mat_green\'];
save('kymo_mat_green', 'kymo_mat_green');
kymo_mat_red = [path_start,'kymo_mat_red\'];
save('kymo_mat_red', 'kymo_mat_red');
% change here
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
%% Extract and save kymograph information
cd(container_path)
clc
clear
kymoname = {}; % initialize
kymofilename = {};
file_names = {};
mat = dir('*.mat'); % matlab structures paths
for q = 1:length(mat)
    load(mat(q).name);
end
clear mat 
clear q
redkymo = dir(kymo_mat_red);
greenkymo = dir(kymo_mat_green);
linetime = dir(linetime_matstruct);
pattern = [".", ".."];
counter = 1;
for i = 1:length(linetime) % Get name of kymos
    if not(startsWith(linetime(i).name, pattern))
    str = linetime(i).name;
    match = [".mat"];
    name = erase(str,match);
%     name = strip(name,'left','_');
    name = name(1:(end-14));
    kymoname{counter,1} = name;
    counter = counter +1;
    end
end
cd(kymo_mat_green)
pattern = [".", ".."];
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
cd(kymo_mat_red)
pattern = [".", ".."];
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
cd(linetime_matstruct)
pattern = [".", ".."];
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

for i = 1:length(kymoname) % Save into a structured array
    data(i).name = kymoname{i,1};
    data(i).green_kymo = kymoname{i,2};
    data(i).red_kymo = kymoname{i,3};
    data(i).line_time = kymoname{i,4};
end
cd(container_path)
save('data.mat', 'data');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  segment the traces, defining the start of the trace and background
...first: click once top left of the trace of interest, click once bottom
...right of the trace of interest.
...Double click on a region that represents background intensity
...Double click at the start of the trace on a pixel that represents trace
...max intensity 
...The workspace for each trace will be saved, indicating which trace it is
...file names
% cd(container_path)
% clc
clear
mat = dir('*.mat'); % matlab structures paths
for q = 1:length(mat)
    load(mat(q).name);
end

for which = 1:length(data)
    kymograph = data(which).red_kymo;
    kymograph = rot90(kymograph, 3);
    imagesc(kymograph, [0,max(kymograph, [], 'all')/4]); % this is artificially bright to let you could the number of traces you want to process
    set(gcf,'Position',[25,25,600,800]); %laptop % set(gcf,'Position',[1150,50,700,1100]); %work computer
    ksize = size(kymograph);
    kw=ksize(:,2);
    kh=ksize(:,1);
    prompt = 'How many traces to extract?';
    num_of_traces = input(prompt);
    for j = 1:num_of_traces
        % select the region containing the trace of your particle of choice by
        % making two clicks: top left and bottom right corners of the future region
        % of interest
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
        close;
        crop_coordinates = [sp(2),sp(4),sp(1),sp(3)];
        figure;
          if max(kymograph,[], 'all') > 30
            imagesc(crop,[0,max(crop, [], 'all')/2])  
            set(gcf,'Position',[25,25,600,800]); %laptop % set(gcf,'Position',[1150,50,700,1100]); % work computer
            disp(['max intensity is ', num2str(max(kymograph,[], 'all'))]);
          else
            imagesc(crop);      
            set(gcf,'Position',[25,25,600,800]); %laptop %set(gcf,'Position',[1150,50,700,1100]); %work computer
          end

        ksize = size(crop);
        kw=ksize(:,2);
        kh=ksize(:,1);
        % choose the background (the region containing no  particles) by
        % making two clicks: top left and bottom right corners 
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

        dsd=2*std(BG); %% this gives you c
        dmin=d-(2*dsd);
        dmax=d+(2*dsd);
        close;

        figure;
        imagesc(crop,[0,max(crop, [], 'all')/2]);
%         set(gcf,'Position',[1150,50,700,1100]); %work computer
        set(gcf,'Position',[25,25,600,800]); %laptop
        % use 2 clicks to select the start of the particle trace; this will be used
        % later for the initial guess of the brightness and position of the
        % particle

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

        close;
        sz=size(crop);
        t=sz(:,1);
        w=sz(:,2);

        b=(sp(11)+sp(13))/2;
        bmin=b-5;
        bmax=b+5;

        wrange=[1:1:w];
        transpose(wrange);
        clf

        filename = [data(which).name, '_tracenum_', num2str(j)];
        cd(segmentation_dir);
        save(filename);
        if j ~= num_of_traces
            imagesc(kymograph,[0,max(kymograph, [], 'all')/2]);
%             set(gcf,'Position',[1150,50,700,1100]); %work computer
            set(gcf,'Position',[25,25,600,800]); %laptop
            ksize = size(kymograph);
            kw=ksize(:,2);
            kh=ksize(:,1);
        else
        end
disp(['this', num2str(j), 'out of', num2str(num_of_traces)])
    end
disp(['you are on trace ', num2str(which), ' out of ', num2str(length(data))])
tic
end

clc
clear
%% Fitting storing fitting results workspace 
% Initialize variables
% clc
% clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path_start = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\2019-11-20 200 mM KCl\';
load('path_start.mat')
load('kymo_matstruct.mat');
load('linetime_matstruct.mat');
load('segmentation_dir.mat');
input_dir2 = [path_start,'segmentation\'];
cd(input_dir2)

A = dir;
file_name_2 = {};
for i = 3: size(A,1)
    file_name_2{1,i-2} = erase(A(i).name,'.mat');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_name_2_location = [path_start,'filename\'];
cd(file_name_2_location)
% save('file_name_2', 'file_name_2')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_dir2 =  [path_start,'fitting\'];   
cd(file_name_2_location)
load('file_name_2')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_time = tic;
% for which2 = 1:length(file_name_2)
for which2 = 100

    %this can get commented out if you remembered to clear the workspace
    %before running the first section
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %     file_name_2_location = 'D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\C_Trap\2019-10-08 no ATP\file_name_2\';
    cd(file_name_2_location)
    load('file_name_2')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %     input_dir2 = 'D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\C_Trap\2019-10-08 no ATP\output\';
    cd(input_dir2)
%    splits the kymograph into lines and fits each line with a gaussian
    load(file_name_2{which2})
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=[1:t]
         subcrop=crop(i,:);
         transpose(subcrop);  
         fitresult = cell( 2, 1 );
         gof = struct( 'sse', cell( 2, 1 ), ...
             'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
         [xData, yData] = prepareCurveData( wrange, subcrop );
         ft = fittype( 'd + (a*exp(-((x-b)/c)^2))', 'independent', 'x', 'dependent', 'y' );
         opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
         opts.Display = 'Off';
         opts.Lower = [amin bmin -2 dmin];
         opts.StartPoint = [a b -1.3 d];
         opts.Upper = [amax bmax -1 dmax];

         [fitresult, gof] = fit( xData, yData, ft, opts );

         fitcoef=coeffvalues(fitresult);
         height(i)=fitcoef(:,1);
         width(i)=fitcoef(:,3);
         pixel(i)=fitcoef(:,2);
         timepix(i)=i;
        intensity(i)= fitcoef(:,1)*sqrt(2*pi*(fitcoef(:,3)^2));

         a=fitcoef(:,1);
         amin=a-dsd;
         amax=a+dsd;

         b=fitcoef(:,2);
         bmin=b-5;
         bmax=b+5;
    disp([num2str(i), '/', num2str(t)]);
    % output_dir2 =  'D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\C_Trap\2019-10-08 no ATP\output2\';   
    cd(output_dir2);
% % % %     filename = (file_name_2{which2});
% % % %     save(filename);
gaussfitting(i).lineintensities = [xData, yData];
gaussfitting(i).gaussfit = fitresult;
    end
disp(['you are on trace ', num2str(which2), ' out of ', num2str(length(file_name_2))])
end
total_time_save = toc(total_time)

%% Experiments making a video of gauss fits
arrayme ={};
for i = 1:length(gaussfitting)
    plot(gaussfitting(i).lineintensities(:,1), gaussfitting(i).lineintensities(:,2))
    hold on
    plot(gaussfitting(i).gaussfit, gaussfitting(i).lineintensities(:,1), gaussfitting(i).lineintensities(:,2))
    xlabel('pixels')
    ylabel('intensity (counts)')
    F = getframe;
    arrayme{i} = F.cdata;
    hold off
    disp(num2str(i))
end

%% Experiments part 2 making a video of gauss fits
videofile = VideoWriter('gaussmovie.avi')
open(videofile)
for i = 1:length(arrayme)
   writeVideo(videofile, imcomplement(arrayme{i}))
end
close(videofile)

%% file name 2 merged with filename
clc
clear
load('path_start.mat')
% path_start = 'D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\2019-11-20 200 mM KCl\';
kymostructdir = [path_start,'kymo_mat\'];
scantimedir = [path_start,'linescan_time_mat\'];
fitdir = [path_start,'fitting\'];
file_name_2_location = [path_start,'filename\'];


tracenames = dir(fitdir);
kymonames = dir(kymostructdir);
scantimenames = dir(scantimedir);
kymonames_cell = {};
file_names = {};

for i = 3:length(kymonames)
kymonames_cell{1,i-2} = erase(kymonames(i).name,'.mat');
scantime_cell{1,i-2} = erase(scantimenames(i).name,'.mat');
end

cd(scantimedir)
for i = 1:length(scantime_cell)
    load(scantime_cell{1,i});
    scantime_cell{2,i} = linescan_time_;
end
    
for i = 3: length(tracenames)
    file_names{1,i-2} = tracenames(i).name;
end

for i = 1:length(kymonames_cell)
    for j = 1:length (scantime_cell)
        if contains(scantime_cell{1,j},kymonames_cell{1,i}) 
            kymonames_cell{2,i} = scantime_cell{2,j};
        else
        end
    end
end

for i = 1:length(file_names)
    for j = 1:length(kymonames_cell)
        if contains(file_names{1,i}, kymonames_cell{1,j})
            file_names{2,i} = kymonames_cell{2,j};
        end
    end
end

cd(file_name_2_location)
save('file_names', 'file_names')

%% select traces to save and which ones to repeat analysis 
%initialize
clc
clear
% load('path_start.mat')
path_start = 'D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\2019-11-25 25 mM KCl\';
% path_start = 'D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\2019-10-23 SWR1 gamma S\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_dir3 =  [path_start,'fitting\'];   
final_results = [path_start,'final\'];
file_name_2_location = [path_start,'filename\'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cd(file_name_2_location);
load('file_names');
keep_these = cell(5,length(file_names));
not_gonna_fit = cell(4,length(file_names));
could_be_better = cell(4,length(file_names));
multiple_molecules = cell(4,length(file_names));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for which3 = 1:length(file_names)
    path_start = 'D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\2019-11-25 25 mM KCl\';
    cd([path_start,'filename\']);
    load('file_names');
    file = file_names{1,which3};
    timestep = file_names{2,which3};
    cd(input_dir3)
    load(file_names{1,which3})
    filename = strrep(filename,'_',' ');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    int=abs(height.*width);
    int=transpose(int);
    intensity = transpose(intensity);
    time=timepix*timestep;
    time=transpose(time);
    position=pixel.*pix;
    position=transpose(position);
    % plot coordinates (blue) and brightness (green) vs time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    yyaxis left
    plot(time, position)
    ylabel('position, \mum');
    yyaxis right
    plot(time, int)
%     set(gcf,'Position',[0,75,700,500]);% for work computer
    set(gcf,'Position',[0,25,400,300]);% for laptop
    title(filename);
    xlabel('time, s');
    ylabel('intensity, a.u.');
    % displays the original cropped image side by side with the traced
    % coordinates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
%     set(gcf,'Position',[650,50,700,1100]); %for work computer
    set(gcf,'Position',[375,25,500,600]);% for laptop
    imagesc(kymograph);
    hold on
    ksize = size(kymograph);
    kw=ksize(:,2);
    kh=ksize(:,1);
    
    rectangle('Position',[crop_coordinates(3),crop_coordinates(1),...
        crop_coordinates(4)-crop_coordinates(3),crop_coordinates(2)-...
        crop_coordinates(1)], 'EdgeColor','y'); 
    hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
%     set(gcf,'Position',[1350,50,700,1100]); % for work computer
    set(gcf,'Position',[825,25,450,600]); % for laptop
    title(filename);
    s1=subplot(1,3,1);imagesc(crop),colormap(gray);
    title('kymograph');
    xlabel('x, px');
    ylabel('y, px');
    s2=subplot(1,3,2);imagesc(crop),colormap(gray);
    hold on
    plot(s2,pixel,timepix,'-y');
    title('overlay');
     xlabel('x, px');
    ylabel('y, px');
    hold off
    s3=subplot(1,3,3);
    plot(s3,pixel,timepix,'ro');
    s3=gca;
    xlim(gca,[0 w]);
    ylim(gca,[0 t]);
    set(gca,'YDir','Reverse');
    title('trace');
    xlabel('x, px');
    ylabel('y, px');
    out=horzcat(time,position,intensity,int);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     choose if you want to save or discard the fitted coordinates
choice = menu('Choose an action','Good Fit','Could-be-better Fit',...
    'multiple molecules', 'Not going to fit, likely', 'Break');
disp(['you are on trace ', num2str(which3), ' out of ', num2str(length(file_name_2))]);
if choice == 1
    keep_these{1,which3} = filename;
    keep_these{2,which3} = out;
    keep_these{3,which3} = crop_coordinates;
    keep_these{4,which3} = crop;
    keep_these{5,which3} = timestep;
    keep_these{6,which3} = kymograph;
elseif choice == 2
    could_be_better{1,which3} = filename;
    could_be_better{2,which3} = out;
    could_be_better{3,which3} = crop_coordinates;
    could_be_better{4,which3} = crop;
    could_be_better{5,which3} = timestep;
elseif choice == 3
    multiple_molecules{1,which3} = filename;
    multiple_molecules{2,which3} = out;
    multiple_molecules{3,which3} = crop_coordinates;
    multiple_molecules{4,which3} = crop;
    multiple_molecules{5,which3} = timestep;
elseif choice == 4
    not_gonna_fit{1,which3} = filename;
    not_gonna_fit{2,which3} = out;
    not_gonna_fit{3,which3} = crop_coordinates;
    not_gonna_fit{4,which3} = crop;
    not_gonna_fit{5,which3} = timestep;
else
    break
end
close all

   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(final_results);
save('keep_these', 'keep_these');
save('not_gonna_fit', 'not_gonna_fit');
save('could_be_better','could_be_better');
save( 'multiple_molecules', 'multiple_molecules');
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

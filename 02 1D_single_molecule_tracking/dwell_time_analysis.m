x = 1;
for i = 1:length(saltdata)
    if saltdata(i).r2point8 == 1
        lengthsalt(x) = length(saltdata(i).frames);
        x = x+1;
    else
    end
end
%%
lowsalt = lengthsalt(saltsalt ==25);
mediumsalt = lengthsalt(saltsalt ==70);
highsalt = lengthsalt (saltsalt == 200);
%%
figure(1)
subplot(3,1,1)
histogram(lowsalt,20,'FaceColor', [ 0.6350    0.0780    0.1840],'BinLimits',[0, 1200])
hold on 
xline(median(lowsalt))
hold off
subplot(3,1,2)
histogram(mediumsalt,20,'FaceColor',[ 0    0.4470    0.7410],'BinLimits',[0, 1200])
hold on 
xline(median(mediumsalt))
hold off
subplot(3,1,3)
histogram(highsalt,20,'FaceColor', [0.4940    0.1840    0.5560],'BinLimits',[0, 1200])
hold on
xline(median(highsalt))
hold off
%%
figure(2)
plot([25,70,200],[median(lowsalt),median(mediumsalt),median(highsalt)])
hold on
plot([25,70,200],[mean(lowsalt),mean(mediumsalt),mean(highsalt)])
hold off
%% Dwell Time Analysis
clc
clear
%matlab kymograph structures
path_start = ('C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\salt concentrations ATP\2019-11-25 25 mM KCl\');
kymo_matstruct = ('C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\salt concentrations ATP\2019-11-25 25 mM KCl\kymo_mat\');
linetime_matstruct = ('C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\salt concentrations ATP\2019-11-25 25 mM KCl\linescan_time_mat\');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_dir = kymo_matstruct;
cd(kymo_matstruct)
A = dir(kymo_matstruct);
cd(linetime_matstruct)
B = dir;
bname = {};
for i = 1:length(B)
bname{i} = B(i).name;
end
file_names = {};
cd(kymo_matstruct)
for i = 3: size(dir,1)
    file_names{1,i-2} = erase(A(i).name,'.mat');
end
cd(linetime_matstruct)
for i = 3: size(dir,1)
    TF = find(contains(bname,[file_names{1,i-2},'_']));
    load (bname{TF})
    file_names{2,i-2} = linescan_time_;
end
A = dir(kymo_matstruct);
% file_name = {};
for i = 3: size(dir,1)
    file_names{1,i-2} = erase(A(i).name,[".mat", "_16_bit"]);
    file_names{3,i-2} = erase(A(i).name,'.mat');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x= 1;
for i = 1:length(file_names)
    kymo_matstruct = input_dir;
    cd(kymo_matstruct)
    dim = 0.1; %um
    file = file_names{3,i};
    timestep = file_names{2,i}; %seconds
    loadpath=strcat(kymo_matstruct,file,'.mat');
    kymograph=load(loadpath);
    kymograph = kymograph.obj_arr;
    kymograph = rot90(kymograph, 3);
    kymograph = double(kymograph);
    imagesc(kymograph, [0,max(kymograph, [], 'all')/4]); % this is artificially bright to let you could the number of traces you want to process
    set(gcf,'Position',[25,25,600,800]); %laptop
    ksize = size(kymograph);
    kw=ksize(:,2);
    kh=ksize(:,1);
    prompt = 'mark another?';
    num_of_traces = input(prompt);
    y =1;
        while num_of_traces == 1
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
            dwelltime(x).name = [file_names{1,i}, '_tracenum_', num2str(y)];
            imagesc(kymograph, [0,max(kymograph, [], 'all')/4]);
            set(gcf,'Position',[25,25,600,800]); %laptop
            ksize = size(crop);
            kw=ksize(:,2);
            dwelltime(x).kh=ksize(:,1);
            y = y+1;
            x = x +1;
            num_of_traces = input(prompt);
        end
    disp(['you are on trace ', num2str(i), ' out of ', num2str(length(file_names))])
end

%% Use this one 

%matlab kymograph structures
path_start = ('D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\salt concentrations ATP\2019-11-20 200 mM KCl\');
kymo_matstruct = ('D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\salt concentrations ATP\2019-11-20 200 mM KCl\kymo_mat\');
linetime_matstruct = ('D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\salt concentrations ATP\2019-11-20 200 mM KCl\linescan_time_mat\');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_dir = kymo_matstruct;
cd(kymo_matstruct)
A = dir(kymo_matstruct);
cd(linetime_matstruct)
B = dir;
bname = {};
for i = 1:length(B)
bname{i} = B(i).name;
end
file_names = {};
cd(kymo_matstruct)
for i = 3: size(dir,1)
    file_names{1,i-2} = erase(A(i).name,'.mat');
end
cd(linetime_matstruct)
for i = 3: size(dir,1)
    TF = find(contains(bname,[file_names{1,i-2},'_']));
    load (bname{TF})
    file_names{2,i-2} = linescan_time_;
end
A = dir(kymo_matstruct);
% file_name = {};
for i = 3: size(dir,1)
    file_names{1,i-2} = erase(A(i).name,[".mat", "_16_bit"]);
    file_names{3,i-2} = erase(A(i).name,'.mat');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x= 1; % comment this out if you need to break midway through analysis otherwise will overwrite previous analysis

choice = 1;
for i = 52:length(file_names)
    kymo_matstruct = input_dir;
    cd(kymo_matstruct)
    dim = 0.1; %um
    file = file_names{3,i};
    timestep = file_names{2,i}; %seconds
    loadpath=strcat(kymo_matstruct,file,'.mat');
    kymograph=load(loadpath);
    kymograph = kymograph.obj_arr;
    kymograph = rot90(kymograph, 3);
    kymograph = double(kymograph);
    imagesc(kymograph, [0,max(kymograph, [], 'all')/4]); % this is artificially bright to let you could the number of traces you want to process
    set(gcf,'Position',[25,25,600,800]); %laptop
    ksize = size(kymograph);
    kw=ksize(:,2);
    kh=ksize(:,1);
%     prompt = 'mark another?';
%     num_of_traces = input(prompt);
    y =1;
    if choice ==3 
        break
    end
    choice = menu('Choose an action','mark another','next','Break');
        while choice == 1
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
            dwelltime(x).name = [file_names{1,i}, '_tracenum_', num2str(y)];
            imagesc(kymograph, [0,max(kymograph, [], 'all')/4]);
            set(gcf,'Position',[25,25,600,800]); %laptop
            ksize = size(crop);
            kw=ksize(:,2);
            dwelltime(x).kh=ksize(:,1);
            choice2 = menu('Did it bind after kymo start','yes','no');
            if choice2 == 1
                dwelltime(x).bindafterstart=1;
            elseif choice2 ==2 
                dwelltime(x).bindafterstart=0;
            end
            choice = menu('Choose an action','mark another','next','Break');
            y = y +1;
            x = x +1;
            if choice ==3 
                break
            end
        end
    disp(['you are on trace ', num2str(i), ' out of ', num2str(length(file_names))])
end

close all
%%
for i = 1:length(dwelltime200)
    bindingwhen200(i) = dwelltime200(i).bindafterstart;
    lenghttime200(i) = dwelltime200(i).seconds;
end
%%
afterstart200=lenghttime200(bindingwhen200==1);
withstart200 =lenghttime200(bindingwhen200==0);
%%
% position
% [1762.142857142857,143.2857142857143,431.4285714285713,1009.714285714286]
binlimits = [0, 85];
figure(1)
subplot(3,1,1)
% histogram(all25, 30,'FaceColor', [ 0.6350    0.0780    0.1840],'BinLimits',binlimits,'Normalization', 'probability','LineWidth', 1)
% hold on
% histogram(withstart25,20,'FaceColor', [ 0.6350    0.0780    0.1840],'BinLimits',[0, 1650],'Normalization', 'probability')
histogram(withstart25, 30,'FaceColor', [ 0.6350    0.0780    0.1840],'BinLimits',binlimits,'Normalization', 'probability','LineWidth', 1)

% histogram(afterstart25, 30,'FaceColor', [ 0.6350    0.0780    0.1840],'BinLimits',binlimits,'Normalization', 'probability','LineWidth', 1)
% hold off

title('dwelltime 25mM KCl')
legend(['with start 25mM n = ', num2str(length(withstart25))]);
% legend(['all 25mM n = ', num2str(length(all25))],['with start 25mM n = ', num2str(length(afterstart25))]);
% legend('with start 25mM','after start 25mM');
% legend('after start 25mM');
xlabel('time (sec)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,2)
% histogram(all70, 30,'FaceColor', [0    0.4470    0.7410],'BinLimits',binlimits,'Normalization', 'probability','LineWidth', 1)
% hold on
% histogram(withstart70,20,'FaceColor', [ 0    0.4470    0.7410],'BinLimits',[0, 1650],'Normalization', 'probability')
histogram(withstart70, 30,'FaceColor', [0    0.4470    0.7410],'BinLimits',binlimits,'Normalization', 'probability','LineWidth', 1)

% histogram(afterstart70, 30,'FaceColor', [0    0.4470    0.7410],'BinLimits',binlimits,'Normalization', 'probability','LineWidth', 1)
% hold off
title('dwelltime 70mM KCl')
% legend('with start 70mM', 'after start 70mM');
% legend('after start 70mM');
xlabel('time (sec)')
legend(['with start 70mM n = ', num2str(length(withstart70))]);

% legend(['all 70mM n = ', num2str(length(all70))],['with start 70mM n = ', num2str(length(afterstart70))]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,3)
% histogram(all200, 30,'FaceColor', [ 0.4940    0.1840    0.5560],'BinLimits',binlimits,'Normalization', 'probability','LineWidth', 1)
% hold on
% histogram(withstart200,20,'FaceColor', [ 0.9290    0.6940    0.1250],'BinLimits',[0, 1650],'Normalization', 'probability')
histogram(withstart200, 30,'FaceColor', [ 0.4940    0.1840    0.5560],'BinLimits',binlimits,'Normalization', 'probability','LineWidth', 1)

% histogram(afterstart200, 30,'FaceColor', [ 0.4940    0.1840    0.5560],'BinLimits',binlimits,'Normalization', 'probability','LineWidth', 1)
% hold off

title('dwelltime 200mM KCl')
% legend(['all 200mM n = ', num2str(length(all200))],['with start 200mM n = ', num2str(length(afterstart200))])
% legend('with start 200mM', 'after start 200mM');
% legend('after start 200mM');
legend(['with start 200mM n = ', num2str(length(withstart200))])

xlabel('time (sec)')


%%
figure(2)
binlimits = [0, 85];
histogram(all25, 30,'FaceColor', [ 0.6350    0.0780    0.1840],'BinLimits',binlimits,'Normalization', 'probability','LineWidth', 1)
hold on
histogram(all70, 30,'FaceColor', [0    0.4470    0.7410],'BinLimits',binlimits,'Normalization', 'probability','LineWidth', 1)
histogram(all200, 30,'FaceColor', [ 0.4940    0.1840    0.5560],'BinLimits',binlimits,'Normalization', 'probability','LineWidth', 1)

hold off

title('dwelltime all measured')
% legend(['all 25mM n = ', num2str(length(all25))],['all 70mM n = ', num2str(length(all70))],['all 200mM n = ', num2str(length(all200))]);
legend(['all 25mM n = ', num2str(length(all25))]);
xlabel('time (sec)')


%%

a = length(withstart25);
b = length(afterstart25);

c = length(withstart70);
d = length(afterstart70);

e = length(withstart200);
f = length(afterstart200);

plot(a,b)
hold on
plot(c,d)
hold on
plot(e,f)
hold off

%% correct blunder 

for i = 1:length(mydata_25mMKCl_noNaN)
    for j = 1:length(dwelltime25)
        l = ["kymo_", "_tracenum_1", "_tracenum_2", "_tracenum_3", "_tracenum_4", "_tracenum_5", "_tracenum_6", "_tracenum_7", "_tracenum_8","_tracenum_9", "_tracenum_10", "_tracenum_11", "_tracenum_12", "_tracenum_13", "_tracenum_14", "_tracenum_15", "_tracenum_16", "_tracenum_17"];
        k = erase(dwelltime25(j).name,l);
        if contains(mydata_25mMKCl_noNaN(i).mol_id, k)
            dwelltime25(j).scanrate = mydata_25mMKCl_noNaN(i).timestep;
        end
    end
end

%% further refine data inadvertent save

for i = 1:length(dwelltime200)
    if isempty(dwelltime200(i).scanrate)
        dwelltime200(i) = [];
    else
    end
end

%% correct blunder part two

for i = 1:length(dwelltime200)
    dwelltime200(i).seconds = (dwelltime200(i).kh * dwelltime200(i).scanrate)/1000;
end

%% total lengths

all200 = [withstart200 afterstart200];
all25 = [withstart25 afterstart25];
all70 = [withstart70 afterstart70];



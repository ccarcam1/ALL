%% Master 1D Single Molecule Tracking Script
% Claudia Carcamo 
% 03 - 14 - 2020 
%% save the directory information here for the rest
% REQUIREMENTS
    % FUNCTION: my_directory_function('color',)
    % Must Start In Directory Of Interest
    % Directory must be made in jupyter notebook script
    
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\Swc2\2020-12-08 Swc2\';
cd(my_directory)
% my_directory_function('green');

%For two color analysis
% my_directory_function('green red');
my_directory_function('green');


%% Extract and save kymograph information
% REQUIREMENTS
    % FUNCTION: my_kymodata_structure --> no input arguments
    % Must Start In "Container" directory
    % Final data structure includes all colors 
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\Swc2\2020-12-08 Swc2\';
cd(my_directory)
cd container
% my_kymodata_structure2('green red')
my_kymodata_structure2('green')
%% segment particles from different color channels
% REQUIREMENTS
    % FUNCTION: structure_name = my_segment_kymos_improvedwithzoom2(start_var, end_var, which_color)
    % Must Start In "Container" directory
    % linescan_time_mat
    % kymo_mat_green
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\Swc2\2020-12-08 Swc2\';
cd(my_directory)
cd([my_directory, 'container'])
load('data.mat')
length(data)
green_segmentation1 = my_segment_kymos_improvedwithzoom2(1 ,10, 'green'); % length(data)
%% Combine separate strcutures 
a=0;% select this for the first one
% a=length(green_segmentation_ALL);% for the first one, don't include this.
for i = 1:length(structure_name)
    green_segmentation_ALL(i+a) = structure_name(i);
end
disp("Done")


%% get rid of crops that are too small to correctly get gaussian fits
counter= 1;
for i = 1:length(green_segmentation_ALL)
    if size(green_segmentation_ALL(i).crop,1)<10
    recordthis(counter) = i ;
    counter = counter +1;
    end
end

disp("Done")
%% Fit gaussians to particles over time
% REQUIREMENTS
    % FUNCTION gaussfitting = my_gaussian_fitting(x,y,segmented_kymos)
    % x = starting point
    % y = ending (can be length(segmented_kymos))
    % segmented_kymos --> results from previous section 
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\sliding\salt concentrations ATP\2020-12-20 SWR1 150mM KCl 1mM ATP no gloxy\';
% % cd(my_directory);
green_gaussfitting = my_gaussian_fitting(1,length(green_segmentation_ALL),green_segmentation_ALL);
cd([my_directory, 'fitting']);
save('green_gaussfitting_second_try.mat', 'green_gaussfitting')

% cd(my_directory);
% cd([my_directory, 'segmentation']);
% load('red_segmentation.mat');
% cd(my_directory)
% cd([my_directory, 'fitting']);
% red_gaussfitting = my_gaussian_fitting(1,length(structure_name),structure_name);
% save('red_gaussfitting.mat', 'red_gaussfitting')
%% Determine Which Fits look good by eye (part 1 organize previous variables)
% REQUIREMENTS
    % FUNCTION keep_these_structure = my_visualize_fits(start_var, end_var, gaussfitting, data, structure_name)
% clear
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\sliding\salt concentrations ATP\2020-12-20 SWR1 150mM KCl 1mM ATP no gloxy\';
% cd([my_directory, 'container']);
% load('data.mat')
% cd([my_directory, 'segmentation']);
% load('red_segmentation_ALL.mat')
% cd([my_directory, 'fitting']);
% load('red_gaussfitting.mat')
cd([my_directory, 'MSD analysis']);
green_fitting_MSD_structure = my_condense_relevant_info('green', data, green_gaussfitting, green_segmentation_ALL);
save('green_fitting_MSD_structure_second.mat', 'green_fitting_MSD_structure')
% cd([my_directory, 'MSD analysis']);
% red_fitting_MSD_structure = my_condense_relevant_info('red', data, red_gaussfitting, structure_name);
% save('red_fitting_MSD_structure.mat', 'red_fitting_MSD_structure')

%% Determine Which Fits look good by eye (part 2 visualize traces and decide which to keep)
% function keep_these_structure = my_visualize_fits(start_var, end_var, pre_MSD)
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\sliding\salt concentrations ATP\2020-12-20 SWR1 150mM KCl 1mM ATP no gloxy\';
cd([my_directory, 'MSD analysis']);
% load('green_fitting_MSD_structure.mat');
disp('..................')
startis =301;
stopis = length(green_fitting_MSD_structure); 
green_keep_these = my_visualize_fits(startis, stopis, green_fitting_MSD_structure);
save(['green_keep_these1 ', num2str(startis), ' to ', num2str(stopis),'.mat'], 'green_keep_these')



% red_keep_these = my_visualize_fits(startis, stopis, red_fitting_MSD_structure);
% save(['red_keep_these ', num2str(startis), ' to ', num2str(stopis),'.mat'], 'red_keep_these')

% % Combine saved structures into one:
% a=length(green_keep_these_structure);
% for i = 1:length(keep_these_structure)
%     green_keep_these_structure(i+a) = keep_these_structure(i);
% end

% % my_visualize_fits2(1, 1, green_fitting_MSD_structure);

%%
% a=0;% select this for the first one
a=length(green_keep_these_All);% for the first one, don't include this.
for i = 1:length(green_keep_these)
    green_keep_these_All(i+a) = green_keep_these(i);
end
disp("Done")

%%
a=0;% select this for the first one
% a=length(green_keep_these_All);% for the first one, don't include this.
for i = 1:length(red_keep_these)
    red_keep_these_All(i+a) = red_keep_these(i);
end
disp("Done")
%% get rid of the blank cells get the MY_DATA structure
counter = 1;
for i = 1:length(green_keep_these_All)
    if ~isempty(green_keep_these_All(i).name)
        MY_DATA(counter)= green_keep_these_All(i);
        counter= counter+1;
    end
end
disp('DONE')

% counter = 1;
% for i = 1:length(red_keep_these_All)
%     if ~isempty(red_keep_these_All(i).name)
%         MY_DATA(counter)= red_keep_these_All(i);
%         counter= counter+1;
%     end
% end
% disp('DONE')

% my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data_Analysis\Raw_Data\2020_09_02_dCas9_SWR1_twocolorsuccess\';
% cd([my_directory, 'MSD analysis']);
% save('Red_MY_DATA.mat', 'MY_DATA')

%% MSD calculation
% %function MSDs = my_MSD_calculator(structure_with_data)
% clc 
% clear
% my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data_Analysis\Projects\SWR1 Project\sliding\ADP\2020-09-22 SWR1 ADP\';
% cd([my_directory, 'MSD analysis']);
% load('Red_MY_DATA.mat');

disp('New Run');
MSDs = my_MSD_calculator2(MY_DATA);
% save('Red_MSDs.mat', 'MSDs')
save('Green_MSDs.mat', 'MSDs')

[rsquared,pvalue,yintercept,slopeval,fit_cutoffis] = myCutOffFinder(MSDs);
for i = 1:length(MY_DATA)
    MY_DATA(i).rsquared = rsquared{i};
    MY_DATA(i).pvalue = pvalue{i};
    MY_DATA(i).yintercept = yintercept{i};
    MY_DATA(i).slopeval = slopeval{i};
    MY_DATA(i).fit_cutoffis = fit_cutoffis(i);
    MY_DATA(i).MSD = MSDs(i).MSD;
    if isnan(fit_cutoffis(i)) || fit_cutoffis(i)==0 % if NaN take first 25% to be fit cut_off
        MY_DATA(i).slope = slopeval{i}(length(slopeval{i}/4));
        MY_DATA(i).r2 = rsquared{i}(length(rsquared{i}/4));
        MY_DATA(i).yinter = yintercept{i}(length(yintercept{i}/4));
        MY_DATA(i).pval = pvalue{i}(length(pvalue{i}/4));
    else
        MY_DATA(i).slope = slopeval{i}(fit_cutoffis(i));
        MY_DATA(i).r2 = rsquared{i}(fit_cutoffis(i));
        MY_DATA(i).yinter = yintercept{i}(length(yintercept{i}/4));
        MY_DATA(i).pval = pvalue{i}(length(pvalue{i}/4));
    end
end
% save('Red_MY_DATA_MSDs.mat', 'MY_DATA')
save('Green_MY_DATA_MSDs.mat', 'MY_DATA')

%% Processing 1D diffusion 
% I will consider anything shorter than 0.75 seconds to be short
% 9-12-2020 I will use 0 seconds as the threshold for red because I want to
% maintain even the shortest traces AND for green traces 

for i = 1:length(MY_DATA)
    length_Seconds(i) = size(MY_DATA(i).particle_tracked(:,2), 1)*(MY_DATA(i).line_time/1000);
end
% 
time_cutoff_seconds = 0;
lookhere = find(length_Seconds < time_cutoff_seconds);
x = 1;
y = 1;
for i = 1:length(MY_DATA)
    if length_Seconds(i) >time_cutoff_seconds
       tempMY_DATA(x) = MY_DATA(i);
       x = x +1;
    else
        shortMY_DATA(y) = MY_DATA(i);
       y = y +1;
    end
end

clear MY_DATA
MY_DATA = tempMY_DATA;
clear tempMY_DATA
% clear shortMY_DATA
clear x
clear y
clear i
clear lookhere
clear length_Seconds
clear time_cutoff_seconds

%% Process NaN values


NaN_negative_which = [];
NaN_which = [];
NaN_positive_which = [];
negative_which_noNaN = [];

x = 1;
y = 1;
z = 1;
a = 1;
for k = 1: length(MY_DATA)
    if isnan(MY_DATA(k).fit_cutoffis)
        NaN_which(x) = k;
        x = x+1;
        if MY_DATA(k).slope<0
            NaN_negative_which(y) = k;
            y = y+1;
        else
            NaN_positive_which(a) = k;
            a = a+1;
        end
    elseif MY_DATA(k).slope<0
        negative_which_noNaN(z) = k;
        z = z+1;
    end  
end

clear x
clear y
clear z
clear a

%% Finalize Data Structures 
% what to make into 0, these aren't moving
% negative_which_noNaN
% and
% NaN_negative_which

for i = 1:length(negative_which_noNaN)
    k = negative_which_noNaN(i);
    MY_DATA(k).slope = 0;
end
for i = 1:length(NaN_negative_which)
    k = NaN_negative_which(i);
    MY_DATA(k).slope = 0;
end
Green_MY_DATA = MY_DATA;
save('Green_MY_DATA_final_structure.mat', 'Green_MY_DATA')
% save('Red_MY_DATA_final_structure.mat', 'MY_DATA')

%% Get rid of poor fitting data
% MY_DATA = Red_MY_DATA;
cut_off_r2 = 0.8;
slope = [];
diffusioncoeff = [];
intercept = [];
r2 = [];
x = 1;
for i = 1:length(MY_DATA)
        slope(i) = MY_DATA(i).slope;
        r2(i) = MY_DATA(i).r2;
        intercept(i) = MY_DATA(i).yinter;
        if MY_DATA(i).r2>cut_off_r2
            MY_DATA(i).r2point8 = 1;
            x = x+1;
        else 
            MY_DATA(i).r2point8 = 0;
        end
end
diffusioncoeff= slope/2;
Dr2less08 = diffusioncoeff(r2<cut_off_r2);
Dr2greater08 = diffusioncoeff(r2>cut_off_r2);

save('Green_diffusioncoeff', 'diffusioncoeff');
save('Green_Dr2less08', 'Dr2less08');
save('Green_Dr2greater08','Dr2greater08');
% 
% save('Red_diffusioncoeff', 'diffusioncoeff');
% save('Red_Dr2less08', 'Dr2less08');
% save('Red_Dr2greater08','Dr2greater08');

%% TEXT BOX FOR HISTOGRAMS
histogram(Dr2greater08, 20)
% hold on
% histogram(redDr2greater08)
% hold off
% dim = [0.7 0.6 0.3 0.3];
% str = ['SWR1 n= ', num2str(length(greenDr2greater08))];%, newline, 'dCas9 n= ', num2str(length(redDr2greater08))];
% annotation('textbox',dim,'String',str,'FitBoxToText','on');
clf
hold on
bins = 20;
upper = .3;
histogram(Dr2greater08_old,bins, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 upper], 'DisplayName',strcat("old SWR1 150mM n = ", num2str(length(Dr2greater08_old))))

histogram(Dr2greater08,bins, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 upper], 'DisplayName',strcat("SWR1 150mM n = ", num2str(length(Dr2greater08))))
legend



%%
counter = 1;
mM70_salt_gammaS = [];
for i = 1:length(gammaS)
%     if gammaS(i).salt == 70
        mM70_salt_gammaS(counter) = gammaS(i).slope/2;
        counter = counter +1;
%     end
end

%%
binlims = [0 0.5];

histogram(Dr2greater08, 30,'BinLimits', binlims,'FaceColor','#7E2F8E', 'Normalization', 'probability')
hold on
histogram(mM70_salt_ATP, 30,'BinLimits', binlims,'FaceColor','#0072BD', 'Normalization', 'probability')
histogram(ATP_diffusion_all2020, 30, 'BinLimits', binlims,'FaceColor','none', 'Normalization', 'probability', 'LineWidth', 1.5)
hold on
% histogram(mM70_salt_gammaS, 30, 'BinLimits', binlims,'FaceColor','black','Normalization', 'probability')
hold off

%%
% % Visualize NaNs to determine which to keep
% 
% for i = 1:length(NaN_which)
% %     k = lookhere(i);
% % for i = 1:length(NaN_negative_which)
%     k = NaN_which(i);
%     if ~isnan(MY_DATA(k).fit_cutoffis)
%         figure
%         set(gcf,'Position',[60,30,1400,1100]);
%         subplot(2,2,1);
%         MSD = MY_DATA(k).MSD;
%         n=size(MY_DATA(k).particle_tracked(:,2), 1);
%         errorbar(MSD(1:n-1, 1), MSD(1:n-1, 2), MSD(1:n-1, 3), 'ok')
%         hold on
%         threshold = MY_DATA(k).fit_cutoffis;
%         x = MSD(1:threshold, 1);
%         y = MSD(1:threshold, 2);
%         mdl = fitlm(x,y);
%         plot(mdl);
%         slope = mdl.Coefficients.Estimate(2);
%         m = MY_DATA(k).name;
%         title(['Best linear fit for ', m],'FontSize', 14)
%         xlabel('time (s)', 'FontSize', 14)
%         ylabel('MSD (um^2)', 'FontSize', 14)
%         legend('off')
%         save_name  = [m,'.png'];
%         hold off
% %         saveas(gcf, ['full length ', save_name])
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         subplot(2,2,3);
%         errorbar(MSD(1:threshold, 1), MSD(1:threshold, 2), MSD(1:threshold, 3), 'ok')
%         hold on
%         mdl = fitlm(x,y);
%         plot(mdl)
%         title(['Best linear fit for ', m],'FontSize', 14)
%         xlabel('time (s)', 'FontSize', 14)
%         ylabel('MSD (um^2)', 'FontSize', 14)
%         annotation('textbox',[.13 .1 .1 .35],'String',['R^2= ', num2str(MY_DATA(k).rsquared(threshold))...
%             , newline, 'p-value= ', num2str(MY_DATA(k).pvalue(threshold))...
%             , newline, 'slope= ', num2str(slope)...
%             , newline, 'D= ',num2str(slope/2), ' {\mu}m^2/sec'...
%             ],'FitBoxToText','on')
%         legend('off')
%         hold off
%         MY_DATA(k).slope = slope;
%         MY_DATA(k).r2 = MY_DATA(k).rsquared(threshold);
%         disp(k)
%         s3=subplot(2,2,[2,4]);
%         imagesc(MY_DATA(k).crop);
%         colormap(gray);
%         hold on
%         plot(s3,MY_DATA(k).particle_tracked(:,2),MY_DATA(k).particle_tracked(:,1),'-y');
%         title('overlay');
%         xlabel('x, px');
%         ylabel('y, px');
%         hold off
%         choice = menu('Menu','Next','Mark','End Session');
%         if choice == 1
%             close all
%         elseif choice == 2
%             prompt = 'Notes:';
%             MY_DATA(k).notes = input(prompt, 's');
%             close all
%         elseif choice == 3
%             close all
%             break
%         end
%     else
%                 figure
%         set(gcf,'Position',[60,30,1400,1100]);
%         subplot(2,2,1);
%         MSD = MY_DATA(k).MSD;
%         n=size(MY_DATA(k).particle_tracked(:,2), 1);
%         errorbar(MSD(1:n-1, 1), MSD(1:n-1, 2), MSD(1:n-1, 3), 'ok')
%         hold on
%         threshold = n/4;
%         x = MSD(1:threshold, 1);
%         y = MSD(1:threshold, 2);
%         mdl = fitlm(x,y);
%         plot(mdl);
%         slope = mdl.Coefficients.Estimate(2);
%         m = MY_DATA(k).name;
%         title(['NaN Best linear fit for ', m],'FontSize', 14)
%         xlabel('time (s)', 'FontSize', 14)
%         ylabel('MSD (um^2)', 'FontSize', 14)
%         legend('off')
%         save_name  = [m,'.png'];
%         hold off
% %         saveas(gcf, ['full length ', save_name])
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         subplot(2,2,3);
%         errorbar(MSD(1:threshold, 1), MSD(1:threshold, 2), MSD(1:threshold, 3), 'ok')
%         hold on
%         mdl = fitlm(x,y);
%         plot(mdl)
%         title(['NaN Best linear fit for ', m],'FontSize', 14)
%         xlabel('time (s)', 'FontSize', 14)
%         ylabel('MSD (um^2)', 'FontSize', 14)
%         annotation('textbox',[.13 .1 .1 .35],'String',['R^2= ', num2str(mdl.Rsquared.Ordinary)...
%             , newline, 'p-value= ', num2str(mdl.Coefficients.pValue(2))...
%             , newline, 'slope= ', num2str(slope)...
%             , newline, 'D= ',num2str(slope/2), ' {\mu}m^2/sec'...
%             ],'FitBoxToText','on')
%         legend('off')
%         hold off
%         MY_DATA(k).slope = slope;
%         MY_DATA(k).r2 = mdl.Rsquared.Ordinary;
%         disp(k)
%         s3=subplot(2,2,[2,4]);
%         imagesc(MY_DATA(k).crop);
%         colormap(gray);
%         hold on
%         plot(s3,MY_DATA(k).particle_tracked(:,2),MY_DATA(k).particle_tracked(:,1),'-y');
%         title('overlay');
%         xlabel('x, px');
%         ylabel('y, px');
%         hold off
%         choice = menu('Menu','Next','Mark','End Session');
%         if choice == 1
%             close all
%         elseif choice == 2
%             prompt = 'Notes:';
%             MY_DATA(k).notes = input(prompt, 's');
%             close all
%         elseif choice == 3
%             close all
%             break
%         end
%     end
% end


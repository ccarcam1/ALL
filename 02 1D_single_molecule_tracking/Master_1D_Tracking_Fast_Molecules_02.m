%% Master 1D Single Molecule Tracking Script
% Claudia Carcamo 
% 03 - 14 - 2020 

%% save the directory information here for the rest
% REQUIREMENTS
    % FUNCTION: my_directory_function('color',)
    % Must Start In Directory Of Interest
    % Directory must be made in jupyter notebook script
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\lambda nucleosome array\2021-01-25 dual color 700-1 array SWR1\';
   
cd(my_directory)
% my_directory_function('green');

%For two color analysis
my_directory_function('green red');
% my_directory_function('green');


%% Extract and save kymograph information
% REQUIREMENTS
    % FUNCTION: my_kymodata_structure --> no input arguments
    % Must Start In "Container" directory
    % Final data structure includes all colors 
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\lambda nucleosome array\2021-01-25 dual color 700-1 array SWR1\';
cd(my_directory)
cd container
my_kymodata_structure2('green red')
% my_kymodata_structure2('green')
load('data.mat')
length(data)
%% segment particles from different color channels
% REQUIREMENTS
    % FUNCTION: structure_name = my_segment_kymos_improvedwithzoom2(start_var, end_var, which_color)
    % Must Start In "Container" directory
    % linescan_time_mat
    % kymo_mat_green
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\lambda nucleosome array\2021-01-25 dual color 700-1 array SWR1\';
cd(my_directory)
cd([my_directory, 'container'])
% load('data.mat')
length(data)
green_segmentation = my_segment_kymos_improvedwithzoom2(11, 15, 'green'); % length(data)
%% Combine separate strcutures 
% a=0;% select this for the first one
a=length(green_segmentation_ALL);% for the first one, don't include this.
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
%% GREEN!~~~~~~~~ Fit gaussians to particles over time
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\lambda nucleosome array\2021-01-25 dual color 700-1 array SWR1\';
disp(strcat('Estimated Time to Completion:_', num2str(floor((2.6*length(structure_name))/60)), '_minutes and_',num2str(mod((2.5*length(structure_name)), 60)), '_seconds.'))
green_gaussfitting = my_gaussian_fitting_fast(1,length(structure_name),structure_name);
cd([my_directory, 'fitting']);
save('green_gaussfitting.mat', 'green_gaussfitting')
cd([my_directory, 'container']);
load('data.mat')
disp('Done with gaussian fitting')
warning('on','signal:findpeaks:largeMinPeakHeight')
%%%%%%%%%%%%%%%%%
cd([my_directory, 'MSD analysis']);
condensed_Struct = my_condense_relevant_info('green', data, green_gaussfitting, structure_name);
save('condensed_Struct_green.mat', 'condensed_Struct')
disp('Done making the condensed structure')
%% RED!~~~~~~~~~~ Fit gaussians to particles over time

% REQUIREMENTS
    % FUNCTION  gaussfitting = my_gaussian_fitting_fast(x,y,segmented_kymos)
    % max_distance_for_connect = 10 (hard coded inside this outer function)
    % x = starting point
    % y = ending (can be length(segmented_kymos))
    % segmented_kymos --> results from previous section 
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\lambda nucleosome array\2021-01-25 dual color 700-1 array SWR1\';

disp(strcat('Estimated Time to Completion:_', num2str(floor((2.6*length(structure_name))/60)), '_minutes and_',num2str(mod((2.5*length(structure_name)), 60)), '_seconds.'))
%%%%%%%%%%%%%%%%%%%%%%%%%
% green_gaussfitting = my_gaussian_fitting_fast(1,length(structure_name),structure_name);
red_gaussfitting = my_gaussian_fitting_fast(1,length(structure_name),structure_name);
%%%%%%%%%%%%%%%%%%%%%%%%%
cd([my_directory, 'fitting']);
%%%%%%%%%%%%%%%%%%%%%%%%%
save('red_gaussfitting.mat', 'red_gaussfitting')
% save('green_gaussfitting.mat', 'green_gaussfitting')
%%%%%%%%%%%%%%%%%%%%%%%%%

cd([my_directory, 'container']);
load('data.mat')
disp('Done with gaussian fitting')
warning('on','signal:findpeaks:largeMinPeakHeight')

% %% Determine Which Fits look good by eye (part 1 organize previous variables)
% % REQUIREMENTS
%     % FUNCTION keep_these_structure = my_visualize_fits(start_var, end_var, gaussfitting, data, structure_name)
cd([my_directory, 'MSD analysis']);
%%%%%%%%%%%%%%%%%%%%%%%%%
condensed_Struct = my_condense_relevant_info('red', data, red_gaussfitting, structure_name);
% condensed_Struct = my_condense_relevant_info('green', data, green_gaussfitting, structure_name);
%%%%%%%%%%%%%%%%%%%%%%%%%
% save('condensed_Struct_green.mat', 'condensed_Struct')
save('condensed_Struct_red.mat', 'condensed_Struct')
%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Done making the condensed structure')
% cd([my_directory, 'MSD analysis']);
% red_fitting_MSD_structure = my_condense_relevant_info('red', data, red_gaussfitting, structure_name);
% save('red_fitting_MSD_structure.mat', 'red_fitting_MSD_structure')


%% GREEN: Determine Which Fits look good by eye (part 2 visualize traces and decide which to keep)
% function keep_these_structure = my_visualize_fits(start_var, end_var, pre_MSD)
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\lambda nucleosome array\2021-01-25 dual color 700-1 array SWR1\';

cd([my_directory, 'MSD analysis']);
% load('green_fitting_MSD_structure.mat');
disp('..................')
startis = 1;
stopis = length(condensed_Struct); 
[keep_these, fit_again_later] = my_visualize_fits(startis, stopis, condensed_Struct);
save(['green_keep_these_good_first_round', num2str(startis), ' to ', num2str(stopis),'.mat'], 'keep_these')
save(['green_fit again later', num2str(startis), ' to ', num2str(stopis),'.mat'], 'fit_again_later')
close all
%% RED: Determine Which Fits look good by eye (part 2 visualize traces and decide which to keep)
% function keep_these_structure = my_visualize_fits(start_var, end_var, pre_MSD)
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\lambda nucleosome array\2021-01-25 dual color 700-1 array SWR1\';

cd([my_directory, 'MSD analysis']);
% load('green_fitting_MSD_structure.mat');
disp('..................')
startis = 1;
stopis = length(condensed_Struct); 
[keep_these, fit_again_later] = my_visualize_fits(startis, stopis, condensed_Struct);
save(['red_keep_these_good_first_round', num2str(startis), ' to ', num2str(stopis),'.mat'], 'keep_these')
save(['red_fit again later', num2str(startis), ' to ', num2str(stopis),'.mat'], 'fit_again_later')
close all
%% BOTH COLORS~~~~~~ COMBINE the information: Keep These*************
% a=0;% select this for the first one
a=length(green_keep_these_All);% for the first one, don't include this.
for i = 1:length(keep_these)
    keep_these_All(i+a) = keep_these(i);
end
disp("Done")

%% BOTH COLORS~~~~~~~ COMBINE the information: Fit Again Later
% fit_again_later_All = fit_again_later;
fit_again_later_All = horzcat(fit_again_later_All, fit_again_later);
disp("Done")


%% GREEN!~~~~~ Mask away regions creating errors and refit 
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\lambda nucleosome array\2021-01-25 dual color 700-1 array SWR1\';
disp('...............................')
masked_crops = my_masked_crops(condensed_Struct, fit_again_later, length(fit_again_later), length(fit_again_later));
save(['green_masked', num2str(length(fit_again_later)), ' to ', num2str(length(fit_again_later)),'.mat'], 'masked_crops')
disp('DONE!')
%% RED!~~~~~ Mask away regions creating errors and refit 
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\lambda nucleosome array\2021-01-25 dual color 700-1 array SWR1\';
disp('...............................')
masked_crops = my_masked_crops(condensed_Struct, fit_again_later, 1, length(fit_again_later));
save(['red_masked', num2str(length(fit_again_later)), ' to ', num2str(length(fit_again_later)),'.mat'], 'masked_crops')
disp('DONE!')
%% RED!~~~~~~~  SECOND GAUSS FITTING WITH MASKED DATA
clc
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\lambda nucleosome array\2021-01-25 dual color 700-1 array SWR1\';
red_gaussfitting_second_fitted_data = my_gaussian_fitting_fast_second_fit(structure_name, fit_again_later, masked_crops);
cd([my_directory, 'fitting']);
save('red_gaussfitting_second_fitted_data.mat', 'red_gaussfitting_second_fitted_data')

cd([my_directory, 'container']);
load('data.mat')
disp('Done with gaussian fitting')
warning('on','signal:findpeaks:largeMinPeakHeight')
%% GREEN!~~~~~~~  SECOND GAUSS FITTING WITH MASKED DATA
clc
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\lambda nucleosome array\2021-01-25 dual color 700-1 array SWR1\';
green_gaussfitting_second_fitted_data = my_gaussian_fitting_fast_second_fit(structure_name, fit_again_later, masked_crops);% structure_name are the segmented traces
cd([my_directory, 'fitting']);
save('green_gaussfitting_second_fitted_data.mat', 'green_gaussfitting_second_fitted_data')

cd([my_directory, 'container']);
load('data.mat')
disp('Done with gaussian fitting')
warning('on','signal:findpeaks:largeMinPeakHeight')
%%
for i = 1:length(fit_again_later)
    j = fit_again_later(i);
    segment_fit_later(i) = structure_name(j);% Insert segmentation ALL structure here
%     segment_fit_later(i).original_position = j;
end
%% RED!~~~~~~~~~ seconded condensed structure 
clc
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\lambda nucleosome array\2021-01-25 dual color 700-1 array SWR1\';

cd([my_directory, 'MSD analysis']);
condensed_Struct_second_fitted_data = my_condense_relevant_info2('red', data, red_gaussfitting_second_fitted_data, segment_fit_later);
save('red condensed_Struct_second_fitted_data.mat', 'condensed_Struct_second_fitted_data')
disp('Done making the condensed structure')
%% GREEN!~~~~~~~~~ seconded condensed structure 
clc
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\01 Data_Analysis\Projects\SWR1 Project\lambda nucleosome array\2021-01-25 dual color 700-1 array SWR1\';

cd([my_directory, 'MSD analysis']);
condensed_Struct_second_fitted_data = my_condense_relevant_info2('green', data, green_gaussfitting_second_fitted_data, segment_fit_later);
save('green condensed_Struct_second_fitted_data.mat', 'condensed_Struct_second_fitted_data')
disp('Done making the condensed structure')

%% RED!~~~~~~~~~~~~ Second Keep These PICKING
disp('..................')
startis =1;
stopis = length(condensed_Struct_second_fitted_data); 
%%%%%%%%%%%%%%%%ASSIGMNET%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make a new visualize fits 2 where fit again later is the original
%position in the array, so that it can be found and integrated easier. 
[keep_these, fit_again_later] = my_visualize_fits2(startis, stopis, condensed_Struct_second_fitted_data);
save(['red_keep_these_good_second_round', num2str(startis), ' to ', num2str(stopis),'.mat'], 'keep_these')
save(['red fit again later_secondround', num2str(startis), ' to ', num2str(stopis),'.mat'], 'fit_again_later')
close all

%% GREEN!~~~~~~~~~~~~ Second Keep These PICKING
disp('..................')
startis =1;
stopis = length(condensed_Struct_second_fitted_data); 
%%%%%%%%%%%%%%%%ASSIGMNET%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make a new visualize fits 2 where fit again later is the original
%position in the array, so that it can be found and integrated easier. 
[keep_these, fit_again_later] = my_visualize_fits2(startis, stopis, condensed_Struct_second_fitted_data);
save(['green_keep_these_good_second_round', num2str(startis), ' to ', num2str(stopis),'.mat'], 'keep_these')
save(['green fit again later_secondround', num2str(startis), ' to ', num2str(stopis),'.mat'], 'fit_again_later')
close all

%% Rescue poor fits:
% c = 1;
% for i = 1:length(keep_these)
%     if isempty(keep_these(i).line_time)
%         number(c) = i;
%         c = c+1;
%     end
% end

% check_these = fit_again_later_All(number);

% [keep_these, fit_again_later] = my_visualize_fits_test(condensed_Struct, check_these);
%% BOTH COLORS~~~~~~~ get rid of the blank cells get the MY_DATA structure
counter = 1;
for i = 1:length(keep_these)
    if ~isempty(keep_these(i).name)
        MY_DATA_2(counter)= keep_these(i);
        counter= counter+1;
    end
end
disp('DONE')
%% BOTH COLORS~~~~~~~ Combine the keep these data structures
% a=0;% select this for the first one
a=length(MY_DATA);% for the first one, don't include this.
for i = 1:length(MY_DATA_2)
    MY_DATA(i+a) = MY_DATA_2(i);
end
disp("Done")

%% RED~~~~~~~~~~~~~~ MSD calculation
% %function MSDs = my_MSD_calculator(structure_with_data)
% clc 
% clear
disp('New Run');
MSDs = my_MSD_calculator2(MY_DATA);
save('Red_MSDs.mat', 'MSDs')
% save('Green_MSDs.mat', 'MSDs')

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
save('Red_MY_DATA_MSDs.mat', 'MY_DATA')
% save('Green_MY_DATA_MSDs.mat', 'MY_DATA')
disp('Done')
%% GREEN~~~~~~~~~~~~~~ MSD calculation
% %function MSDs = my_MSD_calculator(structure_with_data)
% clc 
% clear
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
disp('Done')
%% BOTH COLORS~~~~~ Processing 1D diffusion 
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

% Process NaN values


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

%% RED!~~~~~~~~~~~~ Finalize Data Structures 
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
% Green_MY_DATA = MY_DATA;
Red_MY_DATA = MY_DATA;
save('Red_MY_DATA_final_structure.mat', 'Red_MY_DATA')
% save('Green_MY_DATA_final_structure.mat', 'Green_MY_DATA')
% save('Red_MY_DATA_final_structure.mat', 'MY_DATA')

%% GREEN!~~~~~~~~~~~~ Finalize Data Structures 
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

%% RED!~~~~~~~~~~~~ Get rid of poor fitting data
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

% save('Green_diffusioncoeff', 'diffusioncoeff');
% save('Green_Dr2less08', 'Dr2less08');
% save('Green_Dr2greater08','Dr2greater08');
% 
save('Red_diffusioncoeff', 'diffusioncoeff');
save('Red_Dr2less08', 'Dr2less08');
save('Red_Dr2greater08','Dr2greater08');

%% GREEN!~~~~~~~~~~~~ Get rid of poor fitting data
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
%%
figure(13)
hold on
% input_25 = slope_25(and(pval_25, r2_25))/2; 
% input_70 = slope_70(and(pval_70, r2_70))/2;
% input_150 = slope_150(and(pval_150, r2_150))/2;
input = diffusioncoeff;
what_is_it = 'SWR1 on array n = ';
% m = mean(slope(both)/2)
% input = line_time();
max = 0.14;
min = 0;
bins = 20;
type_of_plot = 'pdf';
type_style = 'bar';
% histogram(Dr2greater0825mM,20, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 2.125], 'DisplayName',strcat("Swc2 25mM KCl n = ", num2str(length(Dr2greater0825mM))))
% hold on
% histogram(diffusioncoeff25mM,20, 'Normalization', 'count','DisplayStyle','bar',"BinLimits",[0 2.125], 'DisplayName', strcat("Swc2 25mM KCl n = ", num2str(length(diffusioncoeff))))
histogram(input,bins, 'Normalization', type_of_plot,'DisplayStyle',type_style,"BinLimits",[min max], 'DisplayName',strcat(what_is_it, num2str(length(input))))
% histogram(diffusioncoeff,20, 'Normalization', 'count','DisplayStyle','bar',"BinLimits",[0 2.125], 'DisplayName', strcat("Swc2 70mM KCl n = ", num2str(length(diffusioncoeff))))
legend
% title("noATP with restrictions")
xlabel("diffusion cofficient um^2/sec")
ylabel("pdf")
% str = {[],[]};
% annotation('textbox', [0.65, 0.65, 0.1, 0.1], 'String', str)
title('')
%% Means and Sdev
input_25 = slope_25(and(pval_25, r2_25))/2; 
input_70 = slope_70(and(pval_70, r2_70))/2;
input_150 = slope_150(and(pval_150, r2_150))/2;
figure(8)
errorbar([25 70 150],[mean(input_25) mean(input_70) mean(input_150)],[std(input_25)/sqrt(length(input_25)) std(input_70)/sqrt(length(input_70)) std(input_150)/sqrt(length(input_150))],'o' )
axis([0 200 0 2])
yticks([0 0.5 1 1.5 2])
box off
% histogram(redDr2greater08)
% hold off
% dim = [0.7 0.6 0.3 0.3];
%% TEXT BOX FOR HISTOGRAMS
histogram(Dr2greater08, 20)
% hold on

% str = ['SWR1 n= ', num2str(length(greenDr2greater08))];%, newline, 'dCas9 n= ', num2str(length(redDr2greater08))];
% annotation('textbox',dim,'String',str,'FitBoxToText','on');




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
figure(1)
green_map = zeros(3);
for i = 1:11
green_map(i,2) = i*0.1-0.1;
end
j = 17;
str= 'kymo_16';
imagesc(flip(MY_DATA(j).kymograph))
colormap(green_map)
caxis([0 6])
hold on
for i = 1:length(MY_DATA)
    if contains(MY_DATA(i).name, str)
        green_chords = MY_DATA(i).crop_coordinates;
        green_position = MY_DATA(i).particle_tracked;
        plot(green_position(:,1)+green_chords(1)-1,green_position(:,2)+green_chords(3)-1,'Color','r','LineWidth',0.5)
    end
end
set(gca,'YTick',0:20:size(MY_DATA(j).kymograph, 1))
set(gca,'YTickLabel',0:20*0.1:size(MY_DATA(j).kymograph, 1)*0.1)
ylabel('Distance(um)')
set(gca,'XTick',0:100:size(MY_DATA(j).kymograph, 2))
set(gca,'XTickLabel',0:round(100*MY_DATA(j).line_time/1000):round(size(MY_DATA(j).kymograph, 2)*MY_DATA(j).line_time/1000))
xlabel('time(seconds)')
title('100nM AIM2-Cy3 diffusion')
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


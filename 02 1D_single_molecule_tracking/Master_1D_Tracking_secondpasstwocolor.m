%% Master 1D Single Molecule Tracking Script
% Claudia Carcamo 
% 03 - 14 - 2020 
%% save the directory information here for the rest
% REQUIREMENTS
    % FUNCTION: my_directory_function('color',)
    % Must Start In Directory Of Interest
    % Directory must be made in jupyter notebook script
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data_Analysis\Projects\SWR1 Project\sliding plus Cas9 marker\2020_09_02_dCas9_SWR1_twocolorsuccess\second pass two color\';
cd(my_directory)
% my_directory_function('green');

%For two color analysis
my_directory_function('green red');

%% Extract and save kymograph information
% REQUIREMENTS
    % FUNCTION: my_kymodata_structure --> no input arguments
    % Must Start In "Container" directory
    % Final data structure includes all colors 
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data_Analysis\Projects\SWR1 Project\sliding plus Cas9 marker\2020_09_02_dCas9_SWR1_twocolorsuccess\second pass two color\';
cd(my_directory)
cd container
my_kymodata_structure2('green red')
%% segment particles from different color channels
% REQUIREMENTS
    % FUNCTION: structure_name = my_segment_kymos(start_var, end_var, which_color)
    % Must Start In "Container" directory
    % linescan_time_mat
    % kymo_mat_green
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data_Analysis\Projects\SWR1 Project\sliding plus Cas9 marker\2020_09_02_dCas9_SWR1_twocolorsuccess\second pass two color\';
cd(my_directory)
cd([my_directory, 'container'])
load('data.mat')
% green_segmentation = my_segment_kymos_improvedwithzoom(1,3, 'green'); % length(data)
% end_here = length(uniq_kymo_names);
second_pass_segment = my_second_pass_two_color_particle_picking(Green_MY_DATA, Red_MY_DATA, 2, 16);
%% Combine separate strcutures 
% % a=0;
% a=length(green_segmentation);
% for i = 1:length(structure_name)
%     green_segmentation(i+a) = structure_name(i);
% end


% a=0;
a=length(red_segmentation);
for i = 1:length(structure_name)
    red_segmentation(i+a) = structure_name(i);
end

%% get rid of crops that are too small to correctly get gaussian fits
counter= 1;
for i = 1:length(green_segmentation)
    if size(green_segmentation(i).crop,1)<10
    recordthis(counter) = i ;
    counter = counter +1;
    end
end

%% Fit gaussians to particles over time
% REQUIREMENTS
    % FUNCTION gaussfitting = my_gaussian_fitting(x,y,segmented_kymos)
    % x = starting point
    % y = ending (can be length(segmented_kymos))
    % segmented_kymos --> results from previous section 
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data_Analysis\Projects\SWR1 Project\sliding plus Cas9 marker\2020_09_02_dCas9_SWR1_twocolorsuccess\second pass two color\';
cd(my_directory);
cd([my_directory, 'segmentation']);
load('two-color second pass segmentation 2 to 16.mat');
second_pass_gaussfitting = my_gaussian_fitting(1,length(second_pass_segment),second_pass_segment);
cd([my_directory, 'fitting']);
save('second_pass_gaussfitting.mat', 'second_pass_gaussfitting')

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
% my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\Projects\SWR1 Project\sliding + Cas9 marker\2020-03-04-Cy5_dCas9_Cy3_SWR1\';
% cd([my_directory, 'container']);
% load('data.mat')
% cd([my_directory, 'segmentation']);
% load('red_segmentation.mat')
% cd([my_directory, 'fitting']);
% load('red_gaussfitting.mat')
% cd([my_directory, 'MSD analysis']);
% red_fitting_MSD_structure = my_condense_relevant_info('red', data, red_gaussfitting, structure_name);
% save('red_fitting_MSD_structure.mat', 'red_fitting_MSD_structure')

clear
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data_Analysis\Projects\SWR1 Project\sliding plus Cas9 marker\2020_09_02_dCas9_SWR1_twocolorsuccess\second pass two color\';
cd([my_directory, 'container']);
load('data.mat')
cd([my_directory, 'segmentation']);
load('two-color second pass segmentation 2 to 16.mat')
cd([my_directory, 'fitting']);
load('second_pass_gaussfitting.mat')
cd([my_directory, 'MSD analysis']);
second_pass = my_condense_relevant_info('green', data, second_pass_gaussfitting, structure_name);
save('second_pass.mat', 'second_pass')

%% Determine Which Fits look good by eye (part 2 visualize traces and decide which to keep)
% function keep_these_structure = my_visualize_fits(start_var, end_var, pre_MSD)
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data_Analysis\Projects\SWR1 Project\sliding plus Cas9 marker\2020_09_02_dCas9_SWR1_twocolorsuccess\second pass two color\';
% cd([my_directory, 'MSD analysis']);
% load('second_pass.mat');
disp('..................')
startis = 1;
stopis= length(second_pass); 
second_pass_fitting_keep = my_visualize_fits(startis, stopis, second_pass);
save(['second_pass_fitting ', num2str(startis), ' to ', num2str(stopis),'.mat'], 'second_pass_fitting_keep')

% % Combine saved structures into one:
% a=length(green_keep_these_structure);
% for i = 1:length(keep_these_structure)
%     green_keep_these_structure(i+a) = keep_these_structure(i);
% end

% % my_visualize_fits2(1, 1, green_fitting_MSD_structure);
%% get rid of the blank cells get the MY_DATA structure
counter = 1;
for i = 1:length(second_pass_fitting_keep)
    if ~isempty(second_pass_fitting_keep(i).name)
        second_pass_MY_DATA(counter)= second_pass_fitting_keep(i);
        counter= counter+1;
    end
end
my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data_Analysis\Projects\SWR1 Project\sliding plus Cas9 marker\2020_09_02_dCas9_SWR1_twocolorsuccess\second pass two color\';
cd([my_directory, 'MSD analysis']);
save('SecondPass_MyData.mat', 'second_pass_MY_DATA')

%% MSD calculation
% %function MSDs = my_MSD_calculator(structure_with_data)
% clc 
% clear
% my_directory = 'C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data_Analysis\Projects\SWR1 Project\sliding plus Cas9 marker\2020_09_02_dCas9_SWR1_twocolorsuccess\second pass two color\';
% cd([my_directory, 'MSD analysis']);
% load('Red_MY_DATA.mat');

disp('New Run');
MSDs = my_MSD_calculator2(MY_DATA);
save('secondpass_MSDs.mat', 'MSDs')

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
save('secondpass_MY_DATA_MSDs.mat', 'MY_DATA')

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

save('secondpass_MY_DATA_final_structure.mat', 'MY_DATA')

%% Get rid of poor fitting data
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
%% TEXT BOX FOR HISTOGRAMS
histogram(Dr2greater08, 20)
hold on
histogram(diffusioncoeff)
hold off
dim = [0.7 0.6 0.3 0.3];
str = ['SWR1 n= ', num2str(length(Dr2greater08)), newline, 'dCas9 n= ', num2str(length(diffusioncoeff))];
annotation('textbox',dim,'String',str,'FitBoxToText','on');




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



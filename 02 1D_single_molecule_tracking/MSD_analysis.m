%% MSD calculation and linear fitting
clc;
clear;
% load('path_start');
path_start = 'D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\cas9\2019-11-06 cas9 3crRNAs ATTO550\';
% path_start = 'D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\2019-10-23 SWR1 gamma S\';
MSD_path = [path_start,'MSD analysis\'];
% load('MSD_path');
RESULTS = [path_start,'final\'];
% load('RESULTS');
cd(RESULTS);
load('keep_these');
cd (MSD_path);

x = 1;
% myexposuretime = 0.05135; % seconds
for i = 1:length(keep_these)
% for i = 2:3 % test size smaller faster to run
    if isempty(keep_these{1,i}) == 0
        L = length(keep_these{4,i});
        mydata(x).frames = [1:L]';
        mydata(x).coords = keep_these{2,i}(1:L,2);
        mydata(x).intensity = keep_these{2,i}(1:L,3); 
        mydata(x).mol_id = keep_these{1,i};
        mydata(x).timestep = keep_these{5,i}; 
        mydata(x).full_kymo = keep_these{6,i};
        x = x+1; 
    else
    end    
end
for i = 1:length(mydata)
    mydata(i).molecule_num = i;
end
% next part
% myexposuretime = 0.05135; % seconds
frames = mydata.frames; % frames from one trajectory
coords = mydata.coords; % coords from one trajectory
% tau = myexposuretime;
pixSize = 0.100; % micrometers

%% comprehensive figure
x = 1;
for i = 1:length(mydata)
    for j = 1:length(keep_these)
        if strcmp(mydata(i).mol_id, keep_these{1,j})
            mydata(x).crop_coords = keep_these{3,j};
            mydata(x).crop = keep_these{4,j};
            x = x+1;
        else
        end
    end
end
close all

%%
x= 1;
for k = 1:length(mydata) % test size smaller faster to run
n=size(mydata(k).coords, 1);
r_all=[];
for i = 1:n-1 % i is the frame lag
    r=zeros(n-i,2);
    for j = 1: n-i % n-i is the total number of lags one can calculate for frame lag t
        r(j, 1)=mydata(k).frames(j+i)-mydata(k).frames(j); % convert frame lag to time lag
        dif = (pixSize) *(mydata(k).coords(j+i, :)-mydata(k).coords(j, :));%pixSize/1000 was removed
        r(j, 2)=sum(dif.^2); % sqaure displacement
    end
    r_all(end+1:end+size(r, 1), :)=r;
end
mydata(x).r_all = r_all; % Not yet ready to implement
        x = x+1; 
disp(k); 
end



% MSD_analyzer - the once master function
% FIRST! calculate the MSD
for k = 1:length(mydata)
    tau = mydata(k).timestep/1000;
    d = mydata(k).r_all;
% set up some empty matrices...
    d_all = [];


% calculates ALL displacements
% % for i = 1:size(mydata(k).frames, 1)
% %     if length(mydata(k).frames) > 1
% % %         d = TrajDispl(coords(:, 1), frames, pixSize);
% %         d_all(end+1:end+size(d, 1), :) = d;
% %     end
% % end

    time_lags = unique(d(:, 1));
    n=size(mydata(k).coords, 1);
    D_all = zeros(n-1,5);
% finds the mean displacement squared from above displacements
    for i = 1:n-1
        ind=find(d(:, 1) == time_lags(i));
        D_all(i, 1) = tau * time_lags(i); % converts to real time lag
        D_all(i, 2) = mean(d(ind, 2)); %MSD in x-direction for timelag i
% %         D_all(i, 3) = std(d(ind, 2))/sqrt(length(ind)); %sem in x-direction for timelag i***********************************D13 
% %         D_all(i, 4) = std(d(ind, 2)); % std
        D_all(i, 5) = length(d(ind, 2)); % num points
    end
    mydata(k).MSD = D_all;
    disp(k)
% MSD = D_all;
% out(k).MSD = MSD;
end


% calculate fits
for k = 1: length(mydata)
MSD = mydata(k).MSD;
record_this = zeros(size(mydata(k).coords, 1)-1);
record_that = zeros(size(mydata(k).coords, 1)-1);
    for i = 3:size(mydata(k).coords, 1)-1
    xx = MSD(1:i, 1);
    yy = MSD(1:i, 2);
    mdl = fitlm(xx,yy);
    record_this(i) = mdl.Rsquared.Ordinary;
    record_that(i) = mdl.Coefficients.pValue(2);
    end
disp(k)
mydata(k).rsquared = record_this;
mydata(k).pvalue = record_that;
end
%%
% Find a cut off point for fitting a line 

for k = 1:length(mydata)
%     x = 1;
    r2_list = (mydata(k).rsquared);
    for i = 2:length(r2_list)-5
        a = r2_list(i+1)-r2_list(i);
        aa = r2_list(i)-r2_list(i-1);
        if a<0 && mydata(k).pvalue(i)<0.01 && aa>0
            b = r2_list(i+2)-r2_list(i+1);
            c = r2_list(i+3)-r2_list(i+2);
            d = r2_list(i+4)-r2_list(i+3);
            e = r2_list(i+5)-r2_list(i+4);
            if b<0 && c<0 && d<0 && e<0
                mydata(k).fit_cutoff2 = i-1;
                break
            else
            end
        else
        end
        mydata(k).fit_cutoff2 = NaN;
    end
end

%% Length relationship to cutoff point
fit_cutoff2 = [];
x= 1;
for i = 1:length(mydata)
%     if ~isnan(mydata(i).fit_cutoff2)
        fit_cutoff2(1,x)= mydata(i).fit_cutoff2;
        fit_cutoff2(2,x)= length(mydata(i).MSD);
        fit_cutoff2(3,x) = mydata(i).fit_cutoff2*(mydata(i).timestep/1000);
        x = x+1;
%     else
%     end
end
%
figure
histogram(fit_cutoff2(3,:),50,'BinLimits', [0 40])
% title('histogram of trace lengths')
xlabel('time(s)')
ylabel('counts')
timepertrace = [];
% mydata = mydata_allATP;
for i = 1:length(mydata)
    timepertrace(i) = length(mydata(i).frames)*(mydata(i).timestep/1000);
end
hold on
histogram(timepertrace,50,'BinLimits', [0 40])'

%%
scatter(fit_cutoff2(3,:),timepertrace)
%% just extract values
close all
for k = 1: length(mydata)
    if ~isnan(mydata(k).fit_cutoff2)
        MSD = mydata(k).MSD;
        threshold = mydata(k).fit_cutoff2;
        x = MSD(1:threshold, 1);
        y = MSD(1:threshold, 2);
        mdl = fitlm(x,y);
        plot(mdl);
        mydata(k).slope = mdl.Coefficients.Estimate(2);
        mydata(k).r2 = mydata(k).rsquared(threshold);
        mydata(k).intercept = mdl.Coefficients.Estimate(1);
        disp(k)
        close all
    else
        MSD = mydata(k).MSD;
        threshold = length(MSD)/4;
        x = MSD(1:threshold, 1);
        y = MSD(1:threshold, 2);
        mdl = fitlm(x,y);
        plot(mdl);
        mydata(k).slope = mdl.Coefficients.Estimate(2);
        mydata(k).r2 = mdl.Rsquared.Ordinary;
        mydata(k).intercept = mdl.Coefficients.Estimate(1);
        disp(k)
        close all
    end
end


%% plot with correct fits 0.01 pvalue **** Careful*****
close all
for k = 1: length(mydata)
    if ~isnan(mydata(k).fit_cutoff2)
        figure(1)
        MSD = mydata(k).MSD;
        n=size(mydata(k).coords, 1);
        errorbar(MSD(1:n-1, 1), MSD(1:n-1, 2), MSD(1:n-1, 3), 'ok')
        hold on
        threshold = mydata(k).fit_cutoff2;
        x = MSD(1:threshold, 1);
        y = MSD(1:threshold, 2);
        mdl = fitlm(x,y);
        plot(mdl);
        slope = mdl.Coefficients.Estimate(2);
        m = mydata(k).mol_id;
        title(['Best linear fit for ', m],'FontSize', 14)
        xlabel('time (s)', 'FontSize', 14)
        ylabel('MSD (um^2)', 'FontSize', 14)
        legend('off')
        save_name  = [m, ' molecule number ', num2str(mydata(k).molecule_num),'.png'];
        hold off
        saveas(gcf, ['full length ', save_name])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(2)
        errorbar(MSD(1:threshold, 1), MSD(1:threshold, 2), MSD(1:threshold, 3), 'ok')
        hold on
        mdl = fitlm(x,y);
        plot(mdl)
        title(['Best linear fit for ', m],'FontSize', 14)
        xlabel('time (s)', 'FontSize', 14)
        ylabel('MSD (um^2)', 'FontSize', 14)
        annotation('textbox',[.14 .3 .4 .6],'String',['R^2= ', num2str(mydata(k).rsquared(threshold))...
            , newline, 'p-value= ', num2str(mydata(k).pvalue(threshold))...
            , newline, 'slope= ', num2str(slope)...
            , newline, 'D= ',num2str(slope/2), ' \mum^2/sec'...
            ],'FitBoxToText','on')
        hold off
        saveas(gcf, ['fit range ', save_name])
        mydata(k).slope = slope;
        mydata(k).r2 = mydata(k).rsquared(threshold);
        disp(k)
        close all
    else
    end
end


%% Screen Shot Interesting Ones, Write Notes Next to Interesting Ones
% Save the structure indicatting that the version has notes 
for k = 1: length(mydata)
% for i = 1:length(fast)
%     k = fast(i);
    if ~isnan(mydata(k).fit_cutoff2)
        figure
        set(gcf,'Position',[60,30,1400,1100]);
        subplot(2,2,1);
        MSD = mydata(k).MSD;
        n=size(mydata(k).coords, 1);
        errorbar(MSD(1:n-1, 1), MSD(1:n-1, 2), MSD(1:n-1, 3), 'ok')
        hold on
        threshold = mydata(k).fit_cutoff2;
        x = MSD(1:threshold, 1);
        y = MSD(1:threshold, 2);
        mdl = fitlm(x,y);
        plot(mdl);
        slope = mdl.Coefficients.Estimate(2);
        match = ["1106", "num"];
        m = erase(mydata(k).mol_id,match);
        title(['', m],'FontSize', 14)        
%         title(['Best linear fit for ', m],'FontSize', 14)
        xlabel('time (s)', 'FontSize', 14)
        ylabel('MSD (um^2)', 'FontSize', 14)
        legend('off')
        save_name  = [m, ' molecule number ', num2str(mydata(k).molecule_num),'.png'];
        hold off
%         saveas(gcf, ['full length ', save_name])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(2,2,3);
        errorbar(MSD(1:threshold, 1), MSD(1:threshold, 2), MSD(1:threshold, 3), 'ok')
        hold on
        mdl = fitlm(x,y);
        plot(mdl)
        title(['Best linear fit for ', m],'FontSize', 14)
        xlabel('time (s)', 'FontSize', 14)
        ylabel('MSD (um^2)', 'FontSize', 14)
        annotation('textbox',[.13 .1 .1 .35],'String',['R^2= ', num2str(mydata(k).rsquared(threshold))...
            , newline, 'p-value= ', num2str(mydata(k).pvalue(threshold))...
            , newline, 'slope= ', num2str(slope)...
            , newline, 'D= ',num2str(slope/2), ' {\mu}m^2/sec'...
            ],'FitBoxToText','on')
        legend('off')
        hold off
        mydata(k).slope = slope;
        mydata(k).r2 = mydata(k).rsquared(threshold);
        disp(k)
        s3=subplot(2,2,[2,4]);
        imagesc(mydata(k).crop);
        colormap(gray);
        hold on
        plot(s3,mydata(k).coords,mydata(k).frames,'-y');
        title('overlay');
        xlabel('x, px');
        ylabel('y, px');
        hold off
        choice = menu('Menu','Next','Mark','End Session');
        if choice == 1
            close all
        elseif choice == 2
            prompt = 'Notes:';
            mydata(k).notes = input(prompt, 's');
            close all
        elseif choice == 3
            close all
            break
        end
    else
                figure
        set(gcf,'Position',[60,30,1400,1100]);
        subplot(2,2,1);
        MSD = mydata(k).MSD;
        n=size(mydata(k).coords, 1);
        errorbar(MSD(1:n-1, 1), MSD(1:n-1, 2), MSD(1:n-1, 3), 'ok')
        hold on
        threshold = n/4;
        x = MSD(1:threshold, 1);
        y = MSD(1:threshold, 2);
        mdl = fitlm(x,y);
        plot(mdl);
        slope = mdl.Coefficients.Estimate(2);
        m = mydata(k).mol_id;
        title(['NaN Best linear fit for ', m],'FontSize', 14)
        xlabel('time (s)', 'FontSize', 14)
        ylabel('MSD (um^2)', 'FontSize', 14)
        legend('off')
        save_name  = [m, ' molecule number ', num2str(mydata(k).molecule_num),'.png'];
        hold off
%         saveas(gcf, ['full length ', save_name])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(2,2,3);
        errorbar(MSD(1:threshold, 1), MSD(1:threshold, 2), MSD(1:threshold, 3), 'ok')
        hold on
        mdl = fitlm(x,y);
        plot(mdl)
        title(['NaN Best linear fit for ', m],'FontSize', 14)
        xlabel('time (s)', 'FontSize', 14)
        ylabel('MSD (um^2)', 'FontSize', 14)
        annotation('textbox',[.13 .1 .1 .35],'String',['R^2= ', num2str(mdl.Rsquared.Ordinary)...
            , newline, 'p-value= ', num2str(mdl.Coefficients.pValue(2))...
            , newline, 'slope= ', num2str(slope)...
            , newline, 'D= ',num2str(slope/2), ' {\mu}m^2/sec'...
            ],'FitBoxToText','on')
        legend('off')
        hold off
        mydata(k).slope = slope;
        mydata(k).r2 = mdl.Rsquared.Ordinary;
        disp(k)
        s3=subplot(2,2,[2,4]);
        imagesc(mydata(k).crop);
        colormap(gray);
        hold on
        plot(s3,mydata(k).coords,mydata(k).frames,'-y');
        title('overlay');
        xlabel('x, px');
        ylabel('y, px');
        hold off
        choice = menu('Menu','Next','Mark','End Session');
        if choice == 1
            close all
        elseif choice == 2
            prompt = 'Notes:';
            mydata(k).notes = input(prompt, 's');
            close all
        elseif choice == 3
            close all
            break
        end
    end
end

%% pvalue and rsquared value changes
figure
k =126;
% set(gcf,'Position',[60,30,1400,1100]);
% subplot(2,2,1);


MSD = mydata(k).MSD;
n=size(mydata(k).coords, 1);
errorbar(MSD(1:n-1, 1), MSD(1:n-1, 2), MSD(1:n-1, 3), 'ok')
hold on
threshold = mydata(k).fit_cutoff2;
x = MSD(1:threshold, 1);
y = MSD(1:threshold, 2);
mdl = fitlm(x,y);
plot(mdl);
slope = mdl.Coefficients.Estimate(2);
m = mydata(k).mol_id;
% title(['Best linear fit for ', m],'FontSize', 14)
xlabel('time (s)', 'FontSize', 14)
ylabel('MSD (um^2)', 'FontSize', 14)
legend('off')
save_name  = [m, ' molecule number ', num2str(mydata(k).molecule_num),'.png'];
% hold on
% yyaxis left
% ylabel('R^2', 'FontSize', 14, 'Color', [0    0.4470    0.7410])
% plot(MSD(1:n-1, 1), mydata(k).rsquared(:,1),'-','LineWidth',2)
% 
% yyaxis right
% 
% 
% plot(MSD(1:n-1, 1), mydata(k).pvalue(:,1),'-','LineWidth',2,  'Color', 'red')
% ylabel('p-value', 'FontSize', 14, 'Color', 'red')
% hold off

%%
slope = [];
diffusioncoeff = [];
intercept = [];
r2 = [];
x = 1;
for i = 1:length(mydata)
%         if mydata(i).r2 > 0.8
        slope(x) = mydata(i).slope;
        r2(x) = mydata(i).r2;
        intercept(x) = mydata(i).intercept;
        x = x+1;
%         else
%         end
end
diffusioncoeff = slope/2;

%%
figure(1)
histogram(diffusioncoeff, 30)
title('diffusion coefficients with ATP');
ylabel('counts');
xlabel('D(\mum^2/sec)');
% saveas(gcf,'diffusion coefficients with ATP.png')
figure(2)
histogram(r2,30)
title('R^2 values with ATP');
ylabel('counts');
xlabel('R^2');
% saveas(gcf,'R^2 values with ATP.png');
%%
r2_25mMKCl = r2;
diffusioncoeff_25mMKCl = diffusioncoeff;
mydata_25mMKCl= mydata;
intercept_25mMKCl = intercept;

save('r2_25mMKCl', 'r2_25mMKCl')
save('diffusioncoeff_25mMKCl', 'diffusioncoeff_25mMKCl')
save('mydata_25mMKCl', 'mydata_25mMKCl')
save('intercept_25mMKCl', 'intercept_25mMKCl')


%%%
%%%%
%%%%%%
%%%%%%%%
%%%%%
%%%%
%%%


%% compare conditions
condition1 = diffusioncoeff_allnoATP;
description1 = 'without ATP';
color1 = 'red';
condition2 = diffusioncoeff_allATP;
description2 = 'with ATP';
color2 = 'blue';
% condition3 = diffusioncoeff_gammaS;
% description3 = 'with ATPgS';
% color3 = 'yellow';
% condition4 = diffusioncoeff_cas9;
% description4 = 'cas9';
% color4 = 'black';

figure
binlims = [-.01 0.1];
histogram(condition1, 30,'FaceColor', color1, 'BinLimits', binlims)
hold on
histogram(condition2, 30,'FaceColor', color2,'BinLimits', binlims)
% histogram(condition3, 30,'FaceColor', color3,'BinLimits', binlims)
% histogram(condition4, 30,'FaceColor', color4,'BinLimits', binlims)
hold off
annotation('textbox',[.63 .3 .4 .6],'String',[description1,' in ', color1...
    , newline, description2,' in ', color2...
...     , newline, description3,' in ', color3...
...     , newline, description4,' in ', color4...
    ],'FitBoxToText','on')
title(['diffusion coefficients ', description1, ' and ', description2...
%     , newline...
...     , ' and ',description3,...
...     ' and ',description4...
    ]);
ylabel('counts');
xlabel('D(um^2/sec)');
% saveas(gcf,['diffusion coefficients', description1, ' and ',description2...
% %     ,' and ',description3...
% %     , ' and ',description4,...
%     ,'.png']);

%% get just R2 and Slope
close all
for k = 1: length(mydata)
MSD = mydata(k).MSD;
n=size(mydata(k).coords, 1);
threshold = mydata(k).fit_cutoff2;
x = MSD(1:threshold, 1);
y = MSD(1:threshold, 2);
% [p,S] = polyfit(x,y,1);
% [y_fit, delta] = polyval(p,x,S);
% plot(x,y_fit,'-')
mdl = fitlm(x,y);
slope = mdl.Coefficients.Estimate(2);
mydata(k).slope = slope;
mydata(k).r2 = mydata(k).rsquared(threshold);
disp(k)
end
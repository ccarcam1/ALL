% %% combine gamma S results 
% MY_DATA= MY_DATA_gammaS;
% 
% %%
% 
% x = length(MY_DATA) +1;
% for i = 1:length(MY_DATA_gammaS)
%     MY_DATA(x) = MY_DATA_gammaS(i)
%     x = x+1;
% end

%% processing lengths of traces 

% I will consider anything shorter than 0.75 seconds to be short 

for i = 1:length(MY_DATA)
    holdingplace(i) = size(MY_DATA(i).particle_tracked(:,2), 1)*(MY_DATA(i).line_time/1000);
end
% 
lookhere = find(holdingplace < 2)
x = 1;
y = 1;
for i = 1:length(MY_DATA)
    if holdingplace(i) >2
       tempMY_DATA(x) = MY_DATA(i);
       x = x +1;
    else
        shortMY_DATA(y) = MY_DATA(i);
       y = y +1;
    end
end
% 
clear MY_DATA
MY_DATA = tempMY_DATA;
%% Process NaN values
% MY_DATA = MY_DATA_cas9;

NaN_negative_which = [];
NaN_which = [];
NaN_positive_which = [];
negative_which_noNaN = [];

x = 1;
y = 1;
z = 1;
a = 1;
for k = 1: length(MY_DATA)
    if isnan(MY_DATA(k).fit_cutoff2)
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
%% Validate NaN values are what I think they are
% Save the structure indicatting that the version has notes 
% check these
% % NaN_which
% % negative_which_noNaN

for i = 1:length(negative_which_noNaN)
%     k = lookhere(i);
% for i = 1:length(NaN_negative_which)
    k = NaN_negative_which(i);
    if ~isnan(MY_DATA(k).fit_cutoff2)
        figure
        set(gcf,'Position',[60,30,1400,1100]);
        subplot(2,2,1);
        MSD = MY_DATA(k).MSD;
        n=size(MY_DATA(k).coords, 1);
        errorbar(MSD(1:n-1, 1), MSD(1:n-1, 2), MSD(1:n-1, 3), 'ok')
        hold on
        threshold = MY_DATA(k).fit_cutoff2;
        x = MSD(1:threshold, 1);
        y = MSD(1:threshold, 2);
        mdl = fitlm(x,y);
        plot(mdl);
        slope = mdl.Coefficients.Estimate(2);
        m = MY_DATA(k).mol_id;
        title(['Best linear fit for ', m],'FontSize', 14)
        xlabel('time (s)', 'FontSize', 14)
        ylabel('MSD (um^2)', 'FontSize', 14)
        legend('off')
        save_name  = [m, ' molecule number ', num2str(MY_DATA(k).molecule_num),'.png'];
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
        annotation('textbox',[.13 .1 .1 .35],'String',['R^2= ', num2str(MY_DATA(k).rsquared(threshold))...
            , newline, 'p-value= ', num2str(MY_DATA(k).pvalue(threshold))...
            , newline, 'slope= ', num2str(slope)...
            , newline, 'D= ',num2str(slope/2), ' {\mu}m^2/sec'...
            ],'FitBoxToText','on')
        legend('off')
        hold off
        MY_DATA(k).slope = slope;
        MY_DATA(k).r2 = MY_DATA(k).rsquared(threshold);
        disp(k)
        s3=subplot(2,2,[2,4]);
        imagesc(MY_DATA(k).crop);
        colormap(gray);
        hold on
        plot(s3,MY_DATA(k).coords,MY_DATA(k).particle_tracked(:,2),'-y');
        title('overlay');
        xlabel('x, px');
        ylabel('y, px');
        hold off
        choice = menu('Menu','Next','Mark','End Session');
        if choice == 1
            close all
        elseif choice == 2
            prompt = 'Notes:';
            MY_DATA(k).notes = input(prompt, 's');
            close all
        elseif choice == 3
            close all
            break
        end
    else
                figure
        set(gcf,'Position',[60,30,1400,1100]);
        subplot(2,2,1);
        MSD = MY_DATA(k).MSD;
        n=size(MY_DATA(k).coords, 1);
        errorbar(MSD(1:n-1, 1), MSD(1:n-1, 2), MSD(1:n-1, 3), 'ok')
        hold on
        threshold = n/4;
        x = MSD(1:threshold, 1);
        y = MSD(1:threshold, 2);
        mdl = fitlm(x,y);
        plot(mdl);
        slope = mdl.Coefficients.Estimate(2);
        m = MY_DATA(k).mol_id;
        title(['NaN Best linear fit for ', m],'FontSize', 14)
        xlabel('time (s)', 'FontSize', 14)
        ylabel('MSD (um^2)', 'FontSize', 14)
        legend('off')
        save_name  = [m, ' molecule number ', num2str(MY_DATA(k).molecule_num),'.png'];
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
        MY_DATA(k).slope = slope;
        MY_DATA(k).r2 = mdl.Rsquared.Ordinary;
        disp(k)
        s3=subplot(2,2,[2,4]);
        imagesc(MY_DATA(k).crop);
        colormap(gray);
        hold on
        plot(s3,MY_DATA(k).coords,MY_DATA(k).particle_tracked(:,2),'-y');
        title('overlay');
        xlabel('x, px');
        ylabel('y, px');
        hold off
        choice = menu('Menu','Next','Mark','End Session');
        if choice == 1
            close all
        elseif choice == 2
            prompt = 'Notes:';
            MY_DATA(k).notes = input(prompt, 's');
            close all
        elseif choice == 3
            close all
            break
        end
    end
end

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

%%
slope = [];
diffusioncoeff = [];
intercept = [];
r2 = [];
x = 1;
for i = 1:length(MY_DATA)
        slope(i) = MY_DATA(i).slope;
        r2(i) = MY_DATA(i).r2;
        intercept(i) = MY_DATA(i).intercept;
        if MY_DATA(i).r2>0.8
            MY_DATA(i).r2point8 = 1;
            x = x+1;
        else 
            MY_DATA(i).r2point8 = 0;
        end
end
diffusioncoeff= slope/2;
Dr2less08 = diffusioncoeff(r2<0.8);
Dr2greater08 = diffusioncoeff(r2>0.8);

save('Dr2less08', 'Dr2less08');
save('Dr2greater08','Dr2greater08');
%%
r2_gammaS= r2;
diffusioncoeff_gammaS = diffusioncoeff;
MY_DATA_gammaS= MY_DATA;
intercept_gammaS = intercept;
Dr2less08_gammaS = Dr2less08;
Dr2greater08_gammaS = Dr2greater08;

save('r2_gammaS', 'r2_gammaS')
save('diffusioncoeff_gammaS', 'diffusioncoeff_gammaS')
save('MY_DATA_gammaS', 'MY_DATA_gammaS')
save('intercept_gammaS', 'intercept_gammaS')
save('Dr2less08_gammaS','Dr2less08_gammaS'); 
save('Dr2greater08_gammaS','Dr2greater08_gammaS');

%% load information on all the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cas9
%nonspecific
% cd('D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\cas9\2019-10-30 summary\MSD analysis\noNaNs\');
% specific
cd('D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\cas9\2019-11-06 cas9 3crRNAs ATTO550\MSD analysis\');
% noATP
cd('D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\noATP\pool no ATP data\MSD analysis\No Nan\');
% ATP
cd('D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\plus ATP\pool of ATP data\MSD analysis\without NaNs\');
cd('C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\salt concentrations ATP\plus ATP 70 mM KCl\pool of ATP data\MSD analysis\without NaNs\');
% gammaS
cd('D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\ATP gammaS\');
% 200mM KCl SWR1
cd('D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\2019-11-20 200 mM KCl\MSD analysis\noNaN');
cd('C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\2019-11-25 25 mM KCl\MSD analysis\no NaN');

% 25 mM KCl SWR1
cd('D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\2019-11-25 25 mM KCl\MSD analysis\no NaN');
cd('C:\Users\carca\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\2019-11-25 25 mM KCl\MSD analysis\no NaN');
%%
clf
figure(1)
binlims = [0 0.15];
% histogram(diffusioncoeff_allATP,20,'BinLimits', binlims, 'FaceColor', '#A2142F','EdgeColor', '#A2142F','FaceAlpha',0.8, 'LineWidth', 3)
histogram(Dr2greater08_allATP, 30, 'BinLimits', binlims,'DisplayStyle', 'stairs','EdgeColor', '#A2142F','EdgeAlpha',1, 'LineWidth', 3)
hold on
% histogram(diffusioncoeff_allnoATP,20,'BinLimits', binlims, 'FaceColor', '#0072BD','EdgeColor', '#0072BD','FaceAlpha',0.6, 'LineWidth', 3)
histogram(Dr2greater08_200mMKCl_noNaN, 30, 'BinLimits', binlims, 'DisplayStyle', 'stairs','EdgeColor', '#4DBEEE','EdgeAlpha',1, 'LineWidth', 3)
% hold off
% histogram(diffusioncoeff_allgammaS,20, 'BinLimits', binlims,'FaceColor', '#77AC30','EdgeColor','#77AC30','FaceAlpha',0.3, 'LineWidth', 3)
histogram(Dr2greater08_25mMKCl_noNaN, 30, 'BinLimits', binlims, 'DisplayStyle', 'stairs','EdgeColor', '#77AC30','EdgeAlpha',1, 'LineWidth', 3)

hold on
% histogram(diffusioncoeff_Cas9_5pN,20, 'BinLimits', binlims, 'FaceColor', 'black','EdgeColor', 'black','FaceAlpha',0.5, 'LineWidth', 3)
% histogram(diffusioncoeff_Cas9_5pN, 30, 'BinLimits', binlims, 'DisplayStyle', 'stairs','EdgeColor', 'black','EdgeAlpha',1, 'LineWidth', 3)

hold off
% title('diffusion coefficients specific cas9');
ylabel('counts');
xlabel('D(\mum^2/sec)');
% saveas(gcf,'stairs specifc cas9 diffusion coefficients with ATP bin limits point2.png');

%%
length(diffusioncoeff_allATP)
length(diffusioncoeff_allnoATP)
length(diffusioncoeff_allgammaS)
length(diffusioncoeff_allCas9)


%%
MY_DATA = MY_DATA_allATP;
for i = 1:length(MY_DATA)
    timepertrace(i) = length(MY_DATA(i).particle_tracked(:,2))*(MY_DATA(i).line_time/1000);
end

%%
histogram(timepertrace, 30)


%%
figure
x = [];
histogram(diffusioncoeff_allCas9, 30,  'BinLimits', ...
    [0, max(diffusioncoeff_allCas9)],'FaceColor','#7E2F8E','FaceAlpha',1)
x = [0:0.001:.045];

pd = fitdist(diffusioncoeff_allCas9','Gamma');
y = pdf(pd,x);
hold on
plot(x,y,'Color', 'm','LineWidth', 1)
% hold off

% figure
x = [];
histogram(diffusioncoeff_allnoATP, 30,  'BinLimits', ...
    [0, max(diffusioncoeff_allnoATP)], 'FaceColor','#0072BD','FaceAlpha',0.75 )
x = [0:0.001:.12];

pd = fitdist(diffusioncoeff_allnoATP','Gamma');
y = pdf(pd,x);
hold on
plot(x,y,'Color', 'b','LineWidth', 1)
% hold off

% figure
x = [];
histogram(diffusioncoeff_allATP, 30,  'BinLimits', ...
    [0, max(diffusioncoeff_allATP)],'FaceColor','#A2142F','FaceAlpha',0.5)
x = [0:0.001:.12];

pd = fitdist(diffusioncoeff_allATP','Gamma');
y = pdf(pd,x);
hold on
plot(x,y, 'Color', 'r','LineWidth', 1)
% hold off

% figure
x = [];
histogram(diffusioncoeff_allgammaS, 30,  'BinLimits', ...
    [0, max(diffusioncoeff_allgammaS)],'FaceColor', '#EDB120','FaceAlpha',0.3)
x = [0:0.001:.12];

pd = fitdist(diffusioncoeff_allgammaS','Gamma');
y = pdf(pd,x);
hold on
plot(x,y,'Color', 'y','LineWidth', 1)
hold off
%%
data = diffusioncoeff_allCas9*1000;
[N,~] = histcounts(data,20);
histogram(data,20)
hold on
pd = fitdist(round(data'),'Poisson');
x = [0: max(round(data'))]
y = pdf(pd,x);
plot(x,y/max(y)*max(N))

%%
pd = fitdist(diffusioncoeff_allCas9','Poisson');
clf
histogram(diffusioncoeff_allCas9, 10,  'BinLimits', ...
    [0, max(diffusioncoeff_allCas9)])
[N,edges] = histcounts(diffusioncoeff_allCas9, 10);
hold on
x = [];
for i = 1:length([0:0.1:max(N)])
x(i) = pd.lambda;
end
plot(x, [0:0.1:max(N)])
hold off

%%
data = diffusioncoeff_allgammaS*1000;
[N,~] = histcounts(data,20);
histogram(data,20)
hold on
pd = fitdist(round(data'),'Normal');
x = [0: max(round(data'))]
y = pdf(pd,x);
plot(x,y/max(y)*max(N))


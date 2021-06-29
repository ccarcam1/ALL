%% Indices
% experimental conditions: 
% 70mM KCl + 1mM ATP,                   ATP_70mMKCl      
% 25mM KCl + 1mM ATP,                   ATP_25mMKCl
% 200mM KCl + 1mM ATP,                  ATP_200mMKCl
% gammaS,                               gammaS
% dCas9,                                dCas9
% ADP,                                  ADP
% noATP,                                noATP
% "SWR1 70mM KCl 1mM ATP dual color",   ATP_70mMKCl_dualcolor
% "cas9 in dual color"                  dCas9_twocolor
% "SWR1 + ZB + noATP"                   ZB_noATP

% experiment = {FINAL_DATA.experiment}.';
% particle_tracked = {FINAL_DATA.particle_tracked}.';
% line_time = [FINAL_DATA.line_time].';

whichone = "200mM KCl + 1mM ATP";% CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_holding = [];

for i = 1:length(experiment)
    if matches(experiment{i}, whichone)
        var_holding(i) = 1; 
    else
        var_holding(i) = 0;
    end
end
var_holding = logical(var_holding);
ATP_highsalt = var_holding; %CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear var_holding
clear whichone
%% Pull out the columns of interest
pvalue = [FINAL_DATA.pvalue].';
experiment = {FINAL_DATA.experiment}.';
pvalue = [FINAL_DATA.pvalue].';
rsquared = [FINAL_DATA.rsquared].';
slope = [FINAL_DATA.slope].';
multiplemolecules = {FINAL_DATA.multiplemolecules}.';

%% restrictions 
pvaluesless = pvalue<0.05;
rsquaredgreated = rsquared>0.9;
restrictions = pvaluesless.*rsquaredgreated;
restrictions = logical(restrictions)
%% Multiple Molecules distinction

for i = 1:length(multiplemolecules)
    if isempty(multiplemolecules{i})
        multimol(i) = 0;
    elseif multiplemolecules{i} == 0;
        multimol(i) = 0;
    elseif multiplemolecules{i} == 1;
        multimol(i) = 1;
    end
end
   
multimol = logical(multimol);

%% restrictions with multimol
pvaluesless = pvalue<0.05;
rsquaredgreated = rsquared>0.9;
restrictions_andmultimol = pvaluesless.*rsquaredgreated*multimol;

restrictions_andmultimol = logical(restrictions)
%% Charged vs uncharged 
ATP_S = slope(ATP_70mMKCl);
ATP_p = pvalue(ATP_70mMKCl);
ATP_r = rsquared(ATP_70mMKCl);
ADP_S = slope(ADP);
ADP_p = pvalue(ADP);
ADP_r = rsquared(ADP);
noATP_S = slope(noATP);
noATP_p = pvalue(noATP);
noATP_r = rsquared(noATP);
ZB_noATP_S = slope(ZB_noATP);
ZB_noATP_p = pvalue(ZB_noATP);
ZB_noATP_r = rsquared(ZB_noATP);
ZB_noATP_mm = multimol(ZB_noATP);

sum(ATP_70mMKCl)
sum(noATP)
sum(ZB_noATP)
sum(ADP)

ATP_pvaluesless = ATP_p<0.05;
ATP_rsquaredgreated = ATP_r>0.8;
noATP_pvaluesless = noATP_p<0.05;
noATP_rsquaredgreated = noATP_r>0.8;
ZB_noATP_pvaluesless = ZB_noATP_p<0.05;
ZB_noATP_rsquaredgreated = ZB_noATP_r>0.8;
ADP_pvaluesless = ADP_p<0.05;
ADP_rsquaredgreated = ADP_r>0.8;

ATP_restr = ATP_pvaluesless.*ATP_rsquaredgreated;
sum(ATP_restr)
noATP_restr = noATP_pvaluesless.*noATP_rsquaredgreated;
sum(noATP_restr)
ZB_noATP_restr = ZB_noATP_pvaluesless.*ZB_noATP_rsquaredgreated;
sum(ZB_noATP_restr)
ZB_noATP_restr_mm = (ZB_noATP_pvaluesless.*ZB_noATP_rsquaredgreated).*ZB_noATP_mm';
sum(ZB_noATP_restr_mm)
ADP_restr = ADP_pvaluesless.*ADP_rsquaredgreated;
sum(ADP_restr)

ATP_restr = logical(ATP_restr);
noATP_restr = logical(noATP_restr);
ZB_noATP_restr = logical(ZB_noATP_restr);
ADP_restr = logical(ADP_restr);
ZB_noATP_restr_mm = logical(ZB_noATP_restr_mm);

close all
figure(1)



% histogram(ATP_S/2)
% hold on
% histogram(ATP_S(ATP_restr)/2,20,'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 0.20])
% hold off
% histogram(ATP_S(ATP_restr)/2,20,'Normalization', 'probability','DisplayStyle','bar',"BinLimits",[0 0.14])

% figure(2)
% histogram(ATP25_S/2)
hold on
% histogram(noATP_S(noATP_restr)/2,20,'Normalization', 'probability','DisplayStyle','bar',"BinLimits",[0 0.14])
% hold off


% figure(3)
% histogram(ATP200_S/2)
hold on
histogram(ZB_noATP_S(ZB_noATP_restr)/2,20,'Normalization', 'count','DisplayStyle','bar',"BinLimits",[0 0.14])
% Look at which at multimolecules
histogram(ZB_noATP_S(ZB_noATP_restr_mm)/2,20,'Normalization', 'count','DisplayStyle','bar',"BinLimits",[0 0.14])
% annotation('textbox', [0.75, 0.7, 0.15, 0.2],'FitBoxToText', 'on', 'String', [strcat("n = ", num2str(length(ADP_S(ADP_restr)))), newline, strcat("n = ", num2str(length(ZB_noATP_S(ZB_noATP_restr))))])
xticks([0 0.04 0.08 0.12 0.16 0.2])
box off 
hold off
% figure(2)
% % plot([25, 70, 200], [mean(ATP25_S(ATP25_restr)/2), mean(ATP_S(ATP_restr)/2), mean(ATP200_S(ATP200_restr)/2)])
% errorbar([25 70 200],[mean(ATP25_S(ATP25_restr)/2) mean(ATP_S(ATP_restr)/2) mean(ATP200_S(ATP200_restr)/2)],[std(ATP25_S(ATP25_restr)/2)/sqrt(length(ATP25_S(ATP25_restr)/2)) std(ATP_S(ATP_restr)/2)/sqrt(length(ATP_S(ATP_restr)/2)) std(ATP200_S(ATP200_restr)/2)/sqrt(length(ATP200_S(ATP200_restr)/2))],'o' )
% axis([0 200 0 0.055])
% yticks([0 0.01 0.02 0.03 0.04 0.05])
% box off
% fl = fitlm([25, 70, 200], [mean(ATP25_S(ATP25_restr)/2), mean(ATP_S(ATP_restr)/2), mean(ATP200_S(ATP200_restr)/2)])
% x = [25, 70, 200];
% y = [mean(ATP25_S(ATP25_restr)/2), mean(ATP_S(ATP_restr)/2), mean(ATP200_S(ATP200_restr)/2)];
% annotation('textbox', [0.65, 0.65, 0.15, 0.15],'FitBoxToText', 'on', 'String', [strcat("adjusted r2 is ", num2str(fl.Rsquared.Adjusted))])
% hold on
% plot(fl)

%% Salt Data
ATP_S = slope(ATP_70mMKCl);
ATP_p = pvalue(ATP_70mMKCl);
ATP_r = rsquared(ATP_70mMKCl);
ATP25_S = slope(ATP_25mMKCl);
ATP25_p = pvalue(ATP_25mMKCl);
ATP25_r = rsquared(ATP_25mMKCl);
ATP200_S = slope(ATP_200mMKCl);
ATP200_p = pvalue(ATP_200mMKCl);
ATP200_r = rsquared(ATP_200mMKCl);

sum(ATP_70mMKCl)
sum(ATP_25mMKCl)
sum(ATP_200mMKCl)

ATP_pvaluesless = ATP_p<0.05;
ATP_rsquaredgreated = ATP_r>0.8;
ATP25_pvaluesless = ATP25_p<0.05;
ATP25_rsquaredgreated = ATP25_r>0.8;
ATP200_pvaluesless = ATP200_p<0.05;
ATP200_rsquaredgreated = ATP200_r>0.8;

ATP_restr = ATP_pvaluesless.*ATP_rsquaredgreated;
sum(ATP_restr)
ATP25_restr = ATP25_pvaluesless.*ATP25_rsquaredgreated;
sum(ATP25_restr)
ATP200_restr = ATP200_pvaluesless.*ATP200_rsquaredgreated;
sum(ATP200_restr)

ATP_restr = logical(ATP_restr);
ATP25_restr = logical(ATP25_restr);
ATP200_restr = logical(ATP200_restr);

close all
figure(1)


% figure(2)
% histogram(ATP25_S/2)
hold on
histogram(ATP25_S(ATP25_restr)/2,'Normalization', 'pdf','DisplayStyle','stairs',"BinLimits",[0 0.2])
% hold off

% histogram(ATP_S/2)
% hold on
histogram(ATP_S(ATP_restr)/2,'Normalization', 'pdf','DisplayStyle','stairs',"BinLimits",[0 0.2])
% hold off
% figure(3)
% histogram(ATP200_S/2)
hold on
histogram(ATP200_S(ATP200_restr)/2,'Normalization', 'pdf','DisplayStyle','stairs',"BinLimits",[0 0.2])
annotation('textbox', [0.75, 0.7, 0.15, 0.2],'FitBoxToText', 'on', 'String', [strcat("n = ", num2str(length(ATP25_S(ATP25_restr)))), newline, strcat("n = ", num2str(length(ATP_S(ATP_restr)))),newline,strcat("n = ", num2str(length(ATP200_S(ATP200_restr))))])
xticks([0 0.04 0.08 0.12 0.16 0.2])
box off 
hold off
figure(2)
% plot([25, 70, 200], [mean(ATP25_S(ATP25_restr)/2), mean(ATP_S(ATP_restr)/2), mean(ATP200_S(ATP200_restr)/2)])
errorbar([25 70 200],[mean(ATP25_S(ATP25_restr)/2) mean(ATP_S(ATP_restr)/2) mean(ATP200_S(ATP200_restr)/2)],[std(ATP25_S(ATP25_restr)/2)/sqrt(length(ATP25_S(ATP25_restr)/2)) std(ATP_S(ATP_restr)/2)/sqrt(length(ATP_S(ATP_restr)/2)) std(ATP200_S(ATP200_restr)/2)/sqrt(length(ATP200_S(ATP200_restr)/2))],'o' )
axis([0 200 0 0.055])
yticks([0 0.01 0.02 0.03 0.04 0.05])
box off
% fl = fitlm([25, 70, 200], [mean(ATP25_S(ATP25_restr)/2), mean(ATP_S(ATP_restr)/2), mean(ATP200_S(ATP200_restr)/2)])
% x = [25, 70, 200];
% y = [mean(ATP25_S(ATP25_restr)/2), mean(ATP_S(ATP_restr)/2), mean(ATP200_S(ATP200_restr)/2)];
% annotation('textbox', [0.65, 0.65, 0.15, 0.15],'FitBoxToText', 'on', 'String', [strcat("adjusted r2 is ", num2str(fl.Rsquared.Adjusted))])
% hold on
% plot(fl)
%% Ploting histograms 
% experimental conditions: 
% 70mM KCl + 1mM ATP,                   ATP_70mMKCl      
% 25mM KCl + 1mM ATP,                   ATP_25mMKCl
% 200mM KCl + 1mM ATP,                  ATP_200mMKCl
% gammaS,                               gammaS
% dCas9,                                dCas9
% ADP,                                  ADP
% noATP,                                noATP
% "SWR1 70mM KCl 1mM ATP dual color",   ATP_70mMKCl_dualcolor
% "cas9 in dual color"                  dCas9_twocolor
% % close all
figure()
histogram(ATP_S(ATP_restr)/2,20, 'Normalization', 'count','DisplayStyle','bar',"BinLimits",[0 0.14])
title("ATP_70mMKCl vs noATP with restrictions")
xlabel("diffusion cofficient um^2/sec")
ylabel("counts")
annotation('textbox', [0.75, 0.75, 0.1, 0.1], 'String', strcat("n = ", num2str(sum(ATP_70mMKCl(restrictions)))))
% % % saveas(gcf,'ATP_70mMKCl with restrictions.png')
hold on
% figure()
histogram(slope(noATP(restrictions))/2,20, 'Normalization', 'count','DisplayStyle','bar',"BinLimits",[0 0.14])
% title("noATP with restrictions")
xlabel("diffusion cofficient um^2/sec")
ylabel("counts")
annotation('textbox', [0.65, 0.65, 0.1, 0.1], 'String', strcat("n = ", num2str(sum(noATP(restrictions)))))
% saveas(gcf,'ATP_70mMKCl with restrictions.png')
hold off

figure()
histogram(slope(ATP_70mMKCl)/2,20, 'Normalization', 'count','DisplayStyle','bar',"BinLimits",[0 0.14])
title("ATP_70mMKCl vs noATP without restrictions")
xlabel("diffusion cofficient um^2/sec")
ylabel("counts")
annotation('textbox', [0.75, 0.75, 0.1, 0.1], 'String', strcat("n = ", num2str(sum(ATP_70mMKCl))))
% saveas(gcf,'ATP_70mMKCl without restrictions.png')

hold on
% figure()
histogram(slope(noATP)/2,20, 'Normalization', 'count','DisplayStyle','bar',"BinLimits",[0 0.14])
% title("noATP without restrictions")
xlabel("diffusion cofficient um^2/sec")
ylabel("counts")
annotation('textbox', [0.65, 0.65, 0.1, 0.1], 'String', strcat("n = ", num2str(sum(noATP))))
% saveas(gcf,'ATP_25mMKCl without restrictions.png')
hold off 


% histogtram(slope(ATP_25mMKCl)/2,20, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 0.14])
% 
% histogram(slope(ATP_200mMKCl)/2,20, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 0.14])
% 
% histogram(slope(noATP)/2,20, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 0.14])
% 
% histogram(slope(gammaS)/2,20, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 0.14])
% 
% histogram(slope(ADP)/2,20, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 0.14])
% 
% histogram(slope(dCas9)/2,20, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 0.14])
% 
% histogram(slope(dCas9_twocolor)/2,20, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 0.14])
% 
% histogram(slope(ATP_70mMKCl_dualcolor)/2,20, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 0.14])
% 


%% pvalues 
% x1= pvalue(ATP_70mMKCl)';
% x2= pvalue(noATP)';
% g = [ones(size(x1))  2*ones(size(x2))];% 3*ones(size(A4)) 4*ones(size(A4))];
% label1 = "1mM ATP, Standard Salt";
% label2 = "No ATP, Standard Salt";
% figure
% boxplot([x1,x2],g,'Notch','on','Labels',{'ATP', 'noATP'})
% title('pvalue for linear fit - comparison between conditions')

% experimental conditions: 
% 70mM KCl + 1mM ATP,                   ATP_70mMKCl      
% 25mM KCl + 1mM ATP,                   ATP_25mMKCl
% 200mM KCl + 1mM ATP,                  ATP_200mMKCl
% gammaS,                               gammaS
% dCas9,                                dCas9
% ADP,                                  ADP
% noATP,                                noATP
% "SWR1 70mM KCl 1mM ATP dual color",   ATP_70mMKCl_dualcolor
% "cas9 in dual color"                  dCas9_twocolor

figure()
histogram(pvalue(noATP),20, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 1])
title("noATP")
xlabel("pvalue")
ylabel("counts")
saveas(gcf,'noATP p-value.png')

figure()
histogram(pvalue(ATP_70mMKCl),20, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 1])
title("ATP_70mMKCl")
xlabel("pvalue")
ylabel("counts")
saveas(gcf,'ATP_70mMKCl p-value.png')

figure()
histogram(pvalue(ADP),20, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 1])
title("ADP")
xlabel("pvalue")
ylabel("counts")
saveas(gcf,'ADP p-value.png')

figure()
histogram(pvalue(gammaS),20, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 1])
title("gammaS")
xlabel("pvalue")
ylabel("counts")
saveas(gcf,'gammaS p-value.png')

figure()
histogram(pvalue(dCas9),20, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 1])
title("dCas9")
xlabel("pvalue")
ylabel("counts")
saveas(gcf,'dCas9 p-value.png')

figure()
histogram(pvalue(dCas9_twocolor),20, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 1])
title("dCas9_twocolor")
xlabel("pvalue")
ylabel("counts")
saveas(gcf,'dCas9_twocolor p-value.png')

figure()
histogram(pvalue(ATP_70mMKCl_dualcolor),20, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 1])
title("ATP_70mMKCl_dualcolor")
xlabel("pvalue")
ylabel("counts")
saveas(gcf,'ATP_70mMKCl_dualcolor p-value.png')

figure()
histogram(pvalue(ATP_25mMKCl),20, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 1])
title("ATP_25mMKCl")
xlabel("pvalue")
ylabel("counts")
saveas(gcf,'ATP_25mMKCl p-value.png')

figure()
histogram(pvalue(ATP_200mMKCl),20, 'Normalization', 'pdf','DisplayStyle','bar',"BinLimits",[0 1])
title("ATP_200mMKCl")
xlabel("pvalue")
ylabel("counts")
saveas(gcf,'ATP_200mMKCl p-value.png')

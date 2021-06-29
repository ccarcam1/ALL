%% First add new data to the end of the data frame
a =length(FINAL_DATA);
% a=0;
for i = 1:length(Green_MY_DATA)
    FINAL_DATA(i+a).name = Green_MY_DATA(i).name;
    FINAL_DATA(i+a).experiment = "Swc2 150mM noATP";
    FINAL_DATA(i+a).line_time = Green_MY_DATA(i).line_time;
    FINAL_DATA(i+a).particle_tracked = Green_MY_DATA(i).particle_tracked;   
    FINAL_DATA(i+a).MSD = Green_MY_DATA(i).MSD;
    FINAL_DATA(i+a).crop_coords = Green_MY_DATA(i).crop_coordinates;
    FINAL_DATA(i+a).crop = Green_MY_DATA(i).crop;
    FINAL_DATA(i+a).kymograph = Green_MY_DATA(i).kymograph;
    FINAL_DATA(i+a).pvalue = Green_MY_DATA(i).pval;
    FINAL_DATA(i+a).rsquared = Green_MY_DATA(i).r2;
    FINAL_DATA(i+a).yintercept = Green_MY_DATA(i).yinter;
    if ~isnan(Green_MY_DATA(i).fit_cutoffis)
         FINAL_DATA(i+a).fitcutoff = Green_MY_DATA(i).fit_cutoffis;
    else 
        FINAL_DATA(i+a).fitcutoff = round(length(Green_MY_DATA(i).particle_tracked)/4);
    end
    FINAL_DATA(i+a).fitcutoff_NaN = Green_MY_DATA(i).fit_cutoffis;
    FINAL_DATA(i+a).slope = Green_MY_DATA(i).slope;
    FINAL_DATA(i+a).multiplemolecules = Green_MY_DATA(i).multiplemolecules;
end
%% Length of traces

% for i = 1:length(FINAL_DATA)
%     lengthis(i) = (size(FINAL_DATA(i).particle_tracked, 1)*FINAL_DATA(i).line_time)/1000;
% end


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
% "SWR1 + ZB + 1mM ATP"                 ZB_ATP
% "SWR1 + ZB + 1mM ATP mixed before"    ZB_ATP_mixed

%% Make logical pointers for data
% experiment = {FINAL_DATA.experiment}.';

% whichone = "70mM KCl + 1mM ATP";% CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% whichone = "25mM KCl + 1mM ATP";% CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% whichone = "200mM KCl + 1mM ATP";% CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% whichone = "gammaS";% CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% whichone = "dCas9";% CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% whichone = "ADP";% CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
whichone = "noATP";% CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% whichone = "SWR1 70mM KCl 1mM ATP dual color";% CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% whichone = "cas9 in dual color";% CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% whichone = "SWR1 + ZB + noATP";% CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% whichone = "SWR1 + ZB + 1mM ATP";% CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% whichone = "SWR1 + ZB + 1mM ATP mixed before";% CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var_holding = [];

for i = 1:length(experiment)
    if matches(experiment{i}, whichone)
        var_holding(i) = 1; 
    else
        var_holding(i) = 0;
    end
end
var_holding = logical(var_holding);
% ATP_70mMKCl = var_holding; %CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATP_25mMKCl = var_holding; %CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATP_200mMKCl = var_holding; %CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gammaS = var_holding; %CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dCas9 = var_holding; %CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADP = var_holding; %CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noATP = var_holding; %CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATP_70mMKCl_dualcolor = var_holding; %CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dCas9_twocolor = var_holding; %CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZB_noATP = var_holding; %CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZB_ATP = var_holding; %CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZB_ATP_mixed = var_holding; %CHANGE THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% pvaluesless = pvalue<0.05;
% rsquaredgreated = rsquared>0.9;
% restrictions = pvaluesless.*rsquaredgreated;
% restrictions = logical(restrictions);

restrictions = pvalue<0.05 & rsquared>0.9;

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
% pvaluesless = pvalue<0.05;
% rsquaredgreated = rsquared>0.9;
% restrictions_andmultimol = pvaluesless.*rsquaredgreated*multimol;
% restrictions_andmultimol = logical(restrictions_andmultimol);

restrictions_andmultimol = pvalue<0.05 & rsquared>0.9 & ~multimol;

%% get rid of short traces 

l_rest = lengthis> 5;
% restrictions_andl = (pvaluesless.*rsquaredgreated).*l_rest';
% restrictions_andl = logical(restrictions_andl);
restrictions_len = pvalue<0.05 & rsquared>0.9 & l_rest;

%% Get Specific Points User Logical Indexing (Part1) 
% ATP_S = slope(ATP_70mMKCl);
% ATP_p = pvalue(ATP_70mMKCl);
% ATP_r = rsquared(ATP_70mMKCl);
% ATP_l = restrictions_andl(ATP_70mMKCl);
% 
% % ADP_S = slope(ADP);
% % ADP_p = pvalue(ADP);
% % ADP_r = rsquared(ADP);
% 
% noATP_S = slope(noATP);
% noATP_p = pvalue(noATP);
% noATP_r = rsquared(noATP);
% noATP_l = restrictions_andl(noATP);
% 
% 
% ZB_noATP_S = slope(ZB_noATP);
% ZB_noATP_p = pvalue(ZB_noATP);
% ZB_noATP_r = rsquared(ZB_noATP);
% ZB_noATP_mm = multimol(ZB_noATP);
% 
% ZB_ATP_S = slope(ZB_ATP);
% ZB_ATP_p = pvalue(ZB_ATP);
% ZB_ATP_r = rsquared(ZB_ATP);
% ZB_ATP_mm = multimol(ZB_ATP);

ZB_ATP_mixed_S = slope(ZB_ATP_mixed);
ZB_ATP_mixed_p = pvalue(ZB_ATP_mixed);
ZB_ATP_mixed_r = rsquared(ZB_ATP_mixed);
ZB_ATP_mixed_mm = multimol(ZB_ATP_mixed);

% sum(ATP_70mMKCl)
% sum(noATP)
% % sum(ADP)
% sum(ZB_noATP)
% sum(ZB_ATP)
sum(ZB_ATP_mixed)
%% Part 2
% ATP_pvaluesless = ATP_p<0.05;
% ATP_rsquaredgreated = ATP_r>0.8;
% 
% noATP_pvaluesless = noATP_p<0.05;
% noATP_rsquaredgreated = noATP_r>0.8;
% 
% ZB_noATP_pvaluesless = ZB_noATP_p<0.05;
% ZB_noATP_rsquaredgreated = ZB_noATP_r>0.8;
% 
% ZB_ATP_pvaluesless = ZB_ATP_p<0.05;
% ZB_ATP_rsquaredgreated = ZB_ATP_r>0.8;

ZB_ATP_mixed_pvaluesless = ZB_ATP_mixed_p<0.05;
ZB_ATP_mixed_rsquaredgreated = ZB_ATP_mixed_r>0.8;

% ADP_pvaluesless = ADP_p<0.05;
% ADP_rsquaredgreated = ADP_r>0.8;

%% Part 3

% ATP_restr = ATP_pvaluesless.*ATP_rsquaredgreated;
% sum(ATP_restr)
% ATP_restr_l = (ATP_pvaluesless.*ATP_rsquaredgreated).*ATP_l;
% sum(ATP_restr_l)
% 
% noATP_restr = noATP_pvaluesless.*noATP_rsquaredgreated;
% sum(noATP_restr)
% noATP_restr_l = (noATP_pvaluesless.*noATP_rsquaredgreated).*noATP_l;
% sum(noATP_restr_l)
% 
% ZB_noATP_restr = ZB_noATP_pvaluesless.*ZB_noATP_rsquaredgreated;
% sum(ZB_noATP_restr)
% 
% ZB_noATP_restr_mm = (ZB_noATP_pvaluesless.*ZB_noATP_rsquaredgreated).*ZB_noATP_mm';
% sum(ZB_noATP_restr_mm)
% 
% ZB_ATP_restr = ZB_ATP_pvaluesless.*ZB_ATP_rsquaredgreated;
% sum(ZB_ATP_restr)
% 
% ZB_ATP_restr_mm = (ZB_ATP_pvaluesless.*ZB_ATP_rsquaredgreated).*ZB_ATP_mm';
% sum(ZB_ATP_restr_mm)

ZB_ATP_mixed_restr = ZB_ATP_mixed_pvaluesless.*ZB_ATP_mixed_rsquaredgreated;
sum(ZB_ATP_mixed_restr)

% ZB_ATP_mixed_restr_mm = (ZB_ATP_mixed_pvaluesless.*ZB_ATP_mixed_rsquaredgreated).*ZB_ATP_mixed_mm';
% sum(ZB_ATP_mixed_restr_mm)
% 
% % ADP_restr = ADP_pvaluesless.*ADP_rsquaredgreated;
% % sum(ADP_restr)
% 
% ATP_restr = logical(ATP_restr);
% noATP_restr = logical(noATP_restr);
% % ADP_restr = logical(ADP_restr);
% ZB_noATP_restr = logical(ZB_noATP_restr);
% ZB_noATP_restr_mm = logical(ZB_noATP_restr_mm);
% ZB_ATP_restr = logical(ZB_ATP_restr);
% ZB_ATP_restr_mm = logical(ZB_ATP_restr_mm);
ZB_ATP_mixed_restr = logical(ZB_ATP_mixed_restr);
% ZB_ATP_mixed_restr_mm = logical(ZB_ATP_mixed_restr_mm);
%% Histogram Plot
% close all
figure(6)
Normalization_is = 'pdf';
Binlimis = [0 0.14]; 
numberofbins = 30;
hold on

% histogram(noATP_S(noATP_restr)/2,numberofbins,'Normalization', Normalization_is,'DisplayStyle','bar',"BinLimits",Binlimis, "DisplayName", ['no ATP || n = ', num2str(length(noATP_S(noATP_restr)))])
% histogram(ATP_S(ATP_restr)/2,numberofbins,'Normalization', Normalization_is,'DisplayStyle','bar',"BinLimits",Binlimis, "DisplayName", ['ATP || n = ', num2str(length(ATP_S(ATP_restr)))])
% histogram(ZB_noATP_S(ZB_noATP_restr)/2,numberofbins,'Normalization', Normalization_is,'DisplayStyle','bar',"BinLimits",Binlimis, "DisplayName", ['ZB noATP || n = ', num2str(length(ZB_noATP_S(ZB_noATP_restr)))])
% histogram(ZB_ATP_S(ZB_ATP_restr)/2,numberofbins,'Normalization', Normalization_is,'DisplayStyle','bar',"BinLimits",Binlimis, "DisplayName", ['ZB 1mM ATP || n = ', num2str(length(ZB_ATP_S(ZB_ATP_restr)))])
histogram(ZB_ATP_mixed_S(ZB_ATP_mixed_restr)/2,numberofbins,'Normalization', Normalization_is,'DisplayStyle','bar',"BinLimits",Binlimis, "DisplayName", ['ZB 1mM ATP mixed || n = ', num2str(length(ZB_ATP_mixed_S(ZB_ATP_mixed_restr)))])

% Look at which at multimolecules
% histogram(ZB_noATP_S(ZB_noATP_restr_mm)/2,numberofbins,'Normalization', Normalization_is,'DisplayStyle','bar',"BinLimits",Binlimis, "DisplayName", ['ZB noATP sm|| n = ', num2str(length(ZB_noATP_S(ZB_noATP_restr_mm)))])
% histogram(ZB_ATP_S(ZB_ATP_restr_mm)/2,numberofbins,'Normalization', Normalization_is,'DisplayStyle','bar',"BinLimits",Binlimis, "DisplayName", ['ZB 1mM ATP sm|| n = ', num2str(length(ZB_ATP_S(ZB_ATP_restr_mm)))])

% Plot without multiple molecules in analysis
% histogram(ZB_noATP_S(~ZB_noATP_restr_mm)/2,numberofbins,'Normalization', Normalization_is,'DisplayStyle','bar',"BinLimits",Binlimis, "DisplayName", ['ZB noATP sm|| n = ', num2str(length(ZB_noATP_S(~ZB_noATP_restr_mm)))])
% histogram(ZB_ATP_S(~ZB_ATP_restr_mm)/2,numberofbins,'Normalization', Normalization_is,'DisplayStyle','bar',"BinLimits",Binlimis, "DisplayName", ['ZB 1mM ATP sm|| n = ', num2str(length(ZB_ATP_S(~ZB_ATP_restr_mm)))])
% histogram(ZB_ATP_mixed_S(~ZB_ATP_mixed_restr_mm)/2,numberofbins,'Normalization', Normalization_is,'DisplayStyle','bar',"BinLimits",Binlimis, "DisplayName", ['ZB 1mM ATP mixed sm|| n = ', num2str(length(ZB_ATP_mixed_S(~ZB_ATP_mixed_restr_mm)))])

% Plot without short traces 
% histogram(noATP_S(~noATP_restr_l)/2,numberofbins,'Normalization', Normalization_is,'DisplayStyle','bar',"BinLimits",Binlimis, "DisplayName", ['noATP l>2|| n = ', num2str(length(noATP_S(~noATP_restr_l)))])
% histogram(ATP_S(~ATP_restr_l)/2,numberofbins,'Normalization', Normalization_is,'DisplayStyle','bar',"BinLimits",Binlimis, "DisplayName", ['ATP l>2|| n = ', num2str(length(ATP_S(~ATP_restr_l)))])

legend
hold off


%%

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
% "SWR1 + ZB + 1mM ATP"                 ZB_ATP
% "SWR1 + ZB + 1mM ATP mixed before"    ZB_ATP_mixed
a = 1;
for i = 1:length(FINAL_DATA2)
    if matches(FINAL_DATA2(i).experiment, "70mM KCl + 1mM ATP")| ...
            matches(FINAL_DATA2(i).experiment, "25mM KCl + 1mM ATP")| ...
            matches(FINAL_DATA2(i).experiment, "200mM KCl + 1mM ATP")| ...
            matches(FINAL_DATA2(i).experiment, "gammaS")| ...
            matches(FINAL_DATA2(i).experiment, "dCas9")| ...
            matches(FINAL_DATA2(i).experiment, "ADP")| ...
            matches(FINAL_DATA2(i).experiment, "noATP")| ...
            matches(FINAL_DATA2(i).experiment, "SWR1 70mM KCl 1mM ATP dual color")| ...
            matches(FINAL_DATA2(i).experiment, "SWR1 + ZB + noATP")| ...
            matches(FINAL_DATA2(i).experiment, "SWR1 + ZB + 1mM ATP")| ...
            matches(FINAL_DATA2(i).experiment, "SWR1 + ZB + 1mM ATP mixed before")
        newone(a) = FINAL_DATA2(i);
        a = a+1;
    end
end





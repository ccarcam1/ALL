% function cutoff = myCutOffFinder(MSDs)
function [rsquared,pvalue,yintercept,slopeval,fit_cutoffis] = myCutOffFinder(MSDs)
tic
l = size(MSDs,2);
rsquared = cell(l, 1);
pvalue = cell(l, 1);
yintercept = cell(l, 1);
slopeval = cell(l, 1);
fit_cutoffis = zeros(1,l);
parfor k = 1:l
    disp(k)
    rsquared_holding = zeros(1,size((3:(size(MSDs(k).coords,1)-1)),2)); % potentially these need to be corrected to be 4 instead of 3
    pvalue_holding = zeros(1,size((3:(size(MSDs(k).coords,1)-1)),2));
    yintercept_holding = zeros(1,size((3:(size(MSDs(k).coords,1)-1)),2));
    slopeval_holding = zeros(1,size((3:(size(MSDs(k).coords,1)-1)),2));
    for i = 4:(size(MSDs(k).coords,1)-1)   
        x = MSDs(k).MSD(1:i, 1);
        y = MSDs(k).MSD(1:i, 2);
        mdl = fitlm(x,y);
        rsquared_holding(i) = mdl.Rsquared.Ordinary;
        pvalue_holding(i) = mdl.Coefficients.pValue(2);
        yintercept_holding(i) = mdl.Coefficients.Estimate(1);
        slopeval_holding(i) = mdl.Coefficients.Estimate(2);
    end
    rsquared{k} = rsquared_holding;
    pvalue{k} = pvalue_holding;
    yintercept{k} = yintercept_holding;
    slopeval{k} = slopeval_holding;
%     cutoff(k).rsquared = rsquared;
%     cutoff(k).pvalue = pvalue;
    for i = 5:length(rsquared{k})-5
        a = rsquared{k}(i+1)-rsquared{k}(i);
        aa = rsquared{k}(i)-rsquared{k}(i-1);
        if a<0 && pvalue{k}(i)<0.01 && aa>0
            b = rsquared{k}(i+2)-rsquared{k}(i+1);
            c = rsquared{k}(i+3)-rsquared{k}(i+2);
            d = rsquared{k}(i+4)-rsquared{k}(i+3);
            e = rsquared{k}(i+5)-rsquared{k}(i+4);
            if b<0 && c<0 && d<0 && e<0
                fit_cutoffis(k) = i-1;
%                 cutoff(k).fit_cutoff = i-1;
                break
            else
            end
        else
        fit_cutoffis(k) = NaN;
%         cutoff(k).fit_cutoff = NaN;
        end
    end
end
toc
end




% for k = 1:20 %length(mydata)
% MSD = MSDss(k).MSD;
% record_this = zeros(size(mydata(k).coords, 1)-1);
% record_that = zeros(size(mydata(k).coords, 1)-1);
% cuthere = round(length(mydata(k).coords)/4);
%     for i = 3: cuthere %size(mydata(k).coords, 1)-1
%     x = MSD(1:i, 1);
%     y = MSD(1:i, 2);
%     mdl = fitlm(x,y);
%     record_this(i) = mdl.Rsquared.Ordinary;
%     record_that(i) = mdl.Coefficients.pValue(2);
%     end
% holdthis(k).rsquared = record_this;
% holdthis(k).pvalue = record_that;
% disp(k)
% end
% 
% for k = 1:20%length(mydata)
% %     x = 1;
%     r2_list = (holdthis(k).rsquared);
%     for i = 2:length(r2_list)-5
%         a = r2_list(i+1)-r2_list(i);
%         aa = r2_list(i)-r2_list(i-1);
%         if a<0 && holdthis(k).pvalue(i)<0.01 && aa>0
%             b = r2_list(i+2)-r2_list(i+1);
%             c = r2_list(i+3)-r2_list(i+2);
%             d = r2_list(i+4)-r2_list(i+3);
%             e = r2_list(i+5)-r2_list(i+4);
%             if b<0 && c<0 && d<0 && e<0
%                 cutoff(k).fit_cutoff = i-1;
%                 break
%             else
%             end
%         else
%         cutoff(k).fit_cutoff = NaN;
%         end
%     end
% end




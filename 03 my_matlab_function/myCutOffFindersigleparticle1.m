%% Claudia Carcamo 2020-04-20 
% Description:  This finds the first tau in a particle's MSD vs time lag
% plot where the linear fit is the best.  It does this my maximizing for
% both the p-value and the r-squared value of the fit. It does this by
% finding the first tau in which fits of tau + 1,2,3,4,&5 have worse fits.
% By doing this, the fit is typically over a small percentage of the full
% MSD plot. Likely <10% of the whole trace, which represents the brownian
% region. 
function [cutoff, NaNsAre] = myCutOffFindersigleparticle(MSD_data)
MSD_timelag= MSD_data.MSD(:,1);
MSD_MSD = MSD_data.MSD(:,2);
    for i = 3:(length(MSD_timelag)-1)   
        x = MSD_timelag(1:i);
        y = MSD_MSD(1:i);
        mdl = fitlm(x,y);
        rsquared(i) = mdl.Rsquared.Ordinary;
        pvalue(i) = mdl.Coefficients.pValue(2);
    end
    for i = 2:length(rsquared)-5
        a = rsquared(i+1)-rsquared(i);
        aa = rsquared(i)-rsquared(i-1);
        if a<0 && pvalue(i)<0.01 && aa>0
            b = rsquared(i+2)-rsquared(i+1);
            c = rsquared(i+3)-rsquared(i+2);
            d = rsquared(i+4)-rsquared(i+3);
            e = rsquared(i+5)-rsquared(i+4);
            if b<0 && c<0 && d<0 && e<0
                cutoff = i-1;            
                NaNsAre= "_";
                break
            else
            end
        else
        cutoff = length(MSD_timelag)/4;
        NaNsAre= MSD_data.molID;
        end
    end
end


% 
% function [cutoff, NaNsAre] = myCutOffFindersigleparticle(MSD_timelag, MSD_MSD)
%     for i = 3:(length(MSD_timelag)-1)   
%         x = MSD_timelag(1:i);
%         y = MSD_MSD(1:i);
%         mdl = fitlm(x,y);
%         rsquared(i) = mdl.Rsquared.Ordinary;
%         pvalue(i) = mdl.Coefficients.pValue(2);
%     end
%     for i = 2:length(rsquared)-5
%         a = rsquared(i+1)-rsquared(i);
%         aa = rsquared(i)-rsquared(i-1);
%         if a<0 && pvalue(i)<0.01 && aa>0
%             b = rsquared(i+2)-rsquared(i+1);
%             c = rsquared(i+3)-rsquared(i+2);
%             d = rsquared(i+4)-rsquared(i+3);
%             e = rsquared(i+5)-rsquared(i+4);
%             if b<0 && c<0 && d<0 && e<0
%                 cutoff = i-1;
%                 break
%             else
%             end
%         else
%         cutoff = NaN;
%         end
%     end
% end

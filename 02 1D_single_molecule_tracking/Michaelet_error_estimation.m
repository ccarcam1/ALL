%% Left over script from the Michelet error estimation 


% % MSDs = myAllDispFun(MSDs);
% % MSDs = mysampleSDev(MSDs);
% % MSDs = myMichaletFun_guesses(MSDs);
% % load('MSDs.mat');
% % MSDs = myMichaletFun(MSDs);
% % le = length(MSDs);
% % slope = zeros(le,1);
% % intercept = zeros(le,1);
% % x = zeros(le,1);
% % for i = 1:le
% %     MSDs(i).MSD_MichaeletError(end) = MSDs(i).MSD_MichaeletError(end-1);
% %     weights = 1-(MSDs(i).MSD_MichaeletError/max(MSDs(i).MSD_MichaeletError))';
% %     MSDs(i).weights = weights;
% %     Slope_lscov(i) = lscov(MSDs(i).MSD(:,1),MSDs(i).MSD(:,2), weights);
% % end

%% Weighted Linear Fitting based on Localization Uncertainty and MSD (Michaelet paper 2011) 
% calculate weighted least mean squared linear fit to all data
% x = zeros(length(MSDs),1);
le = length(MSDs);
slope = zeros(le,1);
intercept = zeros(le,1);
x = zeros(le,1);
for i = 1:le
    MSDss(i).MSD_MichaeletError(end) = MSDss(i).MSD_MichaeletError(end-1);
    weights = 1-(MSDss(i).MSD_MichaeletError/max(MSDss(i).MSD_MichaeletError))';
    MSDss(i).weights = weights;
    x_MichaeletError(i) = lscov(MSDss(i).MSD(:,1),MSDss(i).MSD(:,2), weights);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    weights_sampleSDev = 1-(MSDss(i).sampleSDev/max(MSDss(i).sampleSDev))';
    MSDss(i).weights_sampleSDev = weights_sampleSDev;
    x_sampleSDev(i) = lscov(MSDss(i).MSD(:,1),MSDss(i).MSD(:,2), weights_sampleSDev);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    weightss = 1-(MSDs(i).MSD(:,3)/max(MSDs(i).MSD(:,3)))';
    MSDs(i).weightss = weightss;
    x_SEM(i) = lscov(MSDs(i).MSD(:,1),MSDs(i).MSD(:,2), weightss);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mdl = fitlm(MSDss(i).MSD(1:MSDss(i).cutoff_was,1),MSDss(i).MSD(1:MSDss(i).cutoff_was,2));
    slope(i) = mdl.Coefficients.Estimate(2);
    intercept(i) = mdl.Coefficients.Estimate(1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mdl = fitlm(MSDss(i).MSD(:,1),MSDss(i).MSD(:,2));
    slope_full(i) = mdl.Coefficients.Estimate(2);
    intercept_full(i) = mdl.Coefficients.Estimate(1);
    disp(i);
end

% taus = MSDs(1).MSD(:,1);
% MSD = MSDs(1).MSD(:,2);
% MSDs(1).MSD_MichaeletError(end) = MSDs(1).MSD_MichaeletError(end-1);
% weights = 1-(MSDs(1).MSD_MichaeletError/max(MSDs(1).MSD_MichaeletError))';

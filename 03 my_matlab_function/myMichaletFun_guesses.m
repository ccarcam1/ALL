function MSD_data  = myMichaletFun_guesses(MSD_data)
% plot(MSDss(1).MSD(:,1),(MSDss(1).MSD(:,2)))
l = length(MSD_data);
alpha_val = zeros(l,1);
epsilon_val = zeros(l,1);
NaNsWere = {};
parfor i = 1:l
%     end_point = round(length(MSD_data(i))/4);
%     cutoff = myCutOffFindersigleparticle(MSD_data(i).MSD(:,1), MSD_data(i).MSD(:,2));
    [cutoff, NaNsAre] = myCutOffFindersigleparticle(MSD_data(i));
    x = MSD_data(i).MSD(1:cutoff,1);
    y = MSD_data(i).MSD(1:cutoff,2);
    mdl = fitlm(x,y);
    MSD_data(i).alpha_val = mdl.Coefficients.Estimate(2)*(MSD_data(i).timestep/1000);
    MSD_data(i).epsilon_val = mdl.Coefficients.Estimate(1);
    MSD_data(i).cutoff_was = cutoff;
    MSD_data(i).NaNsWere = NaNsAre;
    disp(['guess# ', num2str(i)]);
end
end

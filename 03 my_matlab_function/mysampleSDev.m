%% Claudia Carcamo 2020-04-21
% Description: This function calculates the sample standard deviation for
% the MSD data
function MSDss = mysampleSDev(MSDss)
for i = 1:length(MSDss)
    sampleSDev = [];
    for j = 1:length(MSDss(i).breakdown)
        sampleSDev(j) = sqrt(sum((MSDss(i).breakdown{j}(:,2)- MSDss(i).MSD(j,2)).^2)/(length(MSDss(i).breakdown{j})-1));
    end
    MSDss(i).sampleSDev = sampleSDev;
    disp(i)
end
end
       
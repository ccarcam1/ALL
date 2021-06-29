%% Claudia Carcamo 2020-04-21
% Description: Breaks r_all into the individual time-lag displacements,
% this was necessary for analysis by mysampleSDev function, which finds the
% sample standard deviation for the MSD curve
function MSDss = myAllDispFun(MSDss)
for i = 1:length(MSDss)
    a = cell([1 length(MSDss(i).MSD)]);
    for j = 1:length(MSDss(i).MSD)
       a{j}=MSDss(i).r_all((MSDss(i).r_all(:,1) == j),:);
    end
    MSDss(i).breakdown = a;
    disp(['cell = ', num2str(i)])
end
end
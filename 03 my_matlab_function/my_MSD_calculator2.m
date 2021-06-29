function MSDs = my_MSD_calculator2(mydata)
x= 1;
pixSize = 0.1;
for i = 1:length(mydata)  
    n=size(mydata(i).particle_tracked(:,2), 1);
    r_all=[];
    for j = 1:n-1 % i is the frame lag
        r=zeros(n-j,2);
        for k = 1: n-j % n-i is the total number of lags one can calculate for frame lag t
            r(k, 1)=mydata(i).particle_tracked(k+j,1)-mydata(i).particle_tracked(k,1); % convert frame lag to time lag
            dif = (pixSize) *(mydata(i).particle_tracked(k+j, 2)-mydata(i).particle_tracked(k, 2));%pixSize/1000 was removed
            r(k, 2)=sum(dif.^2); % sqaure displacement
        end
        r_all(end+1:end+size(r, 1), :)=r;
    end
    MSDs(x).r_all = r_all; 
            x = x+1; 
    disp(i); 
end

for i = 1:length(mydata)
    tau = mydata(i).line_time/1000;
    d = MSDs(i).r_all;
    time_lags = unique(d(:, 1));
    n=size(mydata(i).particle_tracked, 1);
    D_all = zeros(n-1,2);
% finds the mean displacement squared from above displacements
    for j = 1:n-1
        ind=find(d(:, 1) == time_lags(j));
        D_all(j, 1) = tau * time_lags(j); % converts to real time lag
        D_all(j, 2) = mean(d(ind, 2)); %MSD in x-direction for timelag i
        D_all(j, 3) = std(d(ind, 2))/sqrt(length(ind)); %sem in x-direction for timelag i***********************************D13 
% %         D_all(i, 4) = std(d(ind, 2)); % std
%         D_all(j, 5) = length(d(ind, 2)); % num points
    end
    MSDs(i).MSD = D_all;
    MSDs(i).coords = mydata(i).particle_tracked(:,2);
    MSDs(i).molID = mydata(i).name;
    MSDs(i).timestep = mydata(i).line_time;
    
end
end

tracking_info = particle_tracked(noATP);
scan_time = line_time(noATP);
%%
x= 1;
pixSize = 0.1;
for i = 1:length(tracking_info)  
    n=size(tracking_info{i}(:,2), 1);
    r_all=[];
    for j = 1:n-1 % i is the frame lag
        r=zeros(n-j,2);
        for k = 1: n-j % n-i is the total number of lags one can calculate for frame lag t
            r(k, 1)=tracking_info{i}(k+j,1)-tracking_info{i}(k,1); % convert frame lag to time lag
            dif = (pixSize) *(tracking_info{i}(k+j, 2)-tracking_info{i}(k, 2));%pixSize/1000 was removed
            r(k, 2)=sum(dif.^2); % sqaure displacement
        end
        r_all(end+1:end+size(r, 1), :)=r;
    end
    MSDs_hold{x} = r_all; 
            x = x+1; 
    disp(i); 
end

for i = 1:length(tracking_info)
    tau = scan_time(i)/1000;
    d = MSDs_hold{i};
    time_lags = unique(d(:, 1));
    n=size(tracking_info{i}, 1);
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
    MSDs{i} = D_all;
%     coords{i} = mydata(i).particle_tracked(:,2);
%     MSDs(i).molID = mydata(i).name;
%     MSDs(i).timestep = mydata(i).line_time;
    
end

%%
holdingthat = [];
for i = 1:length(MSDs_hold)
    holdingthis{i}=MSDs_hold{i}((MSDs_hold{i}(:, 1) == 5),2);
    holdingthat = vertcat(holdingthat, holdingthis{i});
%     holdingthat = vertcat(holdingthat, 1000*holdingthis{i}/(2*scan_time(i)));
end
    

%%

tracking_info = particle_tracked(noATP);
scan_time = line_time(noATP);
pixSize = 0.1;
counter = 1;
a=[];;
for i = 1:length(tracking_info)
    test_track = tracking_info{i};
    test_time = line_time(i);
    for j = 1:length(test_track)-5
        a(counter, 1)=test_track(j+5,1)-test_track(j,1); 
        dif = (pixSize) *(test_track(j+5, 2)- test_track(j, 2));
        a(counter, 2)=sum(dif.^2);
        counter = counter +1;
    end
end

%%
tracking_info = particle_tracked(noATP);
scan_time = line_time(noATP);
for xyz = 1:length(tracking_info)
test_track = tracking_info{xyz};
test_time = line_time(xyz);

pixSize = 0.1;
n=size(test_track, 1);
windows = n-1;
storage = cell(windows, 1);
storage(:,1)={0};
storage2 = zeros(15, 2);
subtrack = {};
counter = 1;
counter2 =1;
for k = 1:n-5
    for m = 1:5
        subtrack{k,1} = test_track(k:k+m,:);
    end
end


for lmnop = 1:length(subtrack)
    counter = 1;
    storage2 = zeros(15, 2);
for j = 1:5
    for i = 1:6-j
        storage2(counter, 1)=subtrack{lmnop, 1}(i+j,1)-subtrack{lmnop, 1}(i,1); 
        dif = (pixSize) *(subtrack{lmnop, 1}(i+j, 2)-subtrack{lmnop, 1}(i, 2));
        storage2(counter, 2)=sum(dif.^2);
        counter = counter +1;
    end
end
subtrack{lmnop, 2} = storage2;
end

for lmnop = 1:length(subtrack)
for i = 1:5
    subtrack{lmnop, 3}(i) = mean(subtrack{lmnop, 2}((subtrack{lmnop, 2}(:,1) == i),2));
end
end
a{xyz} = subtrack;
end

%%
% figure(1)
hold on
for j = 1:length(a)
for i = 1:length(a{j})
%     plot([1,2,3,4,5]*test_time,subtrack{i,3})
%     plot([1,2]*test_time, subtrack{i,3}(1:2), "LineWidth", 2)
    holdingplace(i) = ((a{j}{i,3}(2)-a{j}{i,3}(1))*pixSize)/(1*(test_time/1000));
    holdingplacex(i) = ((a{j}{i,3}(2)-0)*pixSize)/(2*(test_time/1000));

%     holdingplace2(i) = ((a{j}{i,3}(3)-a{j}{i,3}(1))*pixSize)/(2*(test_time/1000));
%     holdingplace3(i) = ((a{j}{i,3}(4)-a{j}{i,3}(1))*pixSize)/(3*(test_time/1000));

%     p = polyfit([1,2],subtrack{i,3}(1:2),1);
%     f = polyval(p, [1,2]);
%     plot([1,2],f, "LineWidth", 2)
end
end
figure(1)
% hold on
hold on
histogram(holdingplace/2)
figure(2)
hold on
histogram(holdingplacex/2)
% histogram(holdingplace2/2)
% histogram(holdingplace3/2)

% hold off

%%

tracking_info = particle_tracked(noATP);
scan_time = line_time(noATP);

test_track = tracking_info{1};
test_time = line_time(1);

%%
pixSize = 0.1;
n=size(test_track, 1);
windows = n-1;
storage = cell(windows, 1);
storage(:,1)={0};
storage2 = zeros(15, 2);

% for i = 1:windows
%     for j = i: i + 5
%             for k = 1:5
%                 storage2(i,k) = j(i);
%             end    
%     end
% end

% counter = 1;
% for j = 1:5
%     for i = 1:n-j
%         storage2(counter, 1)=test_track(i+j,1)-test_track(i,1); 
%         dif = (pixSize) *(test_track(i+j, 2)-test_track(i, 2));
%         storage2(counter, 2)=sum(dif.^2);
%         counter = counter +1;
%     end
% end
% 
% for i = 1:5
%     storage3(i) = mean(storage2((storage2(:,1) == i),2));
% end


subtrack = {};
counter = 1;
counter2 =1;
for k = 1:n-5
    for m = 1:5
        subtrack{k,1} = test_track(k:k+m,:);
    end
end


for lmnop = 1:length(subtrack)
    counter = 1;
    storage2 = zeros(15, 2);
for j = 1:5
    for i = 1:6-j
        storage2(counter, 1)=subtrack{lmnop, 1}(i+j,1)-subtrack{lmnop, 1}(i,1); 
        dif = (pixSize) *(subtrack{lmnop, 1}(i+j, 2)-subtrack{lmnop, 1}(i, 2));
        storage2(counter, 2)=sum(dif.^2);
        counter = counter +1;
    end
end
subtrack{lmnop, 2} = storage2;
end
for lmnop = 1:length(subtrack)
for i = 1:5
    subtrack{lmnop, 3}(i) = mean(subtrack{lmnop, 2}((subtrack{lmnop, 2}(:,1) == i),2));
end
end
hold on
%%
figure(1)
hold on
for i = 1:length(subtrack)
    plot([1,2,3,4,5]*test_time,subtrack{i,3})
    plot([1,2]*test_time, subtrack{i,3}(1:2), "LineWidth", 2)
    holdingplace(i) = ((subtrack{i,3}(2)-subtrack{i,3}(1))*pixSize)/(1*(test_time/1000));
    holdingplace2(i) = ((subtrack{i,3}(3)-subtrack{i,3}(1))*pixSize)/(2*(test_time/1000));
    holdingplace3(i) = ((subtrack{i,3}(4)-subtrack{i,3}(1))*pixSize)/(3*(test_time/1000));

%     p = polyfit([1,2],subtrack{i,3}(1:2),1);
%     f = polyval(p, [1,2]);
%     plot([1,2],f, "LineWidth", 2)
end
hold off
figure(2)
hold on
histogram(holdingplace/2)
histogram(holdingplace2/2)
histogram(holdingplace3/2)

hold off
%%
function MSDs = my_instantaneous_DifCoef(mydata)
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

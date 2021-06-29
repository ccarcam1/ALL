%% plot both kymo and cropped area for presentations use mydata + kymos
for i = [70:100]
    %reset before the next figure comes up
    close all
    figure;
    %Change pop-up figure window size to fit screen
    set(gcf,'Position',[375,185,600,800]);% for laptop
    %Display the original kymograph in grayscale
    imagesc(mydata(i).full_kymo,[0,max(mydata(i).crop, [], 'all')/2]);
    colormap(gray);
    hold on
% potentially unneeded
%     ksize = size(mydata(i).full_kymo);
%     kw=ksize(:,2);
%     kh=ksize(:,1);
    clear crop_coordinates
    crop_coordinates = mydata(i).crop_coords;
    % box off the segmented area on the full kymograph
    rectangle('Position',[crop_coordinates(3),crop_coordinates(1),...
        crop_coordinates(4)-crop_coordinates(3),crop_coordinates(2)-...
        crop_coordinates(1)], 'EdgeColor','y'); 
    title([mydata(i).mol_id, newline])
    %x axis dimensions
    xticklabels = [(0:1:(size(mydata(i).full_kymo,2)*0.1))];
    xticks = linspace(1, size(mydata(i).full_kymo, 2), numel(xticklabels));
    set(gca, 'XTick', xticks, 'XColor','w', 'XTickLabel', xticklabels,'FontSize',16)
    %y axis dimensions
    zz = (size(mydata(i).full_kymo,1)*mydata(i).timestep)/1000;
    yy = round(zz/10);
    yticklabels = ([0:yy:zz]);
    yticks = linspace(1, size(mydata(i).full_kymo,1), numel(yticklabels));
    set(gca, 'YTick', yticks, 'YColor','w', 'YTickLabel', yticklabels,'YDir','reverse','FontSize',16)
    ylabel('time (s)')
    xlabel('distance (\mum)')
    hold off
%%%
%%% 
    figure;
    crop = mydata(i).crop;
    %Change pop-up figure window size to fit screen
    set(gcf,'Position',[1025,185,200,800]); % for laptop
    imagesc(crop,[0,max(crop, [], 'all')/2]),colormap(gray);
    hold on
    plot(mydata(i).coords,mydata(i).frames,'-y', 'LineWidth', 2.2);
    %x axis dimensions    
    xticklabels = [(0:0.5:(size(crop,2)*0.1))];
    xticks = linspace(1, size(crop, 2), numel(xticklabels));
    set(gca, 'XTick', xticks,'XColor','w', 'XTickLabel', xticklabels, 'FontSize',16)
    %y axis dimensions
    zz = (size(crop,1)*mydata(i).timestep)/1000;
    yy = round(zz/10);
    yticklabels = ([0:yy:zz]);
    yticks = linspace(1, size(crop,1), numel(yticklabels));
    set(gca, 'YTick', yticks,'YColor','w', 'YTickLabel', yticklabels,'YDir','reverse','FontSize',16)
    ylabel('time (s)')
    xlabel('distance (\mum)')
    hold off
    % interactive element for switching to the next kymograph segment
    choice = menu('Choose an action','next','Break');
    if choice == 2
        break 
    end
end



%%

clf
for i = 1:length(tracks_test)
    plot(tracks_test{i,1}(:,1),tracks_test{i,1}(:,2));
    hold on
end

hold off
% xlim([0 10])
xlabel('time (sec)');
ylabel('distance (um)');
titlestr = 'cas9';
title(titlestr);
% saveas (gcf, strcat(titlestr, '.dpng'), 'png')



%% Plot cropped area plus trajectory
for k = 1:3
fig = figure;
imagesc(mydata(k).crop);
colormap(gray);
hold on
plot(mydata(k).coords,mydata(k).frames,'-y');
% title('overlay');
% xlabel('x, px');
% ylabel('y, px');
axis off;
hold off
iptsetpref('ImshowBorder','tight')
print(fig, ['test', num2str(k)], '-dpng', '-r1200')

end

%% order of operations (start run here) name mydata

for i = 1:length(mydata)
    mydata(i).length = size(mydata(i).crop,1);
end

for i = 1:length(mydata)
    mydata(i).width = size(mydata(i).crop,2);
end

for i = 1:length(mydata)
    lengthoflength(i) = mydata(i).length;
end

for i = 1:length(mydata)
    lengthofwidth(i) = mydata(i).width;
end

[len2tracenum,I] = sort(lengthoflength,'descend');

for i = 1:length(mydata)
    len2tracenum(2,i) = mydata(I(i)).width;
end

% sumwidth = 0;
% 
% fig = figure;
% for k = 1:3
% imagesc(mydata(k).crop);
% colormap(gray);
% hold on
% plot(mydata(k).coords,mydata(k).frames,'-y');
% % title('overlay');
% % xlabel('x, px');
% % ylabel('y, px');
% axis off;
% hold off
% iptsetpref('ImshowBorder','tight')
% print(fig, ['test', num2str(k)], '-dpng', '-r1200')
% sumwidth = sumwidth + mydata(k).width;
% end


%% Plot tracks of similar lengths next to one another

middleof = abs(len2tracenum(1,:)- mean(len2tracenum(1,:)));
find(middleof == min(middleof))

% minimumstart = 1 %min(find(B(1,:)<(mean(B(1,:)))+20))
% maxstop = 19%max(find(B(1,:)>(mean(B(1,:)))-20))
minimumstart = find(middleof == min(middleof))-4;
maxstop = find(middleof == min(middleof))  +15;
% new2(1:2,:) = B(1:2,new)
% new2(3,:) = I(new)
startingpt = minimumstart;
cutoffpoint = maxstop;
figure_manytraces = zeros(mydata(I(startingpt)).length,sum(len2tracenum(2,startingpt:cutoffpoint)));
 
widthtonow = 1;
for i = startingpt:cutoffpoint
    figure_manytraces(1:(len2tracenum(1,i)),(widthtonow:(len2tracenum(2,i)+widthtonow-1))) = mydata(I(i)).crop;
    widthtonow = widthtonow + len2tracenum(2,i);
end

fig=figure;
imagesc(figure_manytraces,[0,max(figure_manytraces, [], 'all')/3])
colormap(gray);
axis off;
iptsetpref('ImshowBorder','tight')
% print(fig, ['test for cutoffpoint equal ' , num2str(startingpt), ' to ', num2str(cutoffpoint)], '-dpng', '-r1200')


%% determine window to extract traces from 
% input a number which is the length of the trace you want all your traces
% to be of similar length to. Pick the trace length that seems most
% frequent based on the histogram of trace lengths

histogram(lengthoflength,50);
prompt = "what is the most frequent length?";
mostfreqlen = input(prompt);
[~,nearestone] = (min(abs(len2tracenum(1,:) - mostfreqlen)));

startingpt = nearestone -10;
cutoffpoint = nearestone +9;
figure_manytraces = zeros(mydata(I(startingpt)).length,sum(len2tracenum(2,startingpt:cutoffpoint)));

widthtonow = 1;
for i = startingpt:cutoffpoint
    figure_manytraces(1:(len2tracenum(1,i)),(widthtonow:(len2tracenum(2,i)+widthtonow-1))) = mydata(I(i)).crop;
    widthtonow = widthtonow + len2tracenum(2,i);
end

figure(3)
for i = 1:size(figure_manytraces, 1)
    for j = 1:size(figure_manytraces, 2)
        if figure_manytraces(i,j) > 1 && figure_manytraces(i,j) <5
            figure_manytraces(i,j) = figure_manytraces(i,j)*2;
        else
            figure_manytraces(i,j) = figure_manytraces(i,j)*1;
        end
    end
end

imagesc(figure_manytraces,[0,max(figure_manytraces, [], 'all')/3])
% colormap(parula);
colormap(flipud(gray));
axis off;
iptsetpref('ImshowBorder','tight')



%% Track summary plus scalebars
middleof = abs(len2tracenum(1,:)- mean(len2tracenum(1,:)));
find(middleof == min(middleof))
minimumstart = find(middleof == min(middleof))-4;
maxstop = find(middleof == min(middleof))  +15;
% new2(1:2,:) = B(1:2,new)
% new2(3,:) = I(new)
startingpt = minimumstart;
cutoffpoint = maxstop;
figure_manytraces = zeros(mydata(I(startingpt)).length,sum(len2tracenum(2,startingpt:cutoffpoint)));
widthtonow = 1;
for i = startingpt:cutoffpoint
    figure_manytraces(1:(len2tracenum(1,i)),(widthtonow:(len2tracenum(2,i)+widthtonow-1))) = mydata(I(i)).crop;
    widthtonow = widthtonow + len2tracenum(2,i);
end

fig=figure;
hold on
 %x axis dimensions    
    xticklabels = [(0:3:(size(figure_manytraces,2)*0.1))];
    xticks = linspace(1, size(figure_manytraces, 2), numel(xticklabels));
    set(gca, 'XTick', xticks,'XColor','k', 'XTickLabel', xticklabels, 'FontSize',12)
%y axis dimensions
    zz = (size(figure_manytraces,1)*mydata(i).timestep)/1000;
    yy = round(zz/10);
    yticklabels = ([0:yy:zz]);
    yticks = linspace(1, size(figure_manytraces,1), numel(yticklabels));
    set(gca, 'YTick', yticks,'YColor','k', 'YTickLabel', yticklabels,'YDir','reverse','FontSize',12)
    ylabel('time (s)')
    xlabel('distance (\mum)')

imagesc(figure_manytraces,[0,max(figure_manytraces, [], 'all')/3])
   
hold off
colormap(gray);






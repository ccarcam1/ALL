% green map
green_map = zeros(3);
for i = 1:11
green_map(i,2) = i*0.1-0.1;
end
% red map
red_map = zeros(3);
for i = 1:11
red_map(i,1) = i*0.1-0.1;
end
%%
fig = figure();
imagesc(data(20).green_kymo)
colormap(fig,green_map)
caxis([0 16]) 
hold on
axis off
fig = figure();
imagesc(data(20).red_kymo)
colormap(fig,red_map)
caxis([0 8])
hold off
axis off

%%
for i = 1:length(data)
%     figure(i)
    RGB = rand(size(data(i).red_kymo,1), size(data(i).red_kymo, 2), 3);
    RGB(:,:,1) = data(i).red_kymo;
    RGB(:,:,2) = data(i).green_kymo;
    RGB(:,:,3) = 0;
%     imagesc(RGB)
    cellhere{i} = RGB;
end
%%

C = imfuse(data(2).green_kymo,data(2).red_kymo, 'Scaling', 'independent');
% C = imresize(C, 25);
imshow(C)

%%
imshowpair(data(2).green_kymo,data(2).red_kymo, 'Scaling', 'independent')

function new_image = my_threshold_function(image,threshold)
[row, col] = find(image>threshold);
for i = 1:length(row)
    image(row(i), col(i))= 0;
end
new_image = image;
end
function gaussfitting = my_gaussian_fitting(x,y,segmented_kymos)

total_time = tic;
parfor i = x:y
    sz=size(segmented_kymos(i).crop);
    t=sz(:,1);
    gaussfitting(i).name = segmented_kymos(i).name;
    gaussfitting(i).data = my_internal_fitting_function(segmented_kymos(i), t);
disp(['you are on trace ', num2str(i), ' out of ', num2str(y)])
end
toc(total_time)
end

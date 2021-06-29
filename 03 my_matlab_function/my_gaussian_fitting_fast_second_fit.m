function gaussfitting = my_gaussian_fitting_fast_second_fit(segmented_kymos, fit_again_later, masked_crops)
total_time = tic;
for i = 1:length(fit_again_later)
    j = fit_again_later(i);
    sz=size(segmented_kymos(j).crop);
    t=sz(:,1);
    gaussfitting(i).name = segmented_kymos(j).name;
    gaussfitting(i).data = my_fast_particle_finding_plus_masked(segmented_kymos(j), t, 10, masked_crops{j});
    gaussfitting(i).original_position = j;
disp(['you are on trace ', num2str(i), ' out of ', num2str(length(fit_again_later))])
end
disp('DONE!')
toc(total_time)
end
%%
% subplot(2,1,1)
% figure(1)
i = 3;
scatter(x_SEM,slope)
title('SEM weights vs unweighted initial')
xlabel('slope using SEM weights')
ylabel('slope using unweighted initial')
% subplot(2,1,2)
print('SEM weights vs unweighted initial intercepts through zero', '-dpng')
% figure(2)
x = [0 max(MSDss(i).MSD(:,1))/4];
% y = slope(i)*x+intercept(i);
y = slope(i)*x;
plot(x,y)
hold on
y = x_SEM(i)*x;
plot(x,y)
hold on
errorbar(MSDss(i).MSD(:,1), MSDss(i).MSD(:,2),MSDs(i).MSD(:,3))
hold off
legend('initial linear portion fit (unweighted)', 'weighted by SEM')
title('SEM weights vs unweighted initial')
xlabel('time lag(s)')
ylabel('MSD(\mum^{2})')
ylim([0 max(MSDss(i).MSD(:,2))+max(MSDs(i).MSD(:,3))])
print('SEM weights vs unweighted initial MSD with slopes intercepts through zero', '-dpng')
%%
% subplot(2,1,1)
% figure(1)
i = 3;
scatter(slope_full,slope)
title('unweighted full vs unweighted initial')
xlabel('slope using unweighted full')
ylabel('slope using unweighted initial')
% subplot(2,1,2)
print('unweighted full vs unweighted initial intercepts through zero', '-dpng')
% figure(2)
x = [0 max(MSDss(i).MSD(:,1))/4];
% y = slope(i)*x+intercept(i);
y = slope(i)*x;
plot(x,y)
hold on
y = slope_full(i)*x;
plot(x,y)
hold on
scatter(MSDss(i).MSD(:,1), MSDss(i).MSD(:,2))
hold off
legend('initial linear portion fit (unweighted)', 'full fit (unweighted)')
title('unweighted full vs unweighted initial')
xlabel('time lag(s)')
ylabel('MSD(\mum^{2})')
ylim([0 max(MSDss(i).MSD(:,2))+max(MSDs(i).MSD(:,3))])
print('unweighted full vs unweighted initial MSD with slopes intercepts through zero', '-dpng')

%%
i = 3;
scatter(x_sampleSDev,slope)
title('Sample SDev vs unweighted initial')
xlabel('slope using SampleSDev weights')
ylabel('slope using unweighted initial')
% subplot(2,1,2)
print('Sample SDev vs unweighted initial intercepts through zero', '-dpng')
% figure(2)
x = [0 max(MSDss(i).MSD(:,1))/4];
% y = slope(i)*x+intercept(i);
y = slope(i)*x;
plot(x,y)
hold on
y = x_sampleSDev(i)*x;
plot(x,y)
hold on
errorbar(MSDss(i).MSD(:,1), MSDss(i).MSD(:,2),MSDss(i).sampleSDev)
hold off
legend('initial linear portion fit (unweighted)', 'weighted by Sample SDev')
title('Sample SDev vs unweighted initial')
xlabel('time lag(s)')
ylabel('MSD(\mum^{2})')
ylim([0 max(MSDss(i).MSD(:,2))+max(MSDs(i).MSD(:,3))])
print('Sample SDev vs unweighted initial MSD with slopes intercepts through zero', '-dpng')

%%
i = 3;
scatter(x_MichaeletError,slope)
title('Michaelet weights vs unweighted initial')
xlabel('slope using Michaelet weights')
ylabel('slope using unweighted initial')
% subplot(2,1,2)
print('Michaelet weights vs unweighted initial intercepts through zero', '-dpng')
% figure(2)
x = [0 max(MSDss(i).MSD(:,1))/4];
% y = slope(i)*x+intercept(i);
y = slope(i)*x;
plot(x,y)
hold on
y = x_MichaeletError(i)*x;
plot(x,y)
hold on
errorbar(MSDss(i).MSD(:,1), MSDss(i).MSD(:,2),MSDss(i).MSD_MichaeletError)
hold off
legend('initial linear portion fit (unweighted)', 'Michaelet weights')
title('Michaelet weights vs unweighted initial')
xlabel('time lag(s)')
ylabel('MSD(\mum^{2})')
ylim([0 max(MSDss(i).MSD(:,2))+max(MSDs(i).MSD(:,3))])
print('Michaelet weights vs unweighted initial MSD with slopes intercepts through zero', '-dpng')
%%
x = [0 max(MSDss(i).MSD(:,1))/4];
y = slope(i)*x;
plot(x,y)
hold on
% y = slope_full(i)*x;
% plot(x,y)
% hold on
% y = x_SEM(i)*x;
% plot(x,y)
% hold on
% y = x_sampleSDev(i)*x;
% plot(x,y)
% hold on
y = x_MichaeletError(i)*x;
plot(x,y)
hold on
y = B(i)*x+A(i);
plot(x,y)
legend('initial linear portion fit (unweighted)', 'Michaelet weights', 'Michaelet weights my fit function')
% legend('initial linear portion fit (unweighted)','full fit (unweighted)', 'weighted by SEM', 'weighted by Sample SDev', 'Michaelet weights')
scatter(MSDss(i).MSD(:,1), MSDss(i).MSD(:,2))
title('Weighted fits using Matlab lscov function')
xlabel('time(s)')
ylabel('MSD(\mum^{2})')
hold off
% print('Weighted fits using Matlab lscov function', '-dpng')
%%

i = 174;
% % scatter(Slope_lscov,B)
% % title('lscovfunction vs myfunction')
% % xlabel('lscovfunction weighted slope')
% % ylabel('myfunction weighted slope')
% % % subplot(2,1,2)
% % print([num2str(i),'lscovfunction vs myfunction'], '-dpng')
% figure(2)
x = [0 max(MSDs(i).MSD(:,1))/4];
% y = slope(i)*x+intercept(i);
y = B(i)*x+A(i);
plot(x,y)
hold on
y = Slope_lscov(i)*x;
plot(x,y)
hold on
errorbar(MSDs(i).MSD(:,1), MSDs(i).MSD(:,2),MSDs(i).MSD_MichaeletError)
hold off
legend('myfunction', 'lscovfunction')
title('myfunction vs lscovfunction')
xlabel('time lag(s)')
ylabel('MSD(\mum^{2})')
% ylim([0 max(MSDss(i).MSD(:,2))+max(MSDs(i).MSD(:,3))])
print([num2str(i),'lscovfunction vs myfunction MSD with slopes intercepts through zero'], '-dpng')
%%
% D_myfun = B/2;
% D_lscov = Slope_lscov/2;
figure(3)
histogram(D_lscov,20)%,'BinLimits',[-0.01,0.14])
hold on
histogram(D_myfun,20)%,'BinLimits',[-0.01,0.14])
legend('forced 0', 'not forced')
hold off
xlabel('Diffusion Coefficient (\mum^{2}/s)')
ylabel('count')

% figure(2)
% histogram(A,20,'BinLimits',[-0.3,0.4])
% hold on
% histogram(epsilon_val,20,'BinLimits',[-0.3,0.4])
% xlabel('Y intercept (\mum^{2}/s)')
% ylabel('count')
% legend('forced 0', 'not forced')
% hold off
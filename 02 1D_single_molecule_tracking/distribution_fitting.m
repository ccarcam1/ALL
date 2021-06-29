% lognormal distribution play 
cd('D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\20191203\');

clf
figure(1)
set(gcf,'Position',[60,30,500,1100])
subplot(3,1,1)
dist = Dsalt(saltsalt ==200);
h = histogram(dist,30);
binedges = h.BinEdges;
binwidthhalf = h.BinWidth/2;
valuesare = h.Values;
x = binedges(1:end-1)+binwidthhalf;
logx = log(x);
y = valuesare;
hold on
pd = fitdist(dist','Lognormal')
peakis = exp(pd.mu);
xline(peakis,'--r','DisplayName','lognormal','LineWidth',2);
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
pd = fitdist(dist','Normal')
peakis = pd.mu;
xline(peakis,'--b','DisplayName','normal','LineWidth',2);
peakis = median(dist);
xline(peakis,'--g','DisplayName','median','LineWidth',2);
ylabel('counts')
title('200mM KCl + ATP')
legend('show')
hold off

subplot(3,1,2)
plot(x,y)
xlabel('x = um^2/sec')
ylabel('y (counts)')
title('bin values')

subplot(3,1,3)
plot(logx,y)
hold on
pd = fitdist(dist','Lognormal')
peakis = pd.mu;
xline(peakis,'--r','DisplayName','lognormal','LineWidth',2);
hold off
xlabel('log(x)')
ylabel('y (counts)')
title('bin values')


%%

% pd = fitdist(dist','Lognormal')
% x_values = 50:1:250;
% y = pdf(pd,x_values);
% plot(x_values,y,'LineWidth',2)


exp(pd.mu)
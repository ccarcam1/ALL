
x = 1;
for i = 1:length(saltdata)
    if saltdata(i).r2point8 == 1
        Dsalt(x) = saltdata(i).slope/2;
        saltsalt(x) = saltdata(i).salt;
        x = x+1;
    else
    end
end

%%
hold on
x = DnoATP(saltnoATP ==70);
% x = Dsalt(saltsalt ==70);

mean(x)
histogram(x,30)
hold on
pd = fitdist(x', 'Lognormal');
x_values = 0:.001:.14;
y = pdf(pd,x_values);
plot(x_values,y,'LineWidth',2)
%%


% boxplot(list,name)

boxplot(D,salt)
title(['SWR1 1mM ATP Diffusion at', newline, 'different monovalent salt concentrations'])
xlabel('mM KCl')
ylabel('Diffusion Coefficients (\mum^2/sec)')
%% 

cd('D:\OneDrive - Johns Hopkins University\Ha_CCarcamo\Data\SWR1\20191203\');

clf
figure(1)
dist = Dsalt(saltsalt ==25);
D25mean = mean(Dsalt(saltsalt ==25));
D25median = median(Dsalt(saltsalt ==25));
D25SEM = std(Dsalt(saltsalt ==25))/sqrt(length(Dsalt(saltsalt ==25)));
D25STD = std(Dsalt(saltsalt ==25));
pd = fitdist(dist','Lognormal')
D25lognorm = exp(pd.mu);

dist = Dsalt(saltsalt ==70);
D70mean = mean(Dsalt(saltsalt ==70));
D70median = median(Dsalt(saltsalt ==70));
D70SEM = std(Dsalt(saltsalt ==70))/sqrt(length(Dsalt(saltsalt ==70)));
D70STD = std(Dsalt(saltsalt ==70));
pd = fitdist(dist','Lognormal')
D70lognorm = exp(pd.mu);

dist = Dsalt(saltsalt ==200);
D200mean = mean(Dsalt(saltsalt ==200));
D200median = median(Dsalt(saltsalt ==200));
D200SEM = std(Dsalt(saltsalt ==200))/sqrt(length(Dsalt(saltsalt ==200))); 
D200STD = std(Dsalt(saltsalt ==200)); 
pd = fitdist(dist','Lognormal')
D200lognorm = exp(pd.mu);

dist = DnoATP(saltnoATP ==70);
noATPD70mean = mean(DnoATP(saltnoATP ==70));
noATPD70median = median(DnoATP(saltnoATP ==70));
noATPD70SEM = std(DnoATP(saltnoATP ==70))/sqrt(length(DnoATP(saltnoATP ==70)));
noATPD70STD = std(DnoATP(saltnoATP ==70));
pd = fitdist(dist','Lognormal')
noATPD70lognorm = exp(pd.mu);

dist = DgammaS(saltgammaS ==70);
gammaSD70mean = mean(DgammaS(saltgammaS ==70));
gammaSD70median = median(DgammaS(saltgammaS ==70));
gammaSD70SEM = std(DgammaS(saltgammaS ==70))/sqrt(length(DgammaS(saltgammaS ==70)));
gammaSD70STD = std(DgammaS(saltgammaS ==70));
pd = fitdist(dist','Lognormal')
gammaSD70lognorm = exp(pd.mu);

x = [25 70 200];

y = [D25mean D70mean D200mean];
yy = [D25median D70median D200median];
yyy = [D25lognorm D70lognorm D200lognorm];
std = [D25STD D70STD D200STD];
SEM = [D25SEM D70SEM D200SEM];



% e = errorbar(x, y, SEM, 'DisplayName','ATP','LineWidth',1);
% e.Marker = 'o';
% e.LineStyle = '-';
% hold on

% e = errorbar(x, yy+(y(1)-yy(1)), SEM);
e = errorbar(x, yy, SEM, 'DisplayName','ATP','LineWidth',1);

e.Marker = 'o';
e.LineStyle = '-';

% % e = errorbar(x, yyy+(y(1)-yyy(1)), SEM);
% e = errorbar(x, yyy, SEM, 'DisplayName','ATP','LineWidth',1);
% e.Marker = 'o';
% e.LineStyle = '-';

% hold off
hold on
x = [70];
y = [noATPD70mean];
yy = [noATPD70median];
yyy = [noATPD70lognorm];
std = [noATPD70STD];
SEM = [noATPD70SEM];
errorbar(x, yy, SEM, 'o','DisplayName','noATP','LineWidth',1)
e.Marker = 'o';
e.LineStyle = '-';

x = [70];
y = [gammaSD70mean];
yy = [gammaSD70median];
yyy = [gammaSD70lognorm];
std = [gammaSD70STD];
SEM = [gammaSD70SEM];
errorbar(x, yy, SEM, 'o','DisplayName','gammaS','LineWidth',1)
e.Marker = 'o';
e.LineStyle = '-';

axis([15 220 0.010 0.055])
% title(['SWR1 1mM ATP Diffusion at', newline, 'different monovalent salt concentrations'])
title(['SWR1 1mM ATP Diffusion all conditions'])

xlabel('mM KCl','FontSize',14)
ylabel(['median Diffusion Coefficients', newline,'(\mum^2/sec) with SEM'],'FontSize',14)
legend('show','Location','northwest')

clear std
clear x
clear y
clear D200mean
clear D200STD
clear D25mean
clear D25STD
clear D70mean
clear D70STD
clear noATPD70mean
clear noATPD70STD
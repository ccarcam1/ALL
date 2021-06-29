x = 1;
for i =5:7
leftmost(x) = mydata(i).crop_coords(3);
peakposition = mydata(i).coords(:,1);
avgpeakposition(x) = mean(peakposition);
position(x) = leftmost(x) + avgpeakposition(x);
x = x+1;
end
position = position*0.1;
%%

realdistance1 =30520 ;
realdistance2 = 18070;
realdistance3 = 12130;

distance1 = position(1);
distance2 = position(2);
distance3 = position(3);

x = [distance1 distance2 distance3];
y = [realdistance1 realdistance2 realdistance3];


[p,S] = polyfit(x,y,1); 
% f = polyval(p,x); 
% plot(x,y,'o',x,f,'-') 
% legend('data','linear fit')
% 
% 
% % plot(x, y, '+')


[y_fit,delta] = polyval(p,x,S);
plot(x,y,'bo')
hold on
plot(x,y_fit,'r-')
plot(x,y_fit+2*delta,'m--',x,y_fit-2*delta,'m--')

xlabel('um distance')
ylabel('bp value from 1:48502')
title('um position vs real bp value')

%%


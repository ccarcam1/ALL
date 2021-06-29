file = '20191010-204131 Kymograph 63.h5';
% h5disp(file);


% See contents
info = h5info(file);
% info = h5info(file, '/Force HF');
% info = h5info(file, '/Distance');




% h5read(file);
% h5disp(file, '/Distance/Distance 1')
% distance = h5read(file, '/Distance/Distance 1');
% ForceHF = h5read(file, '/Force HF/Force 1x');
% photoncountgreen = h5read(file, '/Photon count/Green');
%%
forcedata(1).forceHF =  h5read(file, '/Force HF/Force 1x');
forcedata(1).samplerate =  h5readatt(file, '/Force HF/Force 1x', 'Sample rate (Hz)');
forcedata(1).Stop_time_ns = h5readatt(file, '/Force HF/Force 1x', 'Stop time (ns)');
forcedata(1).Start_time_ns = h5readatt(file, '/Force HF/Force 1x', 'Start time (ns)');
forcedata(1).distance_time = h5read(file, '/Distance/Distance 1');


forcedata(1).time = forcedata(1).distance_time.Timestamp-forcedata(1).distance_time.Timestamp(1);
forcedata(1).time = double(forcedata(1).time);
forcedata(1).time = forcedata(1).time/1e9;

forcedata(1).tHF = forcedata(1).Stop_time_ns -forcedata(1).Start_time_ns;
forcedata(1).tHF = double(forcedata(1).tHF);
forcedata(1).tHF = forcedata(1).tHF/1e9;
stepsize = (1/(forcedata(1).samplerate));

forcedata(1).timeHF = [0:stepsize:forcedata(1).tHF-stepsize]';

forcedata(1).distance = forcedata(1).distance_time.Value;
forcedata = rmfield(forcedata, 'distance_time');


%% Downsampling
Htzwanted = 100;
bythis = forcedata.samplerate/Htzwanted;
roundedthis = round(bythis);
forcedata.forceLF = downsample(forcedata.forceHF, roundedthis);
forcedata.timeLF = downsample(forcedata.timeHF, roundedthis);

% distanceDSnumb = length(forcedata.forceHF)/length(forcedata.distance);
% forcedata.forcevLF = downsample(forcedata.forceHF, distanceDSnumb); 

%%
% % whenhappensit = zeros(length(forcedata.time), 1);
% % for i = 1:length(forcedata.time)
% %     whenhappensit(i)= find(forcedata.timeHF == forcedata.time(i));
% % end
% % Use the information from the previous loop to determine the step size for
% % the parsing of the time and force data to match the distance data
for i = 1:length(forcedata.time)
    downsampleontime(i) = 5208*i +1 -5208;
end

forcedata.forcevLF = forcedata.forceHF(downsampleontime);
forcedata.timevLF = forcedata.timeHF(downsampleontime);



%% Compare the downsampling success: Force vs time FD curve
% % figure(1)

clear ind1
clear ind2
clear ind3
clear ind4

BLUE = '#0072BD';
LIGHTGREY = '#B8B8B8';
DARKGREY = '#474746';
GREEN = '#77AC30';
BLACK = 'k';
YELLOW= '#EDB120';

subplot(2, 1, 1);
plot(forcedata.timeHF, forcedata.forceHF, 'Color', LIGHTGREY)
hold on
plot(forcedata.timeLF, forcedata.forceLF,'Color' , DARKGREY)
hold on

[x,y] = ginput
ind1 = interp1(forcedata.timevLF,1:length(forcedata.timevLF),x(1),'nearest');
ind2 = interp1(forcedata.timevLF,1:length(forcedata.timevLF),x(2),'nearest');
ind3 = interp1(forcedata.timevLF,1:length(forcedata.timevLF),x(3),'nearest');
ind4 = interp1(forcedata.timevLF,1:length(forcedata.timevLF),x(4),'nearest');

% endis = length(forcedata.distance);
% incrementis = [600:endis-100];

plot(forcedata.time(ind1:ind2), forcedata.forcevLF(ind1:ind2),'LineWidth',...
    1, 'Color' , BLUE )% BLUE
plot(forcedata.time(ind3:ind4), forcedata.forcevLF(ind3:ind4),'LineWidth',...
    1, 'Color', GREEN )% GREEN

hold off


xlabel('time (s)');
ylabel('Force (pN)');


% %% Compare the downsampling success: Force vs distance
% BLUE = '#0072BD';
% GREEN = '#77AC30';
% BLACK = 'k';
% YELLOW= '#EDB120';

% % figure(2)
% endis = length(forcedata.distance);
% incrementis = [600:endis-100];
subplot(2, 1, 2);
plot(forcedata.distance(ind1:ind2), forcedata.forcevLF(ind1:ind2), 'LineWidth', 1,'Color' , BLUE )
hold on
plot(forcedata.distance(ind3:ind4), forcedata.forcevLF(ind3:ind4), 'LineWidth', 1,'Color', GREEN)
hold off
xlabel('Distance (um)');
ylabel('Force (pN)');




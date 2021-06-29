function bool_plot(i, data, bool)
    subplot(3,1,1)
    plot(data(i).distance(bool), data(i).force(bool))
    xlabel('Distance')
    ylabel('Force')
    subplot(3,1,2)
    plot(data(i).time(bool), data(i).distance(bool))
    xlabel('Time')
    ylabel('Distance')
    subplot(3,1,3)
    plot(data(i).time(bool), data(i).force(bool))
    xlabel('Time')
    ylabel('Force')
end
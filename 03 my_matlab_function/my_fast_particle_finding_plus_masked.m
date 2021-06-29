function save_these_results = my_fast_particle_finding_plus_masked(trace,t, max_distance_for_connect, masked_crop)
fittingcell = {};
gof_results = {};
datacell = {};
height = [];
width = [];
pixel = [];
timepix = [];
intensity = [];
for j=[1:t]
        subcrop=masked_crop(j,:);
%         subcrop=trace.crop(j,:);
        transpose(subcrop);  
        fitresult = cell( 2, 1 );
        gof = struct( 'sse', cell( 2, 1 ), ...
        'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
        [xData, yData] = prepareCurveData(trace.wrange, subcrop);
        ft = fittype( 'd + (a*exp(-((x-b)/sqrt(2)*c)^2))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.Lower = [trace.amin trace.bmin -2 trace.dmin];
        opts.StartPoint = [trace.a trace.b -1.3 trace.d];
        opts.Upper = [trace.amax trace.bmax -1 trace.dmax];
        [fitresult, gof] = fit( xData, yData, ft, opts );
        warning('off','signal:findpeaks:largeMinPeakHeight')
        [vals, inds] = findpeaks(yData,xData,'SortStr','descend','MinPeakHeight',2);
%         [msgStr,warnId] = lastwarn
        if gof.rsquare < 0.05 && ~isempty(inds) && not(abs(inds(1)- trace.b) > max_distance_for_connect)
            trace.b = inds(1);
            trace.bmin=trace.b-5;
            trace.bmax=trace.b+5;
            ft = fittype( 'd + (a*exp(-((x-b)/sqrt(2)*c)^2))', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.Lower = [trace.amin trace.bmin -2 trace.dmin];
            opts.StartPoint = [trace.a trace.b -1.3 trace.d];
            opts.Upper = [trace.amax trace.bmax -1 trace.dmax];
            [fitresult, gof] = fit( xData, yData, ft, opts );
        end
        fitcoef=coeffvalues(fitresult);
        height(j)=fitcoef(:,1);
        width(j)=abs(2*fitcoef(:,3));
        pixel(j)=fitcoef(:,2);
        timepix(j)=j;
        intensity(j)= height(j)*width(j); %fitcoef(:,1)*2*(fitcoef(:,3);
        trace.a=fitcoef(:,1);
        dsd=(trace.d - trace.dmin)/2;
        trace.amin=trace.a-dsd;
        trace.amax=trace.a+dsd;
        trace.b=fitcoef(:,2);
        trace.bmin=trace.b-5;
        trace.bmax=trace.b+5;
        fittingcell{j} = fitresult;
        gof_results{j} = gof;
        datacell{j} = [xData, yData];  
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%temp section%%%%%%%%%%%%%%%%%%%
%         findpeaks(yData,xData)
%         hold on
%         plot(fitresult, xData, yData)
%         hold off
%         title(num2str(gof.rsquare))
%         ylabel(strcat("y slice = ", num2str(j)))
%         xlabel('x')
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end of temp section%%%%%%%%%%%
% disp(j)
end
disp([num2str(j), '/', num2str(t)])
save_these_results.fitting = fittingcell;
save_these_results.gof = gof_results;
save_these_results.data = datacell;
save_these_results.height = height;
save_these_results.width = width;
save_these_results.pixel = pixel;
save_these_results.timepix = timepix;
save_these_results.intensity = intensity;
end



function save_these_results = my_internal_fitting_function(trace,t)
fittingcell = {};
gof_results = {};
datacell = {};
height = [];
width = [];
pixel = [];
timepix = [];
intensity = [];
for j=[1:t]
        subcrop=trace.crop(j,:);
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
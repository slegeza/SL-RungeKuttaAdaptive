function y0 = slrk5(t0, tMax, y0, h)
%t0 = initial time, tMax = final time, y0= initial condition, h = initial
%timestep


func = @(t, y) exp(t)*sin(y); %problem 1
%func = @(t,v) 9.8 - (10e-4/10e-2)*v.^2; %problem 2
%func = @(t,x) 1/(t.^2); %problem 3

tn =t0;
tArr = [];
yArr = [];
tArr(1) = t0;
yArr(1) = y0;

%hmin = 1e-9;
hmin=0.001;
nSteps = 0;
while tn <= tMax
    

    
    f0 = func(tn, y0);
    f1 = func(tn+h/4, y0+(h/4)*f0);
    f2 = func(tn+(3*h/8), y0+(3*h/32)*f0 + (9*h/32)*f1);
    f3 = func(tn+(12*h/13), y0+(1932*h/2197)*f0 - (7200*h/2197)*f1 + (7296*h/2197)*f2);
    f4 = func(tn+h, y0+(439*h/216)*f0 - 8*h*f1 + (3680*h/513)*f2 - (845*h / 4104)*f3);
    f5 = func(tn+h/2, y0-(8*h/27)*f0 + 2*h*f1 - (3544*h/2565)*f2 + (1859*h/4104)*f3 - (11*h/40)*f4);
    
    %y(i) = y0 + (h/6)*(f0 + 2*f1 + 2 *f2 + f3);
    %yPrime = y0 + h*((25/216)*f0+(1408/2565)*f2 + (2197/4104)*f3 - (1/5)*f4);
    yHat = y0 + h*((16/135)*f0 + (6656/12825)*f2 + (28561/56430)*f3 - (9/50)*f4 + (2/55)*f5);
    
    err = h*((f0/360) - (128/4275)*f2 - (2197/75240)*f3 + f4/50 + (2/55)*f5 );
    eps = 10e-5;
    %err = abs(y0 - yHat);
    hNew = 0.9 * h * (eps / ((abs(err) / yHat))).^(1.0/5.0)
    
    if hNew >= h
        tn = tn + h
        nSteps = nSteps + 1;
        h = hNew;
        y0 = yHat

    else
        h = 0.9 * hNew;
    end
    if(t0 + h) > tMax
        h = tMax - t0;
    end
    if((tMax - t0)/tMax) < 1.e-10
        error('timestep extends beyond range');
    end
    if hNew < hmin
        error('timestep too small!');
    end
    tArr = [tArr tn];
    yArr = [yArr y0];

end

plot(tArr, yArr, 'Linewidth', 2)
hold on
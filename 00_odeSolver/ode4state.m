function yout = ode4state(F,t0,h,tfinal,y0)
% ODE4  Classical Runge-Kutta ODE solver.
%   yout = ODE4(F,t0,h,tfinal,y0) uses the classical
%   Runge-Kutta method with fixed step size h on the interval
%      t0 <= t <= tfinal
%   to solve
%      dy/dt = F(t,y)
%   with y(t0) = y0.

%   Copyright 2014 - 2015 The MathWorks, Inc.

global stateout a0vec Qvec fix_p sch_p
a0p = zeros(fix_p.ns,1);
qoutp = zeros(1,4);

y = y0;
yout = y;
tag = 0;
d = 0;
for t = t0 : h : tfinal-h
    if tag == 0
        try
            s1 = F(t,y,1);
            a0p = [a0p, a0vec];
            qoutp = [qoutp; Qvec'];
            s2 = F(t+h./2, y+h.*s1./2,0);
            s3 = F(t+h./2, y+h.*s2./2,0);
            s4 = F(t+h, y+h.*s3,0);
            y = y + h.*(s1 + 2.*s2 + 2.*s3 + s4)/6;
            yout = [yout; y];
            if t >= t0+d
                timestamp = datestr(now(),'dd_mm_HHMMSS');
                save(['partial_' timestamp],'yout','fix_p','sch_p','t','a0p')
                d = d+0.2;
            end
        catch ME
            warning('Solution diverged...')
            disp(ME)
            tag = 1;
            yout = [yout; zeros(size(y,1),1)];
            timestamp = datestr(now(),'dd_mm_HHMMSS');
            save(['error_' timestamp],'yout','fix_p','sch_p','t','a0p','ME')
        end
    else
        yout = [yout; zeros(size(y,1),1)];
    end
end
stateout.a0p = a0p;
stateout.qoutp = qoutp;

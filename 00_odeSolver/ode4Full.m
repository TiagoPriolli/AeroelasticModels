function yout = ode4Full(F,t0,h,tfinal,y0)
% ODE4  Classical Runge-Kutta ODE solver.
%   yout = ODE4(F,t0,h,tfinal,y0) uses the classical
%   Runge-Kutta method with fixed step size h on the interval
%      t0 <= t <= tfinal
%   to solve
%      dy/dt = F(t,y)
%   with y(t0) = y0.

%   Copyright 2014 - 2015 The MathWorks, Inc.

global a0vec Qvec stateout fix_p_uvlm sch_p_uvlm clvec
n_step = numel(t0:h:tfinal-h)+1;
Nstates = numel(y0);
yout = zeros(n_step,Nstates);
a0 = zeros(n_step,fix_p_uvlm.N);
qout = zeros(n_step,4);
clout = zeros(n_step,fix_p_uvlm.N);
circout = zeros(fix_p_uvlm.M,fix_p_uvlm.N);
Xgridout = zeros(fix_p_uvlm.M+1,fix_p_uvlm.N+1);
Ygridout = zeros(fix_p_uvlm.M+1,fix_p_uvlm.N+1);
Zgridout = zeros(fix_p_uvlm.M+1,fix_p_uvlm.N+1);
Xwakegridout = zeros(200,fix_p_uvlm.N+1);
Ywakegridout = zeros(200,fix_p_uvlm.N+1);
Zwakegridout = zeros(200,fix_p_uvlm.N+1);
NVRING = zeros(n_step,1);

i_step = 1;
y = y0;
yout(i_step,:) = y;
tag = 0;
d = 0;
try
for t = t0 : h : tfinal-h
    i_step = i_step + 1;
    if tag == 0
       
            s1 = F(t,y,1);
            %
            a0(i_step,:) = a0vec';
            qout(i_step,:) = Qvec';
            clout(i_step,:) = clvec';
            circout(:,:) = sch_p_uvlm.circd;
            Xgridout(:,:) = sch_p_uvlm.Xgrid;
            Ygridout(:,:) = sch_p_uvlm.Ygrid;
            Zgridout(:,:) = sch_p_uvlm.Zgrid;
            Xwakegridout(:,:) = sch_p_uvlm.Xwakegrid(1:200,:);
            Ywakegridout(:,:) = sch_p_uvlm.Ywakegrid(1:200,:);
            Zwakegridout(:,:) = sch_p_uvlm.Zwakegrid(1:200,:);
            NVRING(i_step) = sch_p_uvlm.NVRING;
            %
            s2 = F(t+h./2, y+h.*s1./2,0);
            s3 = F(t+h./2, y+h.*s2./2,0);
            s4 = F(t+h, y+h.*s3,0);
            y = y + h.*(s1 + 2.*s2 + 2.*s3 + s4)/6;
            yout(i_step,:) = y';
            if t >= t0+d
                timestamp = datestr(now(),'dd_mm_HHMMSS');
                save(['partial_' timestamp],'yout','a0','qout','fix_p_uvlm','sch_p_uvlm','t')
                fprintf('Saved partial results t=%d\n',t)
                d = d+0.2/8;
            end
      
          
    else
        yout(i_step,:) = zeros(1,size(y,1));
    end
end
catch ME
warning('Solution diverged...')
disp(ME)
save('ErrorLogStates','a0','yout','qout','clout','ME')
end    
stateout.a0 = a0;
stateout.qout = qout;
stateout.circout = circout;
stateout.Xgrid = Xgridout;
stateout.Ygrid = Ygridout;
stateout.Zgrid = Zgridout;
stateout.Xwakegrid = Xwakegridout;
stateout.Ywakegrid = Ywakegridout;
stateout.Zwakegrid = Zwakegridout;
stateout.NVRING = NVRING;

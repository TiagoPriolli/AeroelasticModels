function [cl,cd,cm,LocalState] = LDVM_FlatPlate(alpha, alphadot, h, hdot,...
                                        LocalState, Parameters)

%% Parameters                                  
c = Parameters.c;
u_ref = Parameters.Uref;
lesp_crit = Parameters.lesp_crit;
dt = Parameters.dt;
pvt = Parameters.pvt;

v_core = Parameters.v_core;
expect_vort = Parameters.expect_vort;
%n_aterm = Parameters.n_aterm;
del_dist = Parameters.del_dist;
%iter_max = Parameters.iter_max;

x = Parameters.x;
theta = Parameters.theta;
dtheta = Parameters.dtheta;
%n_div = fix_p.n_div;

%cam = zeros(n_div,1);
%cam_slope = zeros(n_div,1);

bound = zeros(70,3);
gamma = zeros(70,1);
bound_int = zeros(70,3);

%Convert to dimensoinal quantities
dt = dt*c/u_ref;
v_core = v_core*c;
del_dist = del_dist*c;

%% States
ntev = LocalState.ntev;
nlev = LocalState.nlev;
tev = LocalState.tev;
lev = LocalState.lev;
dist_wind = LocalState.dist_wind;
kelv_enf = LocalState.kelv_enf;
aterm_prev = LocalState.aterm_prev;
levflag = LocalState.levflag;


%% Calculate bound vortex positions at this time step

dist_wind = dist_wind + (u_ref*dt);

bound(:,2) = -((c-pvt*c)+((pvt*c-x)*cos(alpha))+dist_wind);
bound(:,3) = h + ((pvt*c-x).*sin(alpha));

%% TEV shed at every time step
ntev = ntev + 1;

if (ntev==1)
    tev(ntev,2)=bound(70,2)+(0.5*u_ref*dt);
    tev(ntev,3)=bound(70,3);
else
    tev(ntev,2)=bound(70,2)+((1./3.)*(tev(ntev-1,2)-bound(70,2)));
    tev(ntev,3)=bound(70,3)+((1./3.)*(tev(ntev-1,3)-bound(70,3)));
end

%% Calculate Distances

% Size ntev x 70
TevBoundDistx = tev(1:ntev,2)*ones(1,70) - ones(ntev,1)*bound(:,2)';
TevBoundDistz = tev(1:ntev,3)*ones(1,70) - ones(ntev,1)*bound(:,3)';

% Size ntev x ntev
tevt = tev';
TevTevDistx = tev(1:ntev,2)*ones(1,ntev) - ones(ntev,1)*tevt(2,1:ntev);
TevTevDistz = tev(1:ntev,3)*ones(1,ntev) - ones(ntev,1)*tevt(3,1:ntev);

if nlev > 0
   % Size nlev+1 x 70 (extra space for possible lev shedding)
   LevBoundDistx = lev(1:nlev+1,2)*ones(1,70) - ones(nlev+1,1)*bound(:,2)';
   LevBoundDistz = lev(1:nlev+1,3)*ones(1,70) - ones(nlev+1,1)*bound(:,3)';
   
   % Size nlev+1 x nlev+1 (extra space for possible lev shedding)
   LevLevDistx = lev(1:nlev+1,2)*ones(1,nlev+1) - ones(nlev+1,1)*lev(1:nlev+1,2)';
   LevLevDistz = lev(1:nlev+1,3)*ones(1,nlev+1) - ones(nlev+1,1)*lev(1:nlev+1,3)';
    
   % Size nlev+1 x ntev+1 (extra space for possible lev shedding)
   LevTevDistx = lev(1:nlev+1,2)*ones(1,ntev) - ones(nlev+1,1)*tev(1:ntev,2)';
   LevTevDistz = lev(1:nlev+1,3)*ones(1,ntev) - ones(nlev+1,1)*tev(1:ntev,3)';
end

%% Find A0 and TEV assuming no LEV is formed

eps=10e-6;                     

tevIter_prev = 0;
tevIter = -0.01;
kelv_prev = 0;
aterm = zeros(45,1);
bound_circ = 0;
downwash = zeros(70,1);
uind = zeros(70,1);
wind = zeros(70,1);

for iter = 2:100
    if iter == 100
        error('1D iteration failed')
    end
    
    tev(ntev,1) = tevIter;
    
    calcDownwash
    
    kelv = kelv_enf + bound_circ + sum(lev(1:nlev,1)) + sum(tev(1:ntev,1));
    
    if (abs(kelv) < eps)
        break
    end
    
    dkelv=(kelv - kelv_prev)/(tevIter-tevIter_prev);
    tevIter_prev = tevIter;
    kelv_prev = kelv;
    
    tevIter = tevIter-kelv/dkelv;
end
    
aterm(3:4)=0;

aterm(3) = sum((((downwash(2:70).*cos(2*theta(2:70)))+...
    (downwash(1:69).*cos(2*theta(1:69))))./2).*dtheta);

aterm(4) = sum((((downwash(2:70).*cos(3*theta(2:70)))+...
    (downwash(1:69).*cos(3*theta(1:69))))./2).*dtheta);

aterm(3) = (2./(u_ref*pi))*aterm(3);
aterm(4) = (2./(u_ref*pi))*aterm(4);

adot = zeros(4,1);
adot(1)=(aterm(1)-aterm_prev(1))/dt;
adot(2)=(aterm(2)-aterm_prev(2))/dt;
adot(3)=(aterm(3)-aterm_prev(3))/dt;
adot(4)=(aterm(4)-aterm_prev(4))/dt;

le_vel_x = u_ref-(alphadot*sin(alpha)*pvt*c)+uind(1);
le_vel_y = -(alphadot*cos(alpha)*pvt*c)-hdot+wind(1);
lesp=aterm(1);

%% If LESP_crit is exceeded

if (abs(lesp)>lesp_crit)
    if (lesp>0)
        lesp_cond=lesp_crit;
    else
        lesp_cond=-lesp_crit;
    end
    
    nlev = nlev+1; %Added Before
    
    if (levflag==0)
        lev(nlev,2)=bound(1,2)+(0.5*le_vel_x*(dt));
        lev(nlev,3)=bound(1,3)+(0.5*le_vel_y*(dt));
    else
        lev(nlev,2)=bound(1,2)+((1./3.)*(lev(nlev-1,2)-bound(1,2)));
        lev(nlev,3)=bound(1,3)+((1./3.)*(lev(nlev-1,3)-bound(1,3)));
    end
    levflag = 1;
    
    %Updating the xdist and zdist arrays
    
    LevBoundDistx(nlev,:) = lev(nlev,2) - bound(:,2)';
    LevBoundDistz(nlev,:) = lev(nlev,3) - bound(:,3)';
    
    LevTevDistx(nlev,:) = (lev(nlev,2)-tev(1:ntev,2))';
    LevTevDistz(nlev,:) = (lev(nlev,3)-tev(1:ntev,3))';
    
    LevLevDistx(:,nlev) = lev(nlev,2)-lev(1:nlev,2);
    LevLevDistz(:,nlev) = lev(nlev,3)-lev(1:nlev,3);
    LevLevDistx(nlev,:)=lev(1:nlev,2)-lev(nlev,2);
    LevLevDistz(nlev,:)=lev(1:nlev,3)-lev(nlev,3);
    
    tevIter_prev = 0;
    tevIter = -0.01;
    levIter_prev = 0;
    levIter = -0.01;
    kelv_prev = 0;
    kutta_prev = 0;
    
    for iter=2:100
        if (iter == 100)
        error('2D NR iteration failed')
        end
        
        lev(nlev,1) = levIter_prev;
        tev(ntev,1) = tevIter;
        
        calcDownwash
        
        kelv_tev = kelv_enf + bound_circ + sum(lev(1:nlev,1)) + sum(tev(1:ntev,1));
        kutta_tev = aterm(1) - lesp_cond;
        dkelv_tev = (kelv_tev-kelv_prev)/(tevIter-tevIter_prev);
        dkutta_tev = (kutta_tev-kutta_prev)/(tevIter-tevIter_prev);
    
        %Advancing with lev strength
        lev(nlev,1)=levIter;
        tev(ntev,1)=tevIter_prev;
        
        calcDownwash
        
        kelv_lev = kelv_enf + bound_circ + sum(lev(1:nlev,1)) + sum(tev(1:ntev,1));
        kutta_lev = aterm(1) - lesp_cond;
        dkelv_lev = (kelv_lev-kelv_prev)/(levIter-levIter_prev);
        dkutta_lev = (kutta_lev-kutta_prev)/(levIter-levIter_prev);
        
        %Advancing with both
        lev(nlev,1)=levIter;
        tev(ntev,1)=tevIter;
        
        calcDownwash
        
        kelv = kelv_enf + bound_circ + sum(lev(1:nlev,1)) + sum(tev(1:ntev,1));
        kutta = aterm(1) - lesp_cond;
        
        if (abs(kelv) < eps && abs(kutta) < eps)
            break
        end
        
        tevIter_prev = tevIter;
        levIter_prev = levIter;
        
        tevIter = tevIter-((1/(dkelv_tev*dkutta_lev-dkelv_lev*dkutta_tev))*...
                          ((dkutta_lev*kelv)-(dkelv_lev*kutta)));
    
        levIter = levIter-((1/(dkelv_tev*dkutta_lev-dkelv_lev*dkutta_tev))*...
                          ((-dkutta_tev*kelv)+(dkelv_tev*kutta)));
        kelv_prev = kelv;
        kutta_prev = kutta;
    end
else
    levflag = 0;
end

%% Calculate Fourier Terms and Bound Vorticity

for i = 3:45
    aterm(i) = sum(((downwash(2:70).*cos((i-1).*theta(2:70)) + downwash(1:69).*cos((i-1).*theta(1:69)))./2).*dtheta);
    aterm(i)=(2./(u_ref*pi))*aterm(i);
end

%% Set previous values of aterm to be used for derivatives in next time step
aterm_prev(1:4)=aterm(1:4);

%% Calculate bound_vortex strengths

for i=1:70
   gamma(i) = aterm(1)*(1+cos(theta(i))) +...
              sum(aterm(2:45).*sin((1:44)'.*theta(i)).*sin(theta(i)));
   gamma(i) = gamma(i)*u_ref*c;

end

bound_int(2:70,1) = ((gamma(2:70)+gamma(1:69))./2).*dtheta;
bound_int(2:70,2) = (bound(2:70,2)+bound(1:69,2))./2;
bound_int(2:70,3) = (bound(2:70,3)+bound(1:69,3))./2;

TevBintDistx = tev(1:ntev,2)*ones(1,70) - ones(ntev,1)*bound_int(:,2)';
TevBintDistz = tev(1:ntev,3)*ones(1,70) - ones(ntev,1)*bound_int(:,3)';
if nlev>0
    LevBintDistx = lev(1:nlev,2)*ones(1,70) - ones(nlev,1)*bound_int(:,2)';
    LevBintDistz = lev(1:nlev,3)*ones(1,70) - ones(nlev,1)*bound_int(:,3)';
end

%% Wake rollup

uind_tev = zeros(expect_vort,1);
wind_tev = zeros(expect_vort,1);
uind_lev = zeros(expect_vort,1);
wind_lev = zeros(expect_vort,1);

for i=1:ntev
    % TEV-TEV
    %sel = [1:i-1,i+1:ntev];
    dist = 2*pi*sqrt(v_core.^4 + ...
                (TevTevDistx(:,i).^2 + TevTevDistz(:,i).^2).^2);
    uind_int = (tev(1:ntev,1).*-TevTevDistz(:,i))./dist;
    uind_int(i) = [];
    uind_tev(i) = sum(uind_int);
    wind_int = (-tev(1:ntev,1).*-TevTevDistx(:,i))./dist;
    wind_int(i) = [];
    wind_tev(i) = sum(wind_int);
    
    %TEV-BOUND
    dist = 2*pi*sqrt(v_core.^4 + ...
                (TevBintDistx(i,:)'.^2 + TevBintDistz(i,:)'.^2).^2);
    uind_tev(i) = uind_tev(i) + sum((bound_int(:,1).*TevBintDistz(i,:)')./dist);
    wind_tev(i) = wind_tev(i) + sum((-bound_int(:,1).*TevBintDistx(i,:)')./dist);
    
    %TEV-LEV
    if nlev > 0
    dist = 2*pi*sqrt(v_core.^4 + ...
                (LevTevDistx(1:nlev,i).^2 + LevTevDistz(1:nlev,i).^2).^2);
    uind_tev(i) = uind_tev(i) + sum((lev(1:nlev,1).*-LevTevDistz(1:nlev,i))./dist);
    wind_tev(i) = wind_tev(i) + sum((-lev(1:nlev,1).*-LevTevDistx(1:nlev,i))./dist);
    end
end
    
for i=1:nlev
    % LEV-LEV
    dist = 2*pi*sqrt(v_core.^4 + ...
                (LevLevDistx(1:nlev,i).^2 + LevLevDistz(1:nlev,i).^2).^2);
    uind_int = (lev(1:nlev,1).*-LevLevDistz(1:nlev,i))./dist;
    uind_int(i) = [];
    uind_lev(i) = sum(uind_int);
    wind_int = (-lev(1:nlev,1).*-LevLevDistx(1:nlev,i))./dist;
    wind_int(i) = [];
    wind_lev(i) = sum(wind_int);
    
    % LEV-TEV
    dist = 2*pi*sqrt(v_core.^4 + ...
                (LevTevDistx(i,:)'.^2 + LevTevDistz(i,:)'.^2).^2);
    uind_lev(i) = uind_lev(i) + sum((tev(1:ntev,1).*LevTevDistz(i,:)')./dist);
    wind_lev(i) = wind_lev(i) + sum((-tev(1:ntev,1).*LevTevDistx(i,:)')./dist);
    
    % LEV-BOUND
    dist = 2*pi*sqrt(v_core.^4 + ...
                (LevBintDistx(i,:)'.^2 + LevBintDistz(i,:)'.^2).^2);
    uind_lev(i) = uind_lev(i) + sum((bound_int(:,1).*LevBintDistz(i,:)')./dist);
    wind_lev(i) = wind_lev(i) + sum((-bound_int(:,1).*LevBintDistx(i,:)')./dist);
end
       
%% Consider Newmark beta evolution here (to advect vortices)

tev(1:ntev,2) = tev(1:ntev,2)+(uind_tev(1:ntev,1).*dt);
tev(1:ntev,3) = tev(1:ntev,3)+(wind_tev(1:ntev,1).*dt);

if nlev>0
    lev(1:nlev,2)=lev(1:nlev,2)+(uind_lev(1:nlev,1).*dt);
    lev(1:nlev,3)=lev(1:nlev,3)+(wind_lev(1:nlev,1).*dt);
end

%% Remove TEV and LEV that have crossed a certain distance

if (tev(1,2)-bound(70,2)>del_dist)
    kelv_enf=kelv_enf+tev(1,1);
    tev(1:ntev-1,:)=tev(2:ntev,:);
    ntev=ntev-1;
end
if (nlev>0 && lev(1,2)-bound(70,2)>del_dist)
    kelv_enf=kelv_enf+lev(1,1);
    lev(1:nlev-1,:)=lev(2:nlev,:);
    nlev=nlev-1;
end

%% Load coefficient calculation (nondimensional units)
cnc=(2*pi*((u_ref*cos(alpha)/u_ref)+(hdot*sin(alpha)/u_ref))*...
    (aterm(1)+(aterm(2)/2)));
cnnc=(2*pi*((3*c*adot(1)/(4*u_ref))+(c*adot(2)/(4*u_ref))+...
    (c*adot(3)/(8*u_ref))));
cs=2*pi*aterm(1)*aterm(1);
%The components of normal force and moment from induced velocities are
%calulcated in dimensional units and nondimensionalized later
nonl = sum((((uind(2:70)*cos(alpha))-(wind(2:70)...
    *sin(alpha))).*bound_int(2:70,1)));
nonl_m = sum((((uind(2:70)*cos(alpha))-(wind(2:70)...
    *sin(alpha))).*(x(2:70)).*bound_int(2:70,1)));
nonl=nonl*(2/(u_ref*u_ref*c));
nonl_m=nonl_m*(2/(u_ref*u_ref*c*c));

cn=cnc+cnnc+nonl;
cl=cn*cos(alpha)+cs*sin(alpha);
cd=cn*sin(alpha)-cs*cos(alpha);

%Pitching moment is clockwise or nose up positive
cm=-((2*pi*((u_ref*cos(alpha)/u_ref)+(hdot*sin(alpha)/u_ref))...
    *(((aterm(1)/4)+(aterm(2)/4)-(aterm(3)/8))...
    -(pvt*(aterm(1)+(aterm(2)/2)))))+...
    ((pi*c/u_ref)*(((7*adot(1)/8)+(3*adot(2)/8)+(adot(3)/8)-(adot(4)/32))...
    -(pvt*((3*adot(1)/2)+(adot(2)/2)+(adot(3)/4)))))+(nonl_m));

%% Update States
LocalState.ntev = ntev;
LocalState.nlev = nlev;
LocalState.tev = tev;
LocalState.lev = lev;
LocalState.dist_wind = dist_wind;
LocalState.kelv_enf = kelv_enf;
LocalState.aterm_prev = aterm_prev;
LocalState.levflag = levflag;
LocalState.a0 = aterm(1);

%% Local Functions

%% calcDownwash
function calcDownwash
    if nlev == 0
        disttev = TevBoundDistx.^2 + TevBoundDistz.^2;
        tevMatrix = tev(1:ntev,1)*ones(1,70);
        
        
        uind = sum(tevMatrix.*(-TevBoundDistz)./...
                   (2*pi*sqrt(v_core.^4 + disttev.^2)),1)';
        wind = sum(-tevMatrix.*(-TevBoundDistx)./...
                   (2*pi*sqrt(v_core.^4 + disttev.^2)),1)';
    else
        
        disttev = TevBoundDistx.^2 + TevBoundDistz.^2;
        distlev = LevBoundDistx(1:nlev,:).^2 + LevBoundDistz(1:nlev,:).^2;
        tevMatrix = tev(1:ntev,1)*ones(1,70);
        levMatrix = lev(1:nlev,1)*ones(1,70);
        
        
        uind = sum(tevMatrix.*(-TevBoundDistz)./...
                   (2*pi*sqrt(v_core.^4 + disttev.^2)),1)' + ...
               sum(levMatrix.*(-LevBoundDistz(1:nlev,:))./...
                   (2*pi*sqrt(v_core.^4 + distlev.^2)),1)';
        wind = sum(-tevMatrix.*(-TevBoundDistx)./...
                   (2*pi*sqrt(v_core.^4 + disttev.^2)),1)' + ...
               sum(-levMatrix.*(-LevBoundDistx(1:nlev,:))./...
                   (2*pi*sqrt(v_core.^4 + distlev.^2)),1)';
    end
        
    downwash = (-u_ref*sin(alpha))+(-uind*sin(alpha))...
          +(hdot*cos(alpha))+(-wind*cos(alpha))...
          +(-alphadot*(x-pvt*c));
    
    aterm(1) = sum(((downwash(2:70) + downwash(1:69))./2)*dtheta);
    aterm(2) = sum(((downwash(2:70).*cos(theta(2:70))...
               + downwash(1:69).*cos(theta(1:69)))./2).*dtheta);    

    aterm(1) = (-1./(u_ref*pi)).*aterm(1);
    aterm(2) = (2./(u_ref*pi)).*aterm(2);
    bound_circ = u_ref*c*pi*(aterm(1)+(aterm(2)/2)); 
end

end
% linear dynamics function
% X = [eta; etad]
function [dXdt] = Dynamics_LDVM_flatplate(t,X,tag,q,nm,ns,mmod,U,rho,Phi,vec_b,vec_dy)
global fix_p sch_p clvec cmvec a0vec Qvec state

eta = X(1:nm,1);
etad = X(nm+1:2*nm,1);
sch_p_loc = struct(sch_p);

if tag == 1
    clvecloc = zeros(ns,1);
    cmvecloc = zeros(ns,1);
    a0vecloc = zeros(ns,1);
    stateloc = zeros(ns,2);
    
    Qvecloc = zeros(nm,ns); % vector of generalised forces
    pdyn = 0.5*rho*U^2; % dynamic pressure
    
    sch_p_in = sch_p;
    fix_p_in = fix_p;
    
    parfor jj = 1:ns
        vecphys = Phi(:,:,jj)*eta;
        vecphys_dot = Phi(:,:,jj)*etad;
        h = -vecphys(1);
        alpha = +vecphys(2) + q;
        hdot = -vecphys_dot(1);
        alphadot = +vecphys_dot(2);
        
        [cl, ~, cm,sch_p_loc(jj)] = LDVM_FlatPlate(alpha, alphadot, h, hdot,...
                                                    sch_p_in(jj),fix_p_in);
        Qvecloc(:,jj) = Phi(:,:,jj)'*pdyn*[-cl*vec_dy(jj)*(2*vec_b(jj));...
                                          cm*vec_dy(jj)*(2*vec_b(jj))^2;
                                          0];
        
        clvecloc(jj) = cl;
        cmvecloc(jj) = cm;
        a0vecloc(jj) = sch_p_loc(jj).a0;
        stateloc(jj,:) = [alpha, h];
    end
    
    sch_p = sch_p_loc;
    Qvec = sum(Qvecloc,2);
    clvec = clvecloc;
    cmvec = cmvecloc;
    a0vec = a0vecloc;
    state = stateloc;
    
end

if t == 0, Qvec = 0*Qvec; end

etadd = (mmod.m_eta_s_etadd)\...
    ((- mmod.m_eta_s_etad)*etad + ...
    (- mmod.m_eta_s_eta)*eta + ...
    Qvec);


dXdt(1:nm,1) = etad;
dXdt(nm+1:2*nm,1) = etadd;



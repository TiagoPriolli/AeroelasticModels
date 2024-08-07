% linear dynamics function
% X = [eta; etad; lbda]
function [dXdt] = dynamics_lin(t,X,q,nm,ns,mmod)

global qvec

eta = X(1:nm,1);
etad = X(nm+1:2*nm,1);
lbda = X(2*nm+1:2*nm+2*ns,1);
etadd = (mmod.m_eta_s_etadd - mmod.m_eta_a_etadd)\...
              ((mmod.m_eta_a_etad - mmod.m_eta_s_etad)*etad + ...
               (mmod.m_eta_a_eta - mmod.m_eta_s_eta)*eta + ...
               mmod.m_eta_a_lbda*lbda + ...
               mmod.m_eta_q*q);

lbdad = mmod.m_lbda_etadd*etadd + ... 
        mmod.m_lbda_etad*etad + ...
        mmod.m_lbda_eta*eta + ...
        mmod.m_lbda_lbda*lbda + ...
        mmod.m_lbda_q*q;

qvec = mmod.m_eta_a_etad*etad + mmod.m_eta_a_eta*eta +...
       mmod.m_eta_a_etadd*etadd + mmod.m_eta_a_lbda*lbda + ...
       mmod.m_eta_q*q;

dXdt(1:nm,1) = etad;
dXdt(nm+1:2*nm,1) = etadd;
dXdt(2*nm+1:2*nm+2*ns,1) = lbdad;
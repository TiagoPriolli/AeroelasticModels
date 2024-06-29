% function for finding linear dynamics equilibrium point
function f = trimming_lin(eta,q,nm,ns,mmod)
X = [eta;zeros(nm,1);zeros(2*ns,1)];
dXdt = dynamics_lin(0,X,q,nm,ns,mmod);
f = dXdt(nm+1:2*nm,1);

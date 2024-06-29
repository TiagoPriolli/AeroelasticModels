% function for finding spanwise lift and pitch moment distribution per unit
% span
% - X = [eta;etad;lbda]
% - dXdt = dX/dt (you may obtain via dynamics_lin function)
function [f_tot,f_str,f_str_nc,f_str_c,f_alpha] = aed_dist_lin(X,dXdt,alpha,nm,ns,maed,Phi)
eta = X(1:nm,1);
etad = X(nm+1:2*nm,1);
lbda = X(2*nm+1:2*nm+2*ns,1);
etadd = dXdt(nm+1:2*nm,1);

f_tot = zeros(3,ns);
f_str = zeros(3,ns);
f_str_nc = zeros(3,ns);

for jj = 1:ns
    f_tot(:,jj) = maed.A1(:,:,jj)*Phi(:,:,jj)*etadd + ...
        maed.A2(:,:,jj)*Phi(:,:,jj)*etad + ...
        maed.A3(:,:,jj)*Phi(:,:,jj)*eta + ...
        maed.A4(:,:,jj)*lbda(2*jj-1:2*jj) + ...
        maed.A3(:,:,jj)*[0;1;0]*alpha;
    
    f_str(:,jj) = maed.A1(:,:,jj)*Phi(:,:,jj)*etadd + ...
        maed.A2(:,:,jj)*Phi(:,:,jj)*etad + ...
        maed.A3(:,:,jj)*Phi(:,:,jj)*eta + ...
        maed.A4(:,:,jj)*lbda(2*jj-1:2*jj);
   
    f_str_nc(:,jj) = maed.A1_nc(:,:,jj)*Phi(:,:,jj)*etadd + ...
           maed.A2_nc(:,:,jj)*Phi(:,:,jj)*etad + ...
           maed.A3_nc(:,:,jj)*Phi(:,:,jj)*eta;
       
       
end
f_str_c = f_str - f_str_nc;
f_alpha = f_tot - f_str;
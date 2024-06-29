% Finite-Elements Model
n_elem = 70;
n_nodes = n_elem + 1;

% Strip Theory Model
ns = 20;

% Wing
l = 0.40;

% Elements/Strips Positions
Struc_dy = l/n_elem;
Aero_dy = l/ns;

Struc_vec_y = linspace(0,l,n_nodes);
Aero_vec_y = Aero_dy/2:Aero_dy:l-Aero_dy/2;

% Matriz G
% GDLs -3 e 5 são correspondentes à h e alpha

G = zeros(3*ns,6*n_nodes);
jj = 1; %Strip Counter
for ii = 1:3:3*ns
    FiltNext = Aero_vec_y(jj) < Struc_vec_y;
    Inext = n_nodes - nnz(FiltNext) + 1;
    [Vnext] = min(Struc_vec_y(FiltNext));
    FiltBef = Aero_vec_y(jj) > Struc_vec_y;
    Ibef = nnz(FiltBef);
    [Vbef] = max(Struc_vec_y(FiltBef));
    
    PropNext = abs(Vnext - Aero_vec_y(jj))./Struc_dy;
    PropBef = abs(Vbef - Aero_vec_y(jj))./Struc_dy;
    
    G(ii,[(Ibef-1)*6+3 (Inext-1)*6+3]) = [PropBef PropNext]; %h
    G(ii+1,[(Ibef-1)*6+5 (Inext-1)*6+5]) = [PropBef PropNext]; %Alpha
    G(ii+2,[Ibef Inext]) = [0 0]; %Beta
    jj = jj + 1;
end


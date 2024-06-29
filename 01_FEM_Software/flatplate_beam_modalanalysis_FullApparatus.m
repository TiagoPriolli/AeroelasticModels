% Modal analysis of a clamped plate with beam elements
% coordinate system:
%
% ^    ^ z
%y \   |
%   \__|__________________________
%    \ |                          \
%     \|_________ x                \
%      \                            \
%       \____________________________\
%        \
%
% 
% Output (modes) in a reference system -90deg in z
%  - x pointing to the trailing edge
%  - y pointing to the wing tip
%  - z pointing up

%%

clear all                                                                   % exclui todas as vari?veis
close all                                                                   % fecha todas as figuras
clc                                                                         % limpa a tela de comando


%% Model Parameters
lastro = 1;                                                                 % par?metro auxiliar que indica exist?ncia (1) ou aus?ncia (0) do lastro
offset_lastro = 5/1000;                                                     % deslocamento do centro de massa do lastro com rela??o ? linha de centro
n_eig = 30;                                                                 % n?mero de autovalores desejado, ou seja, n?mero de modos de vibra??o desejado

%% wing parameters
wing.Lw = 400e-3;%350e-3;                                                                % envergadura da asa [m]
wing.cw = 27e-3;%40e-3;                                                                 % corda da asa [m]
wing.tw = 0.8e-3;                                                                  % espessura da asa [m]
wing.ro = 2780;                                                                  % densidade do material da asa (aluminio 2024)
wing.E = 73.1e9;                                                                 % m?dulo de Young do material
wing.v = 0.33;                                                                   % coeficiente de Poisson do material

wing.mw = wing.ro*wing.Lw*wing.cw*wing.tw;                                                           % massa da asa dada pelo produto da densidade com o volume da asa

wing.n_elem = 70;                                                                % n?mero de "elementos de viga", divis?o realizado ao longo da envergadura
wing.L_e = wing.Lw/wing.n_elem;                                                       % comprimento de cada elemento

%% balance parameters
n_mode = 1;

balance.Lw = 475e-3+25e-3;                                                                % envergadura da asa [m]
balance.cw = ones(100,1).*45e-3*sqrt(pi/4);                                                                 % corda da asa [m]
balance.cw(1:13) = 48.5e-3;
balance.cw(14:16) = 51.5e-3*sqrt(pi/4);
balance.cw(17:25) = 60.5e-3*sqrt(pi/4);
balance.cw([42:45 82:85]) = 6.5e-3; 
balance.tw = ones(100,1).*45e-3*sqrt(pi/4);                                                                  % espessura da asa [m]
balance.tw(1:13) = 18.5e-3;
balance.tw(14:16) = 51.5e-3*sqrt(pi/4);
balance.tw(17:25) = 60.5e-3*sqrt(pi/4);
balance.tw([35:38 75:78]) = 6.5e-3;
balance.tw(90:100) = 12.8e-3;
balance.ro = 2700;                                                                  % densidade do material da asa (aluminio 2024)
balance.E = 70.1e9*0.1;                                                                 % m?dulo de Young do material
balance.v = 0.3;                                                                   % coeficiente de Poisson do material
balance.n_elem = 95+5; 
balance.mw = balance.ro.*(balance.Lw/balance.n_elem).*balance.cw.*balance.tw;                                                           % massa da asa dada pelo produto da densidade com o volume da asa

                                                               % n?mero de "elementos de viga", divis?o realizado ao longo da envergadura
balance.L_e = balance.Lw/balance.n_elem; 

%% Wing Elements
%Elementos centrais

m_e = wing.mw/wing.n_elem;                                                            % massa de cada elemento central
Ixx_e = 1/12*m_e*(wing.cw^2+wing.tw^2);                                               % momento de in?rcia com rela??o ao eixo x
Iyy_e = 1/12*m_e*(wing.L_e^2+wing.tw^2);                                              % momento de in?rcia com rela??o ao eixo y
Izz_e = 1/12*m_e*(wing.L_e^2+wing.cw^2);                                              % momento de in?rcia com rela??o ao eixo z

wing.M_e_c = [m_e*eye(3) zeros(3,3)                                              % matriz de massa dos elementos centrais
    zeros(3,3) diag([Ixx_e Iyy_e Izz_e])];

% Elemento da ponta

m_e_t      = wing.mw/wing.n_elem/2;                                                   % massa do elemento da ponta
Ixx_e_t_cg = 1/12*m_e_t*(wing.cw^2+wing.tw^2);                                        % momento de in?rcia com rela??o ao eixo x
Iyy_e_t_cg = 1/12*m_e_t*((wing.L_e/2)^2+wing.tw^2);                                   % momento de in?rcia com rela??o ao eixo y
Izz_e_t_cg = 1/12*m_e_t*((wing.L_e/2)^2+wing.cw^2);                                   % momento de in?rcia com rela??o ao eixo z

M_e_t_cg = [m_e_t*eye(3)                                  zeros(3,3)        % matriz de in?cia dos elementos da ponta
              zeros(3,3)    diag([Ixx_e_t_cg Iyy_e_t_cg Izz_e_t_cg])];

skew = @(vec)([      0   -vec(3)     vec(2)                                 % Produto vetorial
                 vec(3)       0     -vec(1)
                -vec(2)   vec(1)         0]);

T_node_t_2_cg = [       eye(3)     -skew([-wing.L_e/4 0 0].')                    % operador de deslocamento do eixo de rota??o do centro de massa para a dist?ncia
                    zeros(3,3)                    eye(3)];                  % 'wing.L_e/4' localizado na extremidade do elemento.

M_e_t_node = T_node_t_2_cg.'*M_e_t_cg*T_node_t_2_cg;                        % Matriz de massa com centro de rota??o deslocado para o n? localizado na
                                                                            % extremidade do elemento.

% Lastro

m_lastro   = 0.0328;                                                       % massa do lastro
Ixx_lastro = 1.858e-5;                                                      % momento de in?rcia com rela??o ao eixo x
Izz_lastro = 1.858e-5;                                                      % momento de in?rcia com rela??o ao eixo z

M_lastro_cg = [m_lastro*eye(3)                       zeros(3,3)             % matriz de massa do lastro
                    zeros(3,3)  diag([Ixx_lastro 0 Izz_lastro])];

T_node_t_2_lastro = [       eye(3)      -skew([0 offset_lastro 0].')        % Ajuste de offset do lastro ao longo da corda
                        zeros(3,3)                            eye(3)];

M_lastro_node = T_node_t_2_lastro.'*M_lastro_cg*T_node_t_2_lastro;          % Matriz de massa do lastro com offset realizado, o qual foi determinado no bloco
                                                                            % de 'Par?metros f?sicos' dessa rotina.

wing.M_e_t_node = M_e_t_node + lastro*M_lastro_node;                             % Adicionando a massa do lastro ao elemento da ponta da asa. Vale ressaltar que o
                                                                            % par?metro 'lastro' define a exist?ncia ou n?o do lastro, como definido no bloco
                                                                            % 'Par?metros f?sicos' dessa rotina.

%% Balança
% Elemento da raiz

m_e_r = balance.mw(1)/2;                                                        % massa do elemento da ra?z
Ixx_e_r_cg = 1/12*m_e_r*(balance.cw(1).^2+balance.tw(1).^2);                                        % momento de in?rcia com rela??o ao eixo x
Iyy_e_r_cg = 1/12*m_e_r*((balance.L_e/2)^2+balance.tw(1).^2);                                   % momento de in?rcia com rela??o ao eixo y
Izz_e_r_cg = 1/12*m_e_r*((balance.L_e/2)^2+balance.cw(1).^2);                                   % momento de in?rcia com rela??o ao eixo z

M_e_r_cg = [m_e_r*eye(3)                                  zeros(3,3)        % matriz de massa do elemento da ra?z
              zeros(3,3)    diag([Ixx_e_r_cg Iyy_e_r_cg Izz_e_r_cg])];

skew = @(vec)([     0   -vec(3)     vec(2)                                  % Produto vetorial
                vec(3)       0     -vec(1)                                  % a x b = skew(a)*b
               -vec(2)   vec(1)         0]);                                % Esse tipo de comando (@(vec)) ? o mesmo que fazer que:
                                                                            % function skew_mat = skew(vec)
                                                                            % skew_mat = [     0   -vec(3)     vec(2)
                                                                            %              vec(3)       0     -vec(1)
                                                                            %             -vec(2)   vec(1)         0];

T_node_r_2_cg = [     eye(3)    -skew([wing.L_e/4 0 0].')                        % operador de deslocamento do eixo de rota??o do centro de massa para a dist?ncia
                  zeros(3,3)                  eye(3)];                      % 'wing.L_e/4' localizado na extremidade do elemento. Esse operador ? o [d] mencionado
                                                                            % na descri??o desse bloco.
                                                                            
balance.M_e_r_node = T_node_r_2_cg.'*M_e_r_cg*T_node_r_2_cg;                        % Matriz de massa com centro de rota??o deslocado para o n? localizado na
                                                                            % extremidade do elemento.

%Elementos centrais
balance.M_e_c = zeros(6,6,balance.n_elem);
for ii=1:balance.n_elem
    m_e = balance.mw(ii);                                                            % massa de cada elemento central
    Ixx_e = 1/12*m_e*(balance.cw(ii)^2+balance.tw(ii)^2);                                               % momento de in?rcia com rela??o ao eixo x
    Iyy_e = 1/12*m_e*(balance.L_e^2+balance.tw(ii)^2);                                              % momento de in?rcia com rela??o ao eixo y
    Izz_e = 1/12*m_e*(balance.L_e^2+balance.cw(ii)^2);                                              % momento de in?rcia com rela??o ao eixo z

    balance.M_e_c(:,:,ii) = [m_e*eye(3) zeros(3,3)                                              % matriz de massa dos elementos centrais
    zeros(3,3) diag([Ixx_e Iyy_e Izz_e])];
end
                                                                            
                                                                            
                                                                            


%% Matriz de massa da estrutura:

wing.n_nodes = wing.n_elem+1;                                                         % n?mero de elementos presentes no modelo, definido no in?cio dessa rotina. Vale
balance.n_nodes = balance.n_elem+1;                                                                            % ressaltar que h? um elemento adicional por ser o elemento da ra?z.

n_nodes = wing.n_nodes + balance.n_nodes-1;
M = zeros(6*n_nodes,6*n_nodes);                                             % inicializa??o da matriz de massa. Como cada n? apresenta 6 gras de liberdade,
                                                                            % logo a estrutura sera dotada de "6 x n" graus de liberdade, por conter 'n'
                                                                            % elementos  de 6 graus de liberdade cada um: 3 de transla??o e 3 de rota??o
i_node = 1;                                                                 % ?ndice de identifica??o do primeiro n?
M(6*(i_node-1)+(1:6),6*(i_node-1)+(1:6)) = balance.M_e_r_node;                      % inclus?o do primeiro bloco de massa, o qual ? referente ao elemento da ra?z

for i_node = 2:balance.n_nodes-1
    M(6*(i_node-1)+(1:6),6*(i_node-1)+(1:6)) = balance.M_e_c(:,:,i_node);                       % inclus?o dos demais blocos de massa, os quais s?o referente aos elementos do
end                                                                         % meio da asa.

for i_node = balance.n_nodes:n_nodes-1
    M(6*(i_node-1)+(1:6),6*(i_node-1)+(1:6)) = wing.M_e_c;                       % inclus?o dos demais blocos de massa, os quais s?o referente aos elementos do
end                                                                         % meio da asa.

i_node = n_nodes;                                                           % ?ndice de identifica??o do ?ltimo n?
M(6*(i_node-1)+(1:6),6*(i_node-1)+(1:6)) = wing.M_e_t_node;                      % inclus?o do ?ltimo bloco de massa, o qual ? referente ao elemento da ponta

M = sparse(M);                                                              % Compacta uma matriz com um dimens?o muito grande a qual cont?m um alto n?mero
                                                                            % de elemento nulos. Essa compacta??o ? realizada de maneira a guardar os
                                                                            % n?o nulos e sua respectiva posi??o na matriz de linha e coluna.

%% Matriz de rigidez do elemento ASA

G = wing.E/(2*(1+wing.v));                                                            % M?dulo de cisalhamento

A   = wing.cw*wing.tw;                                                                % ?rea de se??o transvers?o da asa
Iy  = 1/12*wing.cw*wing.tw^3;                                                         % Momento de in?rcia de ?rea da asa em torno do eixo y no centro de massa
Iz  = 1/12*wing.tw*wing.cw^3;                                                         % Momento de in?rcia de ?rea da asa em torno do eixo z no centro de massa

J = 1/3*wing.cw*wing.tw^3;                                                             % Momento de in?rcia da asa calculado por J = 1/3*wing.cw*wing.tw^3

a  = wing.E*A/wing.L_e;                                                               % Vari?veis auxiliares para a montagem da matriz de rigidez
bz = wing.E*Iz/wing.L_e^3;                                                            %
by = wing.E*Iy/wing.L_e^3;                                                            %
t  = G*J/wing.L_e;                                                                    %

K_e = zeros(12,12);                                                                   % Inicializa??o da matriz de rigidez

% Ver equa??o 8.29 (Introduction to Finite Elements in Engineering;)

K_e(1,1) = a;
K_e(2,2) = 12*bz;
K_e(3,3) = 12*by;
K_e(4,4) = t;
K_e(5,3) = -6*by*wing.L_e;
K_e(5,5) = 4*by*wing.L_e^2;
K_e(6,2) = 6*bz*wing.L_e;
K_e(6,6) = 4*bz*wing.L_e^2;

K_e(7,1)  = -a;
K_e(8,2)  = -12*bz;
K_e(8,6)  = -6*bz*wing.L_e;
K_e(9,3)  = -12*by;
K_e(9,5)  = 6*by*wing.L_e;
K_e(10,4) = -t;
K_e(11,3) = -6*by*wing.L_e;
K_e(11,5) = 2*by*wing.L_e^2;
K_e(12,2) = 6*bz*wing.L_e;
K_e(12,6) = 2*bz*wing.L_e^2;

K_e(7,7)   = a;
K_e(8,8)   = 12*bz;
K_e(9,9)   = 12*by;
K_e(10,10) = t;
K_e(11,9)  = 6*by*wing.L_e;
K_e(11,11) = 4*by*wing.L_e^2;
K_e(12,8)  = -6*bz*wing.L_e;
K_e(12,12) = 4*bz*wing.L_e^2;

wing.K_e = K_e+tril(K_e,-1).';                                                   % Matriz de rigidez do elemento

%% Matriz de rigidez do elemento BALANCE

G = balance.E/(2*(1+balance.v));                                                            % M?dulo de cisalhamento

A   = balance.cw.*balance.tw;                                                                % ?rea de se??o transvers?o da asa
Iy  = 1/12.*balance.cw.*balance.tw.^3;                                                         % Momento de in?rcia de ?rea da asa em torno do eixo y no centro de massa
Iz  = 1/12.*balance.tw.*balance.cw.^3;                                                         % Momento de in?rcia de ?rea da asa em torno do eixo z no centro de massa

J = 1/3.*balance.cw.*balance.tw.^3;                                                             % Momento de in?rcia da asa calculado por J = 1/3*wing.cw*wing.tw^3

a  = balance.E.*A./balance.L_e;                                                               % Vari?veis auxiliares para a montagem da matriz de rigidez
bz = balance.E.*Iz./balance.L_e^3;                                                            %
by = balance.E.*Iy./balance.L_e^3;                                                            %
t  = G.*J./balance.L_e;                                                                    %

K_e = zeros(12,12);                                                                   % Inicializa??o da matriz de rigidez
balance.K_e = zeros(12,12,balance.n_elem);

% Ver equa??o 8.29 (Introduction to Finite Elements in Engineering;)
for ii=1:balance.n_elem
    K_e(1,1) = a(ii);
    K_e(2,2) = 12*bz(ii);
    K_e(3,3) = 12*by(ii);
    K_e(4,4) = t(ii);
    K_e(5,3) = -6*by(ii)*balance.L_e;
    K_e(5,5) = 4*by(ii)*balance.L_e^2;
    K_e(6,2) = 6*bz(ii)*balance.L_e;
    K_e(6,6) = 4*bz(ii)*balance.L_e^2;

    K_e(7,1)  = -a(ii);
    K_e(8,2)  = -12*bz(ii);
    K_e(8,6)  = -6*bz(ii)*balance.L_e;
    K_e(9,3)  = -12*by(ii);
    K_e(9,5)  = 6*by(ii)*balance.L_e;
    K_e(10,4) = -t(ii);
    K_e(11,3) = -6*by(ii)*balance.L_e;
    K_e(11,5) = 2*by(ii)*balance.L_e^2;
    K_e(12,2) = 6*bz(ii)*balance.L_e;
    K_e(12,6) = 2*bz(ii)*balance.L_e^2;

    K_e(7,7)   = a(ii);
    K_e(8,8)   = 12*bz(ii);
    K_e(9,9)   = 12*by(ii);
    K_e(10,10) = t(ii);
    K_e(11,9)  = 6*by(ii)*balance.L_e;
    K_e(11,11) = 4*by(ii)*balance.L_e^2;
    K_e(12,8)  = -6*bz(ii)*balance.L_e;
    K_e(12,12) = 4*bz(ii)*balance.L_e^2;

    balance.K_e(:,:,ii) = K_e+tril(K_e,-1).';                                                   % Matriz de rigidez do elemento
end

%% Matriz de rigidez da estrutura:

K = zeros(6*n_nodes,6*n_nodes);

for i_elem=1:balance.n_elem
    K(6*(i_elem-1)+(1:12),6*(i_elem-1)+(1:12)) = ...
        K(6*(i_elem-1)+(1:12),6*(i_elem-1)+(1:12))+balance.K_e(:,:,i_elem);
end

for i_elem=balance.n_elem+1:balance.n_elem+wing.n_elem
    K(6*(i_elem-1)+(1:12),6*(i_elem-1)+(1:12)) = ...
        K(6*(i_elem-1)+(1:12),6*(i_elem-1)+(1:12))+wing.K_e;
end

K = sparse(K);

% Modos de vibra??o

[eigvec,eigval] = eigs(K,M,n_eig,0.1);

% Ordenando os modos de vibra??o

[eig_val_ord,eig_order] = sort(diag(eigval));                               % ordena os autovalores do menor para o maior, armazenando o autovalor em si
                                                                            % (eig_val_ord) e sua posi??o no vetor (eig_order).

eigvec = eigvec(:,eig_order);                                               % ordena os autovetores em ordem crescente de acordo com a ordem obtida dos
                                                                            % autovalores, o qual est? armazenado no vetor 'eig_order'. Vale ressaltar que,
                                                                            % enquanto os autovalores s?o ordenados pela coluna, os autovetores s?o
                                                                            % ordenados pela linha. Por isso a ordem foi estabelecida como segundo argumento
                                                                            % do comando.

eigval = eigval(eig_order,eig_order);                                       % ordena os autovalores em ordem crescente em uma matriz diagonal

eigval = diag(eigval);                                                      % organiza os autovalores em um vetor coluna, em ordem crescente

shape_functions = eigvec;                                                   % nomeia os autovetores como "fun??es de forma" ou "shape functions", pois
                                                                            % s?o as formas modais ou formas de vibra??o.

% Normaliza??o

[u_max,v_max] = max(shape_functions);
[u_min,v_min] = min(shape_functions);

sinal = (abs(u_max)>=abs(u_min)).*1+(abs(u_max)<abs(u_min)).*(-1);

% M?ximo deslocamento
shape_functions_max = shape_functions./(ones(size(shape_functions,1),1)*...
    (max(abs(shape_functions)).*sinal));

MHH_max = shape_functions_max.'*M*shape_functions_max;
KHH_max = shape_functions_max.'*K*shape_functions_max;

vec_freq_max = sqrt(diag(KHH_max./MHH_max))/(2*pi);

% Massa
shape_functions_mass = shape_functions./(sqrt(ones(size(K,1),1)*...
    diag(shape_functions.'*M*shape_functions).'));

MHH_mass = shape_functions_mass.'*M*shape_functions_mass;
KHH_mass = shape_functions_mass.'*K*shape_functions_mass;

vec_freq_mass = sqrt(diag(KHH_mass./MHH_mass))/(2*pi);

% reshape(shape_functions_mass(:,7),6,3).'

% Engaste

n_dof = 6*n_nodes;                                                          % obten??o do n?mero de graus de liberdade total da asa

all_dofs = 1:n_dof;                                                         % vetor que cont?m todos os graus de liberdade

U = eye(n_dof);                                                             % matriz diagonal com dimens?o do n?mero de graus de liberdade total

s_node = 1;                                                                 % n? que ser? engastado
s_dofs = 6*(s_node-1)+1:6*(s_node-1)+6;                                     % graus de liberdade do engaste
f_dofs = setdiff(all_dofs,s_dofs);                                          % graus de liberdade livre

Us = U(:,s_dofs);
Us = sparse(Us);
Uf = U(:,f_dofs);
Uf = sparse(Uf);

Mff = Uf.'*M*Uf;    % Mff = M(f_dofs,f_dofs);
% Mff_mat = Uf.'*M*Uf; Mff_index = M(f_dofs,f_dofs); erro_Mff =
% sum(sum(abs(Mff_mat-Mff_index)));
Mfs = Uf.'*M*Us;    % Mfs = M(f_dofs,s_dofs);
Msf = Mfs.';        % Msf = M(s_dofs,f_dofs); Msf = Mfs.';
Mss = Us.'*M*Us;    % Mss = M(s_dofs,s_dofs);

Kff = Uf.'*K*Uf;
Kfs = Uf.'*K*Us;
Ksf = Mfs.';
Kss = Us.'*K*Us;

[eigvec,eigval] = eigs(Kff,Mff,n_eig-6,0.1);

[eig_val_ord,eig_order] = sort(diag(eigval));
eigvec = eigvec(:,eig_order);
eigval = eigval(eig_order,eig_order);
eigval = diag(eigval);

shape_functions = Uf*eigvec;    % shape_functions = zeros(n_dof,n_eig-6); shape_functions(f_dofs,:) = eigvec;
                                % Verificar que size(eigvec) ? n_dof-6
                                % linhas por n_eig-6 colunas

% Normaliza??o

[u_max,v_max] = max(shape_functions);
[u_min,v_min] = min(shape_functions);

sinal = (abs(u_max)>=abs(u_min)).*1+(abs(u_max)<abs(u_min)).*(-1);

% M?ximo deslocamento
shape_functions_max = shape_functions./(ones(size(shape_functions,1),1)*...
    (max(abs(shape_functions)).*sinal));

MHH_max = shape_functions_max.'*M*shape_functions_max;
KHH_max = shape_functions_max.'*K*shape_functions_max;

vec_freq_max = sqrt(diag(KHH_max./MHH_max))/(2*pi);

% Massa
shape_functions_mass = shape_functions./(sqrt(ones(size(K,1),1)*...
    diag(shape_functions.'*M*shape_functions).'));

MHH_mass = shape_functions_mass.'*M*shape_functions_mass;
KHH_mass = shape_functions_mass.'*K*shape_functions_mass;

vec_freq_mass = sqrt(diag(KHH_mass./MHH_mass))/(2*pi);

% rotating axes (-90deg in z)

shape_functions_mass_rot = 0*shape_functions_mass;
shape_functions_mass_rot(1:6:end,:) = -shape_functions_mass(2:6:end,:);
shape_functions_mass_rot(2:6:end,:) = shape_functions_mass(1:6:end,:);
shape_functions_mass_rot(3:6:end,:) = shape_functions_mass(3:6:end,:);
shape_functions_mass_rot(4:6:end,:) = -shape_functions_mass(5:6:end,:);
shape_functions_mass_rot(5:6:end,:) = shape_functions_mass(4:6:end,:);
shape_functions_mass_rot(6:6:end,:) = shape_functions_mass(6:6:end,:);

vec_omega = vec_freq_mass*2*pi; % vector of freqs in rad/s
vec_mu = 0*vec_omega + 1; % vector of modal masses
%vec_y = (0:wing.Lw/wing.n_elem:wing.Lw); % vector of y positions of the nodes
vec_y_balance = 0:balance.Lw/balance.n_elem:balance.Lw;
vec_y_wing = balance.Lw:wing.Lw/wing.n_elem:wing.Lw+balance.Lw;
vec_y = [vec_y_balance vec_y_wing(2:end)];

modalstr.wing.Lw = wing.Lw;
modalstr.wing.cw = wing.cw;
modalstr.wing.tw = wing.tw;
modalstr.wing.mw = wing.mw;
modalstr.shape = shape_functions_mass_rot;
modalstr.omega = vec_omega;
modalstr.mu = vec_mu;
modalstr.y = vec_y;

%  saving

save(['flatplate_',num2str(1000*wing.Lw),'_',num2str(1000*wing.cw),'_',num2str(10000*wing.tw),'_ballast_',num2str(lastro),'_offset_',num2str(1000*offset_lastro),'_nelem_',num2str(wing.n_elem),'.mat'],'modalstr')
%save(['flatplate_',num2str(1000*wing.Lw),'_',num2str(1000*wing.cw),'_',num2str(10000*wing.tw),'_ballast_',num2str(lastro),'_offset_',num2str(1000*offset_lastro),'_nelem_',num2str(wing.n_elem),'_Fixed','.mat'],...
%		'Uf','Us','Kff','Mff','Ksf','Msf','shape_functions_mass_rot','n_elem')


%% plotting modal shapes 
figure
for n_mode =1:5
UD_y = [1;1]*modalstr.y;
UD_x = [-1;1]*(0*modalstr.y+modalstr.wing.cw/2);
UD_z = 0*UD_y;

%n_mode = 6;

% deformed shape according to selected mode
D_x = [1;1]*modalstr.shape(1:6:end,n_mode)';
D_y = [1;1]*modalstr.shape(2:6:end,n_mode)'+...
      ([1;1]*modalstr.shape(6:6:end,n_mode)').*UD_x;
D_z = [1;1]*modalstr.shape(3:6:end,n_mode)'-...
    ([1;1]*modalstr.shape(5:6:end,n_mode)').*UD_x;

fat = 0.1*modalstr.wing.Lw/max(max(abs([D_x,D_y,D_z])));

subplot(1,5,n_mode)
surf(UD_x,UD_y,UD_z,'FaceColor','w','FaceAlpha',0.5,'EdgeColor','k');
hold on
surf(UD_x+fat*D_x,UD_y+fat*D_y,UD_z+fat*D_z,'FaceColor','m','FaceAlpha',0.5,'EdgeColor','k');
axis equal
xlabel('x');ylabel('y');zlabel('z');zlim([-0.04 0.04]);xlim([-0.0135 0.06])
title(['Mode ',num2str(n_mode),', f = ',num2str(modalstr.omega(n_mode)/2/pi),'Hz']);
view(-15,45)
end

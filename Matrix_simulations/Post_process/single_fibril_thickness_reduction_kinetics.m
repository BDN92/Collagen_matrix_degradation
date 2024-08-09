% function file tobe called in post processing code to solve for fibril
% diameters

function [t, df_final] = single_fibril_thickness_reduction_kinetics(df_0, lf_0, NE_0, time)


% lengths are in nm

% collagen structure related 
d_m = 1.5;                       % diameter of one tropocollagen unit in nm
l_m = 300 + 0.54 * 67;           % length of one tropocollagen unit in nm
V_tc = pi * (d_m/2)^2 * l_m;     % volume of one tropocollagen unit
S_tc = d_m * l_m;                % exposed rectangular area of one tropocollagen unit
l_b = 8;                         % length lattice site in nm
n_tc = 300/l_b;                  % number of binding sites in one tropocollagen unit 


% nimber of binding sites estimation on the surface
Ns_0 = n_tc * (pi * df_0 * lf_0)/S_tc;
Nv_0 = n_tc * (pi * df_0^2/2 * lf_0)/V_tc; 
NsV_0 = (pi * df_0 * lf_0)/S_tc;
NsR_0 = Ns_0 - NsV_0;


% time span and initial condition
tspan = [0 time];
y0 = [NE_0 NsR_0 NsV_0 0 0 0 0 0 0];  


% solve
[t,y] = ode23s(@(t,y) eqns(t,y, Ns_0, Nv_0), tspan, y0);


% checking balance of total enzyme and total sites
%Netotal = y(:,1) + y(:,4) + y(:,5) + y(:,6) + y(:,7) + y(:,8); 

Ns = y(:,2) + y(:,3) + y(:,4) + y(:,5) + 2 * y(:,6) + 2 * y(:,7) + 2 * y(:,8);
%Nstotal = Ns + 2 * y(:,9); 

df_final = Ns * (S_tc/n_tc) * (1/pi) * (1/lf_0);


function dydt = eqns(t, y, Ns_0, Nv_0)
dydt = zeros(9,1);

% rate constants for adsorption-desorption
k1f = 0.0154;  
k1b = 5.65 * 10^(-3);

k2f = k1f;    
k2b = k1b;

% rate constants for hopping related to regular sites
k3f = 3*10^(3);      
k3b = k1b;      

% rate constants for hopping related to vulnarable sites
k4f = k3f;  
k5f = k3f/45; 

% rate constant for local unwinding
k_B = 1.38 * 10^(-23);     % J/K
T = 273 + 37 ;
k_BT = k_B * T;
h_c = 6.626 * 10^(-34);    % J-s
kbThc = (k_B * T)/h_c;
residue = 56;    %(Perumal et al 1000 residue approx 300 nm)
N_A = 6.023 * 10^(23);
Ea = (1.9  * 10^3 * residue)/N_A ;

mul = 1;
%mul = mul +  0.013;
kwf = kbThc * exp(-mul*Ea/k_BT); 

% rate constant for cleavage
kc1 = 0.583;
kc2 = 0.472;
kcf = kc1;         

% definitions of variables

% N_E (1) 
% N_R^s (2)   % N_V^s (3)
% N_ER^s (4)   % N_EV^s (5)
% N_ERER^s (6)   % N_EVER^s (7)    % N_EVER*^s (8)
% N^p (9)

% designing function related to
% local stress assisted regular sites removal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_B = 1.38 * 10^(-23);     % J/K
T = 273 + 37 ;
k_BT = k_B * T;
h_c = 6.626 * 10^(-34);    % J-s
kBThc = (k_B * T)/h_c;
residue = 28;              %(Perumal et al 1054 residue approx 300 nm)
N_A = 6.023 * 10^(23);
Ea = (1.9  * 10^3 * residue)/N_A ;

a = 8 * 10^(-9);        
l_bond = 0.36 * 10^(-9);
gamma = k_BT/(a * 1.5 * 10^(-9));

nsites_remain = Ns_0 - 2 * y(9);

tau1 = 7;
tau2 = 1.4;

tau =  tau1 - tau2 + tau2 * exp(2 * y(9)/Nv_0);
F_min =  (gamma * a * kcf * tau) * (2 * y(7) * kwf * tau2);

f_min = ((F_min)/nsites_remain);
En = f_min * l_bond;

vf = exp( -  (mul*Ea - En)/(k_BT));
vb = exp( -  (mul*Ea + En)/(k_BT));

n_b = 1;
k_net =  n_b * kBThc * (vf - vb);
f =   y(2) * k_net ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ne = y(4) + y(5) + 0.0001;   % adding small number to avoid zero divisibility at t = 0
NN =  (nsites_remain -  y(4) - y(5) - 2 * (y(6) + y(7) + y(8)))/ne;
bb = 25;
enzyme_dia = 10;
b_site_l = 8;
b_site_b = 1.5;
cc = (pi * enzyme_dia^2)/(b_site_l * b_site_b);
phi = cc * NN/(bb+NN);


% ODEs set 3 test
dydt(1) = - k1f * y(1) + k1b * y(4) - k2f * y(1) + k2b * y(5) + kcf * y(8);

dydt(2) = - k1f * y(1) + k1b * y(4) - k3f * y(4) * phi + k3b * y(6) - k4f * y(5) * phi  - f;

dydt(3) = - k2f * y(1) + k2b * y(5) - k5f * y(4)  ;

dydt(4) = k1f * y(1) - k1b * y(4) - k3f * y(4) * phi  + k3b * y(6) - k5f * y(4)  ;

dydt(5) = k2f * y(1) - k2b * y(5) - k4f * y(5) * phi ;

dydt(6) = k3f * y(4) * phi - k3b * y(6);

dydt(7) = k4f * y(5) * phi  + k5f * y(4)   - kwf * y(7);

dydt(8) =  kwf * y(7) - kcf * y(8) ;

dydt(9) = kcf * y(8) + (1/2) * f ;

end


end


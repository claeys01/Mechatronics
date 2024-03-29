%%
clear all
close all



% Vessel Description

L = 0.97;               % Length overall [m]
B = 0.32;               % Beam overall [m]
T_f = 0.12;              % Draft forward [m]
T_a = 0.12;             % Draft aft   [m]
Displacement = 18;      % Displacement  [kg]
bow_thrust_LCG = 0.35;  % Bow Thruster LCG (w.r.t COG) [m]
scale = 33;             % Scale 1:33
ubgtr = 3 ;             % Upper bevel gear teeth ratio 13:39
tgrr = 3;               % Total gear reduction ratio
D_prop = 0.065;         % Propeller Diameter     
Aft_Thrust_LCG = -0.35; % Aft Thruster LCG (w.r.t COG)
Aft_Thrust_TCG = 0.065; % Aft Thruster TCG (port)





% System Inertia matrix
rho = 1000; %Water Density (not sure if 1000)
x_g = 0;
y_g = 0;
m = 18;     % Mass
I_z = (1/120)*pi*rho*L*B*T_f*(B^2+L^2); % Moment of Inertia

X_dotu = 1.2;
X_dotv = 0;
X_dotw = 0;
X_dotp = 0;
X_dotq = 0;
X_dotr = 0;

Y_dotu = 0; 
Y_dotv = 66.1; 
Y_dotw = 0;
Y_dotp = 0;
Y_dotq = 0;
Y_dotr = 0;

N_dotu = 0; 
N_dotv = 0;
N_dotw = 0;
N_dotp = 0;
N_dotq = 0;
N_dotr = 2.3; 



% Added Mass
M_a = -[X_dotu,X_dotv,X_dotr;
       Y_dotu, Y_dotv, Y_dotr;
       N_dotu,N_dotv,N_dotr];

% Rigid Body inertia
M_rb = [m,0,-m*y_g;
       0,m, m*x_g ;
       -m*y_g, m*x_g, I_z];

M = M_a + M_rb;




T1=1; % fixed for now
T2=1;
Tb=0;
alpha1= 0 ; %recht naar achteren
alpha2= 0 ; %recht naar achteren
tau_u = T1*cos(alpha1)+T2*cos(alpha2) ;
tau_v = T1*sin(alpha1)+T2*sin(alpha2)+Tb ;
%tau_r = -T1*cos(alpha1)*0.065+T2*cos(alpha1)*0.35+T2*cos(alpha2)-T2*sin(alpha2)+Tb*0.35;
tau_r = T1*cos(alpha1)*0.065 - T1*sin(alpha1)*0.35 - T2*cos(alpha2)*0.065 + T2*sin(alpha2)*0.35+Tb*0.35;

T_all = [1 0 0; 0 1 1; -Aft_Thrust_TCG -Aft_Thrust_LCG Aft_Thrust_LCG];

a_c = 0; % angle one of the current in our case (2D) this stays 0
b_c = 0; % angle two of the current can be changed in 2D plane
speed_c = 0.6; % speed of the current

vw = 0;
beta_w = 0;

flow=[a_c;b_c;speed_c];
tau = [tau_u;tau_v;tau_r];
%Environmental Force
tau_e = [0; 0; 0];


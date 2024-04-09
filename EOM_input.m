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

M_a = [X_dotu,X_dotv,X_dotr;
       Y_dotu, Y_dotv, Y_dotr;
       N_dotu,N_dotv,N_dotr];


% Rigid Body inertia

M_rb = [m,0,-m*y_g;
       0,m, m*x_g ;
       -m*y_g, m*x_g, I_z];

M = M_a + M_rb;


%Initial Trust allocation

T1=0.2; % fixed for now
T2=0.2;
Tb=0;
alpha1= 0 ; %straight
alpha2= 0 ; %straight
tau_x = T1*cos(alpha1)+T2*cos(alpha2) ;
tau_y = T1*sin(alpha1)+T2*sin(alpha2)+Tb ;
%tau_m = -T1*cos(alpha1)*0.065+T2*cos(alpha1)*0.35+T2*cos(alpha2)-T2*sin(alpha2)+Tb*0.35;
tau_m = T1*cos(alpha1)*0.065 - T1*sin(alpha1)*0.35 - T2*cos(alpha2)*0.065 + T2*sin(alpha2)*0.35+Tb*0.35;


%Trust allocation matrix

T_all = [1 0 0; 0 1 1; -Aft_Thrust_TCG -Aft_Thrust_LCG Aft_Thrust_LCG];


%Debug trust

tau_u = 0.2;
tau_v = 0.2;
tau_r = 0;


%Current

current_magnitude_percentage = 0; % a percentage that gives the current speed in percentage of the mx speed

a_c = 0; % angle one of the current in our case (2D) this stays 0
b_c = 0; % angle two of the current can be changed in 2D plane
speed_c = 1.5; % max speed of the current
flow=[a_c;b_c;speed_c];
tau = [tau_u;tau_v;tau_r];

%Wind

vw = 0;
beta_w = 0;


%Environmental Force

tau_e = [0; 0; 0];


%DC data

Ea=12;
Ia=9;
wm=564.9;
wm0=666.6;
istall=54;
Ra=0.22680;
La=0.28;
Kt=1.7629*10^-2;
Kb=1.7629*10^-2;
Jm=6.24*10^-4;
CM=2.909*10^-5;
ia0=1.1;
Trat=Ia*Kt;

%pid values

Kp=1;
Kl=1;
Kd=1;

%start position

eta_start_local = [1,0,0];

%requested position

eta_req_global = [0,0,0];


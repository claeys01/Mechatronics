%%
clear all
close all

%Data = readtable("Tito Neri - Speed to Force.xlsx");
%Port_Side_Thruster = Data(1:15,1:2);
%Starboard_Thruster = Data(19:33,1:2);
%Bowthruster = Data(38:55,1:2);



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


%Inputs
syms u v w r x y Psi tau_u tau_v tau_r; 
%u = ; %surge
%v = sym('v'); % sway
w = 0; % heave
p = 0; % roll
q = 0; % pitch
%r = sym('r');  % yaw

V =   [u;v;r];            % Velocity vector
eta = [x;y;Psi];          %Position vector ?
tau = [tau_u;tau_v;tau_r]; % Force input vector




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

% Coriolis and centripetal matrix

% Rigid-Body Coriolis and centripetal
C_rb = [0,0, -m*(x_g*r+v);
        0,0, -m*(y_g*r-u);
        m*(x_g*r + v), m*(y_g*r- u),0];

% Hydrodynamic and Centripetal matrix

a_1 =  X_dotu*u;
a_2 =  Y_dotv*v;

C_a = [0,   0,   a_2;
       0,   0,  -a_1;
      -a_2, a_1,  0];

C = C_rb + C_a;

% Linear Damping Matrix

%d11 = alpha_u*u^3 + beta_u*u^2 + gamma_u*u ;
%d22 = alpha_v*v^3 + beta_v*v^2 + gamma_v*v ;




%From file brightspace
d_11 = 0.0069*u^3 - 0.1314*u^2 + 0.8842*u - 2.0431;  %X_drag
d_22 = 0.0289*v^3 - 0.5191*v^2 + 3.3058*v - 2.3601;  %Y_drag
c1 = 122.3;
c2 = 4.8;
d_33 = (1/80)*c1*r^3*L^4 + (1/12)*c2*r*L^2;

D = [d_11, 0,   0;
      0,  d_22, 0;
      0,   0,  d_33];

%From excell
syms rpm
portsidethruster  = 2E-09*rpm^3 + 2E-07*rpm^2 + 0.0003*rpm + 0.0261;
starboardthruster = 1E-09*rpm^3 + 2E-07*rpm^2 + 0.0004*rpm - 0.0368;


R_t = -(4.7278*u^3 +6.4285*u^2+0.2249*u -0.0049);
R_tast= 1.1*R_t ;

T1=2; % fixed for now
T2=2;
Tb=0;
alpha1= 0 ; %recht naar achteren
alpha2= 0 ; %recht naar achteren
tau_u = T1*cos(alpha1)+T2*cos(alpha2) ;
tau_v = T1*sin(alpha1)+T2*sin(alpha2)+Tb ;
tau_r = -T1*cos(alpha1)*0.065+T2*cos(alpha1)*0.35+T2*cos(alpha2)-T2*sin(alpha2)+Tb*0.35;

a_c = 0; % angle one of the current in our case (2D) this stays 0
b_c = 0; % angle two of the current can be changed in 2D plane
speed_c = 0.2; % speed of the current
%v_c_NED = [speed_c*cos(a_c)*cos(b_c) speed_c*sin(b_c) speed_c*sin(a_c)*cos(b_c)];
%R = [cos(Psi) -sin(Psi) 0; sin(Psi) cos(Psi) 0; 0 0 1];
%v_c_BOD = v_c_NED * R;


tau = [tau_u;tau_v;tau_r];
%Environmental Force
tau_e = [0 0 0];
%tau_e = [0;0;0];

%AA = tau*M^(-1);

% hallo











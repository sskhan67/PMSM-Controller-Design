clear all
clc

%% Parameters
R=67.69e-3; % Stator inner  radius [m]
N=32;% Number of turns 
m=3; % phase
P_r=1; % pole pair
W=5;  %  Coil width 
T_p=6;%Pole pitch 
sp= 36; % Total stator slots
q=sp/(2*P_r); % Number of slots per pole per phase 
L=160.22e-3; % machine stack length[m]
n1=80;   % 1 slot= n1 divisions 
del_phi=(2*pi)/(sp*n1);  
u0=pi*4e-7;  % Permeability constant 
L_ind0=u0*R*L;
g_min=.47e-3; % Main airgap length
g_max=34.3e-3;% Intetr polar slot space  airgap length
%% Airgap Function 

angle_step=360;

 for  j=1:1:360
            
      phi1(j)=2*j*pi/angle_step;
 
      for  i=1:1:360
              g(j,i)=0;
         theta(i)=2*i*pi/360;
         g(j,i)=1/g_min+(1/g_max*cos(P_r*2*((phi1(j)-theta(i)))));
      end
         
 end 

Inv_Airgap=g(1,:);
figure(1), clf
plot(Inv_Airgap)

%
xlabel('Rotor position, \theta_r [deg]','FontSize',20);
ylabel('Inverse air gap function [1/m]','FontSize',20,'LineWidth',3);


grid on
%% Phase Winding Calculation    

angle_step = 360;
for k= 1:m
    for i=1:1:angle_step
         phi(i)=i*pi*2/angle_step;
         N_all(k,i) = 0;
         n=1;
%           for n=1:2:5 % No. of harmonics
     % k_wm= Winding factor
%  k_wm= (sin((n*pi)/(2*m))*sin(n*(pi/2)*(W/T_p)))/(q*sin((n*pi)/(2*q*m))); 
 k_wm= 1;
     % winding function
 N_all(k,i) = N_all(k,i) + 4*N/(P_r*2*pi*m)*(k_wm/n)*cos(P_r*n*(phi(i)-(k-1)*2*pi/(P_r*m))); 
%           end
    end
end
% Self winding function 
N_A=N_all(1,:); 
N_B=N_all(2,:);
N_C=N_all(3,:);
figure(2), clf
plot(N_all(1,:))
hold on
plot(N_all(2,:), 'r')
plot(N_all(3,:), 'g')
xlabel('Angular position, \phi,','FontSize',20);
ylabel('Self Winding function','FontSize',20,'LineWidth',3);
grid on
%Mtutual winding function 
N_AB=N_all(1,:).*N_all(2,:);
N_AC=N_all(1,:).*N_all(3,:);
N_BC=N_all(2,:).*N_all(3,:);
figure(3), clf
plot(N_AB, 'b')
hold on
plot(N_AC,'r')
plot(N_BC,'g')
xlabel('Angular position, \phi ','FontSize',20);
ylabel('Mutual Winding function','FontSize',20);
grid on
%% Indunctance Calculation[mH]
L_A=0;   
L_B=0;
L_C=0;
L_AB=0;
L_AC=0;
L_BC=0;
for l=1:1:(n1*sp)
    L_A=L_A+(N_A.^2*del_phi).*Inv_Airgap;  
    L_B=L_B+(N_B.^2*del_phi).*Inv_Airgap;
    L_C=L_C+(N_C.^2*del_phi).*Inv_Airgap;
    L_AB=L_AB+(N_AB*del_phi).*Inv_Airgap;
    L_AC= L_AC+(N_AC*del_phi).*Inv_Airgap;
    L_BC=L_BC+(N_BC*del_phi).*Inv_Airgap;
end
%%
%L_l=1; % Leakage inductance [mH];
C_p=60; % No of conductors per pole per phase
h=20;    % No. of parallel paths
t=15;    % No. of layers of stator windings
l_b=5;  % Mean length of end winding
phi_l=50;%No of flux lines 
B= 500; % Stator coil length
L_l=((2*P_r*C_p^2)/(h^2*1e8))*(l_b/t+phi_l*B)*1000;   %[mH]

%%
L_a=(L_ind0*L_A*1000)+L_l; % self inductance phase a
L_b=(L_ind0*L_B*1000)+L_l; % self inductance phase b
L_c=(L_ind0*L_C*1000)+L_l; % self inductance phase c
figure(4), clf
plot(L_a)
hold on
plot(L_b,'r')
plot(L_c,'g')
xlabel('Rotor angular position [degree] ','FontSize',20);
ylabel('Self inductance [mH]','FontSize',20);
grid on
L_ab=L_ind0*L_AB*1000; % mutual inductance phase ab
L_ac=L_ind0*L_AC*1000;
L_bc=L_ind0*L_BC*1000;
figure(5), clf
plot(L_ab, 'b')
hold on
plot(L_ac, 'g')
plot(L_bc, 'r')
xlabel('Angular position ','FontSize',20);
ylabel('Mutual phase inductance','FontSize',20);
grid on
%%%%%Phase  Inductance Matrix[mH]
L_a_b_c=[L_a L_ab L_ac;L_ab L_b L_bc; L_ac L_b L_c];

%% d-q inductance (for sinusiodal phase inductance)
x=max(L_a);          % Max value of Pahse 'a' self inductance 
y=min(L_a);          % Min value of Pahse 'a' self inductance 
L_AA=(x+y)/2;        % Constant Magnetizing inductance   
L_BB=L_AA-y;         %Amplitude of sinusoidal  varying  magnetizing  inductance 
L_d= (3/2*(L_AA+L_BB))/1000  % [H]
L_q=(3/2*(L_AA-L_BB))/1000   % [H]
figure(6),clf 
plot(theta,L_a) 
hold on  
plot(theta,L_AA,'r') 
grid on
%% d-q inductance using Transformation matrix[mH]
angle_step=360; 
for j=0:1:359
    i=j+1;
    theta=2*pi*i/angle_step;
T=[cos(theta) cos(theta-2*pi/3) cos(theta-4*pi/3); sin(theta) sin(theta-2*pi/3) sin(theta-4*pi/3);1/sqrt(2) 1/sqrt(2) 1/sqrt(2) ]; % transforation matrix
B=inv(T); % Inverse of T
L_a_b_c=[L_a(i) L_ab(i) L_ac(i); L_ab(i) L_b(i) L_bc(i); L_ac(i) L_bc(i) L_c(i)];
L_dq=T*L_a_b_c*B;    % Inductance unit [H]
% figure (100)
% plot(theta,L_dq(1,1))
% hold on
% figure (101)
% plot(theta,L_dq(2,2))
% hold on
end
%% Motor input parameters 
r_s=0.5; % Stator resistance [ohm]
J=.001; % Moment of inrertia
U_max=400;  % Dc volt [V]
I_s_max=100; % Max stator current [A]

%% Current controller design
D = .707; % Damping factor

Ts_iq = L_q/r_s;               % Time constant of the Id-Model [s]
w_0q = (2*pi)/(Ts_iq);         % Eigenfrequency in [1/s]
K_iiq = w_0q^2*L_q ;          % I-Parameter
K_ipq = (2*D*w_0q*L_q) - r_s ;% P-Parameter

% I_d Controller PI
Ts_id = L_d/r_s;               % Time constant of the Id-Model [s]
w_0d = (2*pi)/(Ts_id);         % Eigenfrequency in [1/s]
K_iid = w_0d^2*L_d ;          % I-Parameter
K_ipd = (2*D*w_0d*L_d) - r_s ;% P-Parameter
%% Speed controller design 
w0_w = 0.1/(Ts_iq/10);          % Eigenfrequency in [1/s]
alpha = 10;                     % Pole factor
K_wd = (w0_w*Ts_iq/10*(2*D + alpha) - 1);
K_wp = w0_w^2*Ts_iq/10*(2*alpha*D + 1);
K_wi = w0_w^3*Ts_iq/10*alpha;

   
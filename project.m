%1D unsteady heat conduction problem
clc
clear all;
%input

L=1;                          % wall length(m)
lambda=5.0;                           % conductivity(W/m-K)
rho=2000;                        % density(kg/m^3)
cp=200;                          % specific heat capacity(J/kg-K)
T_ini=293.2;                      % initial temperature(K)
t_fin=273.2;                      %final temperature
alpha=500;                        % heat transfer coefficient(K)

% Setup time steps
M=100;                            % number of time steps
t=40;
DELTA_t=t/(M);                 % time step duration(s)
for j=1:M
    time(j)=(j-1)*DELTA_t;
end

%Setup nodes
N=10;                              % number of nodes 
DELTA_x=L/(N);                   % distance between adjacent nodes 
x=0:DELTA_x:L;                   % position of each node 
%criteria
DELTA_t_crit_N = DELTA_x*rho*cp/(2*(lambda/DELTA_x+alpha));
if DELTA_t > DELTA_t_crit_N
    disp(' ')
    disp('Time step exeeds the limit');
    return
end

%Initial wall temperatures T(i,1)
for i=1:N+1
   T(i,1)=T_ini;
end
% Step trough time
for j=1:(M-1)
    % Heat flux condition(q=n*(-k*dT/dx))[W/m^2]
    T(1,j+1)=T(1,j)+2*lambda*(T(2,j)-T(1,j))*DELTA_t/(rho*cp*DELTA_x^2);
      for i=2:(N)
           T(i,j+1)=T(i,j)+lambda*(T(i-1,j)+T(i+1,j)-2*T(i,j))*DELTA_t/(rho*cp*DELTA_x^2);
      end
      % Heat flux condition(q=n*(-k*dT/dx))[W/m^2] + heat transfer coefficient(hout*(Tfin-T))[W/(m^2*K]
      T(N+1,j+1)=T(N,j)+(2*lambda*(T(N-1,j)-T(N,j))/(rho*cp*DELTA_x^2)+2*alpha*(t_fin-T(N,j))/(rho*cp*DELTA_x))*DELTA_t;
end
%plot
plot(time,T)

figure
plot(x,T)


figure
g=T(:,M);
plot(x,g)
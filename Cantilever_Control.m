% Dynamics and Control of Aerospace Structures with Fundamentals of
% Aeroelasticity
% Cantilevar Beam LQR Control
% Material: Aluminum
% Step tip load
% By Ahmhmd Polimi 2016
tic
clc
clear
%% Material Properties (Aluminum)
E=69e9; % Elastic modulus (Pa)
rho=2.7; % Density (g/cm^3)
rho=rho/(1000/(100)^3); % Density (kg/m^3)
f=1; % 1 N step load at the tip
tf=1.5; % Final time
tspan=linspace(0,tf,500); %The time for the simulation

%% Beam Properties
b=0.2; % Section width (m)
t=0.05; % Section thickness (m)
l=2; % Beam length (m)
A=b*t; % Section area (m^2)
J=(b*t^3)/12; % Second Moment of Inertia (m^4)
m_al=rho*A; % mass per unit length (kg/m)

%% Finite Elements Discretization
m=50; % Number of elements
n=4; % Number of DOFs
NN=m*2+2;
h=l/m; % Element length
xbeam=0:h:l; % Divisions of domain

% Element Matrices
% Save element matrices in a structure
Ke=((E*J)/h^3)*[12 6*h -12 6*h
                6*h 4*h^2 -6*h 2*h^2
                -12 -6*h 12 -6*h
                6*h 2*h^2 -6*h 4*h^2];
Me=((m_al*h)/420)*[156 22*h 54 -13*h
                22*h 4*h^2 13*h -3*h^2
                54 13*h 156 -22*h
                -13*h -3*h^2 -22*h 4*h^2];
fe=[0 0 f 0]';

% Matrix Assembly
M=zeros(NN,NN);
K=zeros(NN,NN);
F=zeros(NN,1);

for i=1:m
    Le=zeros(4,NN); % Matrix for contribution of element in the global matrix
    Le(1,2*i-1)=1;
    Le(2,2*i)=1;
    Le(3,2*i+1)=1;
    Le(4,2*i+2)=1;
    K=K+Le'*Ke*Le;
    M=M+Le'*Me*Le;
    if i==m
        F=F+Le'*fe;
    else
        F=F+Le'*zeros(n,1);
    end
end

%% Boundary conditions
% Eliminating first and 2nd rows and columns of both matrices
% Beam clamped at x = 0
K=K(3:NN,3:NN);
M=M(3:NN,3:NN);
F=F(3:NN);

%% Solving the eigenvalue problem to find
% the approximate modes for the motion
[U,omegas]=eig(K,M);
omega=sqrt(omegas);

% Modal matrices
M_m=U'*M*U; % Modal mass
K_m=U'*K*U; % Modal stiffness
F_m=U'*F; % Modal force

%% Damping matrix
% Damping model based on proportional damping
% C=alpha*M+Beta*k
% The coefficients selected in order to give damping ratio of 0.005 for the
% first and the 2nd modes
xi_1=0.005;
xi_2=0.005;
betac=(2*xi_2*omega(2,2)-2*xi_1*omega(1,1))/(omega(2,2)^2-omega(1,1)^2);
alfac=2*xi_1*omega(1,1)-betac*omega(1,1)^2;
C=alfac*M+betac*K;
C_m=alfac*M_m+betac*K_m;
% Building a vector for xi
nx=NN-2;
xi=zeros(1,nx);
omegad=zeros(1,nx);
for i=1:nx
    xi(i)=C_m(i,i)/(2*M_m(i,i)*omega(i,i));
    omegad(i)=omega(i,i)*sqrt(1-xi(i)^2);
end

% State space (modal coordinates)
Lc=zeros(nx,1);
Lc(nx-1)=1;
A_sm=zeros(2*nx,2*nx);
B_sm=zeros(2*nx,1);
C_sm=zeros(1,2*nx);
C_sm(1)=1;
for i=1:nx
    A_sm(2*i-1:2*i,2*i-1:2*i)=[0 1
        -omega(i,i)^2 -2*xi(i)*omega(i,i)];
    B_sm(2*i-1:2*i,1)=[0 U(:,i)'*Lc/M_m(i,i)]';
    C_sm(1,2*i-1:2*i)=[U(end-1,i) 0];
end
Cantsys_m=ss(A_sm,B_sm,C_sm,0);

% Reduced order model
% Only taking the first n_md modes
n_md=1; % Number of modes to be included in the reudce modeal
A_r=A_sm(1:2*n_md,1:2*n_md);
B_r=B_sm(1:2*n_md,1);
C_r=C_sm(1,1:2*n_md);
Cantsys_r=ss(A_r,B_r,C_r,0);

%% LQR Controller Design
% Designing a controller to reduce the total energy
Wz=1;
Wu=1;
rho=0.01;
% The matrix H for obtaining the total energy
H=zeros(2*n_md,2*n_md);
for i=1:n_md,
H(2*i-1:2*i,2*i-1:2*i) = [omega(i,i)/sqrt(2) 0; 0 1/sqrt(2)];
end
% Finding the matrices Q and R
Q=H'*Wz*H;
R=rho*Wu;
G=lqr(A_r,B_r,Q,R);
Cantsys_con=ss(A_r-B_r*G,B_r,C_r,0);
[Y,T]=step(Cantsys_m,tspan);
[Yc,Tc]=step(Cantsys_con,tspan);
figure(1)
plot(T,Y,Tc,Yc);
title('Tip displacement of the beam')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
grid on
legend('Without LQR',['With LQR \rho= ',num2str(rho)])

%% LQR Controller + Observer
rhoo=0.01;
Qo=eye(2*n_md);
Ro=rhoo;
LT=lqr(A_r',C_r',Qo,Ro);
Lobs=LT';
Observer=ss(A_r-B_r*G-Lobs*C_r,Lobs,-G,0);
Cantsys_obs=feedback(Cantsys_m,Observer,+1);

%% Simulation of controlled system
xo=[0 1 zeros(1,2*(n_md-1))]';
[y,tn,x]=initial(Cantsys_m,[xo;zeros(2*(nx-n_md),1)],tspan);
[yc,tc,xc]=initial(Cantsys_con,xo,tspan);
[yo,to,xo]=initial(Cantsys_obs,[xo;zeros(2*nx,1)],tspan);
E=zeros(1,length(tn));
u=zeros(1,length(tn));
Ec=zeros(1,length(tc));
uc=zeros(1,length(tc));
Eo=zeros(1,length(to));
uo=zeros(1,length(to));
for i=1:length(tn)
    E(i)=x(i,1:2*n_md)*H^2*x(i,1:2*n_md)';
    u(i)=-G*x(i,1:2*n_md)';
end
for i=1:length(tc)
    Ec(i)=xc(i,:)*H^2*xc(i,:)';
    uc(i)=-G*xc(i,:)';
end
for i=1:length(to)
    Eo(i)=xo(i,1:2*n_md)*H^2*xo(i,1:2*n_md)';
    uo(i)=-G*xo(i,1:2*n_md)';
end
figure(2)
subplot(2,1,1)
plot(tn,E,tc,Ec,to,Eo)
title('Total Energy vs time')
xlabel('Time (seconds)')
ylabel('Total Energy Nm')
legend('W/o control',['+LQR \rho= ',num2str(rho)],['+Observer \rho= ',num2str(rhoo)])
grid on
subplot(2,1,2)
plot(tc,uc,to,uo)
title('Control force vs time')
xlabel('Time (seconds)')
ylabel('Control Force N')
legend('LQR Controller','Controller+Observer')
grid on
function [tsd,xsd]=ADM4DSDE(tend,dt,ICOND)
% clc; 
% clear all;
% 
% set(groot,'defaulttextinterpreter','latex');  
% set(groot, 'defaultAxesTickLabelInterpreter','latex');  
% set(groot, 'defaultLegendInterpreter','latex');
% set(groot,'defaultAxesFontSize',18)
% set(groot,'defaultLegendFontSize',14)
% set(groot, 'DefaultLineLineWidth', 4);
%% From paper default values

% D=0.15;
% S1in=50;
% y1=42.14;
% y2=116.5;
% y3=268;
% y4=1.165;
% d1=0.01;
% d2=0.01;
% kappa1=1.2; %kappa in meadows
% H1= 0.1; %(r1  in meadows et al.)
% kappa2=1.64; %m_I in meadows
% K=9.28;
% k1=0.05;
% k2=0.5;
%% 


D=0.015;
S1in=500;
y1=42.14;
y2=116.5;
y3=268;
y4=1.165;
k1=0.1;
k2=0.001;
mu11=120; %kappa in meadows
H1= 10; %(r1  in meadows et al.)
mu22=0.064; %m_I in meadows
H2=9.28;

%% Take3
global alpha beta sigma1
a=1/H1;
b=1/H2;
c=y2*b;
d=y3*mu22*b/mu11;
%  alpha=S1in/H1*D/mu11;
 alpha=4;
 beta=y1*a/c;
% sigma1=k1/mu11;
sigma1=.1;
sigma2=k2/D;
%      epsilon=D/mu11
    epsilon= 3.00e-02;
% omega=mu22/D;
omega=4;


% 
% alpha
%   alpha=.001;
% sigma1
% sigma1=1.1;
 
 if (sigma1>1)||(ICOND(4)<alpha/beta)
     fprintf('error in sigma1 or beta alpha y0 values')
     sigma1
     y0diff=ICOND(4)-alpha/beta
%      return
 end


NDf=@(t,x) [alpha-epsilon*x(1)-beta*x(1)*x(3)/(1+x(1));
    -epsilon*x(2)+x(1)*x(3)/(1+x(1))-x(2)*x(4)/(1+x(2));
    -epsilon*x(3)-sigma1*x(3)+x(1)*x(3)/(1+x(1));
    -epsilon*(1+sigma2)*x(4)+epsilon*omega*x(4)*x(2)/(1+x(2))]; %this is the ODE part

Noise=@(t,x) [0.1*x(1);
    0.1*x(2);
    0.1*x(3);
    0.1*x(4)];

% Noise=@(t,x) 10*[0.1;
%     0.1;
%     0.1;
%     0.01];
% dt=.01;
% tend=50
 tsd=0:dt:tend;
opts=sdeset('SDEType','Ito');%'Randseed',1);
 xsd=sde_euler(NDf,Noise,tsd,ICOND,opts);




%% plots of both dim and non dim model with all params
% 
%    opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
% % % % [t,x]=ode45(f,[0,200],Icond,opts);
% % 
%    [t,x]=ode45(NDf,[0,tend],NDicond,opts);
% % % 
% % %  
%  figure(3)
% subplot(4,1,1);
% plot(t,x(:,1),tsd,xsd(:,1));
% ylabel('$S_1$')
% subplot(4,1,2);
% plot(t,x(:,2),tsd,xsd(:,2));
% ylabel('$S_2$')
% 
% subplot(4,1,3);
% plot(t,x(:,3),tsd,xsd(:,3));
% ylabel('$X_1$')
% subplot(4,1,4);
% plot(t,x(:,4),tsd,xsd(:,4));
% ylabel('$X_2$')
% xlabel('Time')
end
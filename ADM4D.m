function [t,x]=ADM4D(tend,ICOND,al)

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

switch nargin
    case 3
        alpha=al;
    otherwise
        alpha=6;
end
alpha
 beta=y1*a/c
% sigma1=k1/mu11;
sigma1=.1
sigma2=k2/D
%      epsilon=D/mu11
%     epsilon= 3.00e-02;
epsilon=0.000125;
% 
% beta=0.3357;
% sigma1=0.1;
% sigma2=0.66667;
omega=4.266;
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


    opts = odeset('RelTol',1e-5,'AbsTol',1e-5);

    [t,x]=ode45(NDf,[0,tend],ICOND,opts);



end
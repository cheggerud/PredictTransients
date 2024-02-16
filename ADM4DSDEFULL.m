function [tsd,xsd]=ADM4DSDEFULL(tend,dt,ICOND,sig)
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



% D=0.00015;
% S1in=500000;
D=0.055;
S1in=700;
y1=42.14;
y2=116.5;
y3=268;
y4=1.165;
k1=0.1;
k2=0.001;
mu11=0.5;%1.2; %kappa in meadows
H1= 10; %(r1  in meadows et al.)
mu22=0.0064; %m_I in meadows
H2=9.28;


df=@(t,x) [D*(S1in-x(1))-y1*mu11*x(1)*x(3)/(H1+x(1));  %x1 is the simple substrate, x2 the acids, x3 the acidogens, x4 the methanogens. 
    -D*x(2)+y2*mu11*x(1)*x(3)/(H1+x(1))-y3*mu22*x(2)*x(4)/(H2+x(2));
    -D*x(3)-k1*x(3)+mu11*x(1)*x(3)/(H1+x(1));
    -D*x(4)-k2*x(4)+mu22*x(4)*x(2)/(H2+x(2))]; %this is the ODE part



Noise=@(t,x) [sig*x(1);
    sig*x(2);
    sig*x(3);
    sig*x(4)];


 tsd=0:dt:tend;
% opts=sdeset('SDEType','Ito')%,'Randseed',1);
%  xsd=sde_euler(df,Noise,tsd,ICOND,opts);


xsd(1,:)=ICOND;
for i=2:length(tsd)
    xsd(i,:)=xsd(i-1,:)+(dt*df(0,xsd(i-1,:))+dt*(randn(4,1).*Noise(0,xsd(i-1,:))))';

    if xsd(i,1)<0
        xsd(i,1)=0;
    elseif xsd(i,2)<0
        xsd(i,2)=0;
         elseif xsd(i,3)<0
        xsd(i,3)=0;
         elseif xsd(i,3)<0
        xsd(i,3)=0;
         elseif xsd(i,4)<0
        xsd(i,4)=0;
    end
       
end



end
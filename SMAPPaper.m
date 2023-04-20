clc; clear all;


%% preamble

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


%% 



TEND=15;
ICOND= [0.800    0.500    0.500   40.0000];
alpha=4;

%generate SDE timeseries
umin=0;
umax=1;
vmin=0;
vmax=1;
xmin=0;
xmax=1;
ymin=alpha/beta;
ymax=100;

 ustar=sigma1/(1-sigma1)
tic



%% Set up the SMAP paramaters: We are training the data on a FULL time series. Just on one single time series. 




numhist=20;% if using historical data, how many points to use. 
dt=0.01; %how fine of resolution in the SDE realization:
dthistdat=.1; %how much time between data points in the SMAP;
dtdat=dthistdat/dt;
warning('off')
numrels=5000; %the number of randIC stochasitc realization we test. 
predind=1;% which index do we want to predict. 
theta=5;

%make the control time series. 
[tsd,xsd]=ADM4DSDE(TEND,dt,ICOND);
rawdata=xsd(1:dtdat:end,:);
traintime=tsd(1:dtdat:end);


controldata=rawdata(2:end,:);
controltime=traintime(2:end); 


%% Random ICs:

ICOND=[ICOND;rand(numrels-1,4)];
ICOND(2:end,1)=umin+(umax-umin)*ICOND(2:end,1);
ICOND(2:end,2)=vmin+(vmax-vmin)*ICOND(2:end,2);
ICOND(2:end,3)=xmin+(xmax-xmin)*ICOND(2:end,3);
ICOND(2:end,4)=ymin+(ymax-ymin)*ICOND(2:end,4);
tic
 for k=1:numrels
     %% create the SD full data:

     Icond=ICOND(k,:);
[tsd,xsd]=ADM4DSDE(TEND,dt,Icond);
%  xsd=log(xsd);

%% Determine at what time approx the transients end:

for i=10:length(tsd)-1
    nnn(i+1)=norm(xsd(i+1,:)-xsd(i-1,:))/(2*dt);
%     if abs((xsd(i+1)-xsd(i-1)))/(2*dt)<1e-1&& abs(xsd(i,1)-ustar)<0.5%abs((xsd(i+1)-xsd(i-1)))/(2*dt)<0.5&& abs(xsd(i,1)-log(ustar))<1 % if no scaling of xsd then use abs((xsd(i+1)-xsd(i-1)))/(2*dt)<1e-1&& abs(xsd(i,1)-ustar)<0.5
%         % istrans(i)=0;
  if norm(xsd(i+1,1)-xsd(i-1,1))/(2*dt)<1e-1 && abs(xsd(i,1)-ustar)<0.5
    ttrans(k)=tsd(i);
    i;
    break;
    
    else
        ttrans(k)=NaN;
    end
end

%  figure(1)   
%  xlabel('time')
%  ylabel('u(t)')
%  p(k)=plot(tsd,(xsd(:,1)),'linewidth',1);
%  hold on
%  q=plot(ttrans(k),(xsd(i,1)),'Color',p(k).Color);
%  q.Marker='o';
%  q.MarkerSize=10;
%  q.MarkerFaceColor=p(k).Color;
%  hold on

 

%% Take a subset of the full data (or not) and create the embeddings if neccessary:

rawdata=xsd(1:dtdat:end,:);
traintime=tsd(1:dtdat:end);

   
data=rawdata(2:end,:);

datatime=traintime(2:end); %now thats fixed. 

if k==1 %initialize vectors
    pred=NaN*ones(numrels,length(datatime)-1);
    eavg=NaN*ones(numrels,length(datatime)-1);
    prederror=NaN*ones(numrels,length(datatime)-1);
end
%%  Now insert the SMAP prediction.

%construct A and B;
%% If using control series use:

% for tstar=1:length(datatime)-1 %i is what is to be predicted. 
%     
%   
%     
%        normss=vecnorm((controldata-data(tstar,:)),2,2);
%         DM=1/length(normss)*sum(normss);
%      B=[];
%         A=[];
%         
%           for i=1:length(controltime)-1
%              B(i)=exp(-theta*normss(i)/DM)*controldata(i+1,predind);
%                
%                  for j=1:4 % the embedding dimension
%                      A(i,j)=exp(-theta*normss(i)/DM)*controldata(i,j);
%                  end
%          end
%         C=A\B';
%     
%     pred(k,tstar)=data(tstar,:)*C; %this is actually the forcast. SO the time index is one step behind here. 
%     prederror(k,tstar)=abs(pred(k,tstar)-data(tstar+1,predind)); %prediction error
%     prederrorrel(k,tstar)=abs(pred(k,tstar)./data(tstar+1,predind)-1);% $relative prediciton error
    
    %% If using the historical data then use for tind=numhist+1:length(traintime) %i is what is to be predicted. 
 for tind=numhist+1:length(datatime)-1   
    histdata=data(tind-numhist:tind-1,:);
    
       normss=vecnorm(histdata(1:end,:)-data(tind-1,:),2,2);
        DM=1/length(normss)*sum(normss);
     B=[];
        A=[];
        
          for i=1:numhist-1
             B(i)=exp(-theta*normss(i)/DM)*histdata(i+1,predind);
               
                 for j=1:4
                     A(i,j)=exp(-theta*normss(i)/DM)*histdata(i,j);
                 end
         end
        C=A\B';
    
    pred(k,tind)=data(tind-1,:)*C;
    prederror(k,tind)=abs(pred(k,tind)-data(tind,predind));
    prederrorrel(k,tind)=abs(pred(k,tind)./data(tind,predind)-1);%abs(pred(k,tind)-traindata(tind,predind));
    

end




 end
 toc
 %% find transient time from errors:
 %% Using the errors see if we can predict when the transients stop:
%start at the point when predictions were being made: DUH! 

%% USE the error to make the prediction for tstar: 

%smooth out the error curve in some fashion. 

smootherror=smoothdata(prederror,2,'gaussian',[10,0]);

smootherrorrel=smoothdata(prederrorrel,2,'gaussian',[10,0]);

%% Want to look at the first difference being over some threshold for both. 
 etols=logspace(-5,1,2000);
% 
differr=diff(smootherror,1,2);
diffrel=diff(smootherrorrel,1,2);
timediff=datatime(3:end);

%% predict the transient end for first differences in error

% for i=1:length(etols)
%     for k=1:numrels
%         ttranserror(k,i)=NaN;
%         for j=50:length(timediff)
%            
%             if differr(k,j)>etols(i)
%                 ttranserror(k,i)=timediff(j);
%                 break;
% %             elseif differr(k,j)>difmax
% %                 difmax=differr(k,j);
% %                 ttranserror(k,i)=timediff(j);
%             end
%         end
%     end
% end
% prop=ttranserror./ttranserror;
% prop(isnan(prop))=0;
% propserr=sum(prop,1)/numrels;
% 
% corrs=corrcoef([ttrans' ttranserror],'Rows','pairwise');
% figure(3)
% yyaxis left
% semilogx(etols,corrs(1,2:end));
% yyaxis right
% semilogx(etols,propserr);
% title('for using the first differences on abs error')
% 
% %below makes a correlation plot
% % bestcorr=find(max(corrs(1,2:end))==corrs(1,2:end));
% % 
% % figure(13)
% % clf
% % plot(5:15,5:15)
% % hold on
% %  str="corr coef="+corrs(1,bestcorr(1)+1) ;
% % text(10,7,str)
% % plot(ttrans,ttranserror(:,bestcorr(1)),'.','markersize',20)
% % xlabel('Numeric')
% 
% % 
% % %% predict the transient end for rel error
% % 
% for i=1:length(etols)
%     for k=1:numrels
%         ttranserrorrel(k,i)=nan;
%         difmax=0;
%         for j=50:length(timediff)
%            
%             if diffrel(k,j)>etols(i)
%                 ttranserrorrel(k,i)=timediff(j);
%                 break;
% %             elseif diffrel(k,j)>difmax
% %                 difmax=diffrel(k,j);
% %                 ttranserrorrel(k,i)=timediff(j);
%             end
%         end
%     end
% end
% prop=ttranserrorrel./ttranserrorrel;
% prop(isnan(prop))=0;
% propsrel=sum(prop,1)/numrels;
% 
% corrsrel=corrcoef([ttrans' ttranserrorrel],'Rows','pairwise');
% figure(4)
% 
% yyaxis left
% semilogx(etols,corrsrel(1,2:end));
% yyaxis right
% semilogx(etols,propsrel);
% title('for using the first differences on abs rel error')
% 
% % bestcorrrel=find(max(corrsrel(1,2:end))==corrsrel(1,2:end));
% % 
% % figure(14)
% % clf
% % plot(5:15,5:15)
% % hold on
% %  str="corr coef="+corrsrel(1,bestcorrrel(1)+1) ;
% % text(10,7,str)
% % plot(ttrans,ttranserrorrel(:,bestcorrrel(1)),'.','markersize',20) %plots the best correlation plot
% % xlabel('Numeric')
% 
% %% if using just the error 
% figure(9)
% clf
% % for k=1:numrels
% % plot(traintime,smootherror(k,:),'Color',p(k).Color)
% % hold on
% % end
% %  etols=logspace(-3,0,1000);
% for i=1:length(etols)
%     for k=1:numrels
%     for j=1:length(datatime)-1
%        
%         if smootherror(k,j)<etols(i)
%             ttranserror(k,i)=traintime(j);
%             if traintime(j)>5 % the transient must be at least 5 time units long, otherwise the low error comes from the luck of having close initial conditions. 
%             break;
%             end
%         else
%             ttranserror(k,i)=NaN;
%         end
%     end
%     end
% end
% 
% corrs=corrcoef([ttrans' ttranserror],'Rows','pairwise');
% semilogx(etols,corrs(1,2:end))
% title('using  abs error under a threshold ') %initially correlation will be good due to initial conditions starting closish to each other. 


figure(10)
clf
%% Using the metric for error is above a threshold and is a max. 
for i=1:length(etols)
    for k=1:numrels
    for j=1:length(datatime)-2
       
        if smootherrorrel(k,j)>etols(i)& smootherrorrel(k,j+1)<smootherrorrel(k,j)
            ttranserrorrel(k,i)=traintime(j);
            if traintime(j)>5 % the transient must be at least 5 time units long, otherwise the low error comes from the luck of having close initial conditions. 
            break;
            end
        else
            ttranserrorrel(k,i)=NaN;
        end
    end
    end
end
prop=ttranserrorrel./ttranserrorrel;
prop(isnan(prop))=0;
propsrel=sum(prop,1)/numrels;


corrsrel=corrcoef([ttrans' ttranserrorrel],'Rows','pairwise');
yyaxis left
semilogx(etols,corrsrel(1,2:end));
ylabel('correlation coefficient')
yyaxis right
semilogx(etols,propsrel);
ylabel('Proportion of predictions made')
xlim([1e-3,10])
xlabel('Error threshold, $\delta_m$')
%initially correlation will be good due to initial conditions starting closish to each other. 



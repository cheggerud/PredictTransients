clc; clear all;


%% preamble


D=0.055;
S1in=700;
% D=0.00015;
% S1in=50000;
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


figure(22)
clf

df=@(t,x) [D*(S1in-x(1))-y1*mu11*x(1)*x(3)/(H1+x(1));
    -D*x(2)+y2*mu11*x(1)*x(3)/(H1+x(1))-y3*mu22*x(2)*x(4)/(H2+x(2));
    -D*x(3)-k1*x(3)+mu11*x(1)*x(3)/(H1+x(1));
    -D*x(4)-k2*x(4)+mu22*x(4)*x(2)/(H2+x(2))];


%% Find the equilibrium Values
sig=5;
%time step of the SDE
TEND=40;% end time of the SDE
 ICOND= [.5 .464 .03, .675];
% ICOND=[1,1,1,100];

opts = odeset('RelTol',1e-15,'AbsTol',1e-15);
 [t,x]=ode45(df,[0,TEND],ICOND,opts);
% [t,x]=ADM4DSDEFULL(TEND,0.01,ICOND,0);
   [tsd,xsd]=ADM4DSDEFULL(TEND,0.01,ICOND,sig);% SDE simulation ensuring jumps won't go into the negative. If they do, then the bacteria stay at zero. 

% figure(3)
% plot(t,x(:,1),tsd,xsd(:,1))
% hold on
eqs=fsolve(@(x) df(0,x),[1 1 17 672]);
eqs=fsolve(@(x) df(0,x),eqs);% eqs now gives the equilibrium value. 




%% generate SDE timeseries:
numrels= 10;
dt=0.01;
 sigs=5*logspace(-4,0,10);
 % sigs=[0.1];
%this is the vector to play with the variance: % max sig should 5. Anything higher then no transients appear. i.e. more than 10% of sims give all transient dynamics. 
 
 % numtsteps=[10,20,100,200];%integer;
numtsteps=[5,10,20,100];
timestepbetweendata=numtsteps*dt;

 etols=logspace(-2,2,15);
 etolsl=length(etols);

sigl=length(sigs);
gapl=length(numtsteps);


theta=5;
%generate SDE timeseries

xmin=0;
xmax=0.5;
ymin=0;
ymax=1;
S1min=0;
S1max=1;
S2min=0;
S2max=1; % as discussed in the simplification part of the paper. We want our IC's to be small or else we know that transients might not occur. Thus we keep the IC somewhat restricted but still random. 

 ICOND1= [.5 .464 .03, .675];
  % rng(515)
ICOND=[rand(numrels,4)];
ICOND(2:end,1)=S1min+(S1max-S1min)*ICOND(2:end,1);
ICOND(2:end,2)=S2min+(S2max-S2min)*ICOND(2:end,2);
ICOND(2:end,3)=xmin+(xmax-xmin)*ICOND(2:end,3);
ICOND(2:end,4)=ymin+(ymax-ymin)*ICOND(2:end,4);

 difs=ICOND-eqs;
  difsq=sqrt(sum(difs.^2,2));% this gives the distance each IC is away from the EQ. Lets sort it just so its easier later:

  [difsort,I]=sort(difsq);
  difsq=difsq(I);
ICOND=ICOND(I,:);
 % error=NaN*ones([length(sigs),numrels,length(numtsteps),length(numhist),length(tsd)]);
 tenderror=NaN*ones([sigl,numrels,gapl,etolsl]);
   tenderroro=NaN*ones([sigl,numrels,gapl,etolsl]);
 tranend=NaN*ones([numrels,sigl]);





%% Test out the predictions for parameter space

for s=8
      figure(s)
      clf
for j=1:numrels
    clear ts xs
       [ts,xs]=ADM4DSDEFULL(TEND,0.01,ICOND(j,:),sigs(s));
S1=xs(:,1)';
S2=xs(:,2)';
X1=xs(:,3)';
X2=xs(:,4)';


pp=plot(tsd,S1)
 hold on
 drawnow
ssdif=diff(S1);
[m,in]=max(S1);
eqdif=abs(S1-eqs(1));
    for i=in:length(tsd)
                 if eqdif(i)<sigs(s)+5 %5*sqrt(sigs(s))
                     ind=i;
                tranend(j,s)=tsd(i);
                 break
                 else
                    tranend(j,s)=NaN;
                 end
    end   
    plot(tranend(j,s),S1(ind),'o','Color',pp.Color,'MarkerFaceColor',pp.Color,'MarkerSize',15)
    hold on
    xlabel('time')
    ylabel('$S_1(t)$','Interpreter','latex')
end
end


%% 
fprintf("generating SDES ")
warning('off')
% rng(1)
tic
for s=1:length(sigs)
    s
    [ttrain,xtrain]=ADM4DSDEFULL(TEND,0.01,ICOND(4,:),sigs(s));
   
 %       figure(100*s)
 %     clf
 % plot(ttrain,xtrain(:,1))
% hold on
% concat the training data to end 5 units after the transient ends: 
Stra1=xtrain(:,1)';

[m,in]=max(Stra1);
eqdif=abs(Stra1-eqs(1));
    for i=in:length(ttrain)
                 if eqdif(i)<sigs(s)+2 %5*sqrt(sigs(s))
                     tendtrainind=i;
                tendtrain=ttrain(i);
                 break
                 else
                    tendtrain=NaN;
                 end
    end 
    if isnan(tendtrain)
        error('No transient computed in contral data')
    end
% iterate over all other rels:


for j=1:numrels
    clear ts xs
       [ts,xs]=ADM4DSDEFULL(TEND,0.01,ICOND(j,:),sigs(s));
S1=xs(:,1)';
S2=xs(:,2)';
X1=xs(:,3)';
X2=xs(:,4)';

 % figure(s)
 % plot(tsd,S1)
 % hold on
 % drawnow
ssdif=diff(S1);
[m,in]=max(S1);
eqdif=abs(S1-eqs(1));
    for i=in:length(tsd)
                 if eqdif(i)<sigs(s)+2 %5*sqrt(sigs(s))
                     ind=i;
                tranend(j,s)=tsd(i);
                 break
                 else
                    tranend(j,s)=NaN;
                 end
    end   
    % plot(tranend(j,s),S1(ind),'o','MarkerSize',20)
    % hold on

    for gap=1:gapl 
            tsdat=ts(1:numtsteps(gap):end);
            s1dat=S1(1:numtsteps(gap):end);
            s2dat=S2(1:numtsteps(gap):end);
            x1dat=X1(1:numtsteps(gap):end);
            x2dat=X2(1:numtsteps(gap):end);

             tstrain=ttrain(1:numtsteps(gap):tendtrainind+500);
             s1train=xtrain(1:numtsteps(gap):tendtrainind+500,1)';
             s2train=xtrain(1:numtsteps(gap):tendtrainind+500,2)';
             x1train=xtrain(1:numtsteps(gap):tendtrainind+500,3)';
             x2train=xtrain(1:numtsteps(gap):tendtrainind+500,4)';
            
                
                % clear err
                 err=nan*ones(length(tsdat),1);
                 
                   for i=1:length(s1dat)-1% this is now the SMAP loop.
                    A=[];
                    B=[];
                   
                    
                    norm= sqrt((s1train-s1dat(i)).^2+(s2train-s2dat(i)).^2+(x1train-x1dat(i)).^2+(x2train-x2dat(i)).^2);  %compute the weights: distance of each training point to the current point. 
                    DM=1/length(norm)*sum(norm);
                    for k=1:length(s1train)-1 % This loop is to populate the SMAP;
                        B(k)=exp(-theta*norm(k)/DM)*s1train(k+1);
                        A(k,:)=exp(-theta*norm(k)/DM)*[s1train(k),s2train(k),x1train(k),x2train(k)];
                        
                    end
                    C=A\B';
                    proj=[s1dat(i),s2dat(i),x1dat(i),x2dat(i)]*C;
                     %now to plot the predictions they start at tsd(1:gap:end), then tsd(numhist:end);
                    % err(i)=abs(proj-s1dat(i+1));
                    err(i)=abs((proj-s1dat(i+1))/s1dat(i+1));
                       % error(s,j,gap,his,(i)*numtsteps(gap))=err(i);
                    % relerror(s,n,gap,his,(i)*numtsteps(gap))=relerr;
                    % %toggle this for relative error:
                   end
                   %do something with the error curve now:
                     smootherror=smoothdata(err,'gaussian',[4,0]);
                    err=smootherror;


                     % if s==2
                     %                 figure(gap*33)
                     %                 pp=plot(tsdat,err);
                     %                 hold on
                     %                 % plot(trendcheck,err(i),'*','MarkerSize',10,'Color',pp.Color);
                     %                 % hold on
                     %                 plot(tranend(j,s),0,'o','MarkerSize',20,'Color',pp.Color);
                     %                 hold on
                     %                 plot(tsdat,s1dat,':','Color',pp.Color)
                     %                 hold on
                     %                 drawnow
                     %        end

                   % figure(2)
                   % plot(tsdat,err)
                   % hold on
                   % 

                    locmaxerr=islocalmax(err);
                   for e=1:etolsl 
                       
                       for i=1:length(err)
                           if err(i)>etols(e)&locmaxerr(i)
                               trendcheck=tsdat(i);
                             
                               if trendcheck>6 % the transient must be at least 5 time units long, otherwise the low error comes from the luck of having close initial conditions. 
                             tenderror(s,j,gap,e)=trendcheck;

                            % if s==9
                            %          figure(22)
                            %          pp=plot(tsdat,err);
                            %          hold on
                            %          plot(trendcheck,err(i),'*','MarkerSize',10,'Color',pp.Color);
                            %          hold on
                            %          plot(tranend(j,s),0,'o','MarkerSize',20,'Color',pp.Color);
                            %          hold on
                            % end
                                break;
                               end
                           end
                            
                       end



                       %  sp=reshape(error(s,j,gap,his,:),[1],[]);                
                       % for i=1:length(sp)
                       %      if sp(i)<etols(e)
                       %        trendcheck=tsd(i);
                       %      if trendcheck>4 % the transient must be at least 5 time units long, otherwise the low error comes from the luck of having close initial conditions. 
                       %       tenderroro(s,j,gap,his,e)=trendcheck;
                       %          break;
                       %      end
                       % 
                       % 
                       %   end
                       % end
                   end

            end
            
    end



end
toc







 






%% Now try to compute the corrcoefs?" or r2;

%tranend(n,s) is to be compared to all of tenderror(s,n,:,;,:);

correrror=zeros(([length(sigs),length(numtsteps),length(etols)]));
correrrorcheck=zeros(([length(sigs),length(numtsteps),length(etols)]));
rsq=zeros(([length(sigs),length(numtsteps),length(etols)]));
%try plotting to start?

for s=1:sigl
    for gap=1:gapl
        
            for e=1:length(etols)
                correrror(s,gap,e)=corr(tranend(:,s),tenderror(s,:,gap,e)','rows','complete');
                %compute it by hand and compare?
                cc=[tranend(:,s),tenderror(s,:,gap,e)'];
                 cc(any(isnan(cc), 2), :) = [];
                 xxs=cc(:,1);
                 yys=cc(:,2);
                 xbar=mean(xxs);
                 ybar=mean(yys);
                 num=sum((xxs-xbar).*(yys-ybar));
                 den=sqrt(sum((xxs-xbar).^2)*sum((yys-ybar).^2));

                 correrrorcheck(s,gap,e)=num/den;
r=1-sum((xxs-yys).^2)./sum((yys-ybar).^2);
                   rsq(s,gap,e)=r;

                 check=correrror(s,gap,e)-correrrorcheck(s,gap,e);
                  if check>0.001                 
                      check
                      break
                  end
                  if isnan(correrror(s,gap,e))
                      if isnan(correrrorcheck(s,gap,e))
                          
                      else 
                          fprintf('error not both nan')
                      end
                  end
            
        end
    end
end
%% PLot the corr as a heat maps:
 figure(102)
        clf
       
    t0=tiledlayout('flow')
     figure(104)
        clf
    t1=tiledlayout('flow')
    correrror(isnan(correrror))=0;
    rsq(isnan(rsq))=0;



    for gap=1:gapl
        
        ebysig=reshape(rsq(:,gap,:),[sigl],[length(etols)]);
figure(102)
h=nexttile
pcolor(etols,sigs,ebysig);
 xticks([1e-2,1e1])
  yticks([1e-3,1e0])
set(gca, 'XScale', 'log');
 set(gca, 'YScale', 'log');
      % shading interp
title("$\Delta t=$"+timestepbetweendata(gap),'Interpreter','latex','FontSize',10);

 caxis manual
 caxis([0 1]);
hold on


  ebysig=reshape(correrror(:,gap,:),[sigl],[length(etols)]);
figure(104)
h1=nexttile
pcolor(etols,sigs,ebysig);
 xticks([1e-2,1e1])
  yticks([1e-3,1e0])
set(gca, 'XScale', 'log');
 set(gca, 'YScale', 'log');
      % shading interp
title("$\Delta t=$"+timestepbetweendata(gap),'Interpreter','latex','FontSize',10);

 caxis manual
 caxis([0 1]);
hold on


    end


xlabel(t0,'Error threshold, $\delta_r$','Interpreter','latex','FontSize',14)
ylabel(t0,'Noise level, $\sigma$','Interpreter','latex','FontSize',14)
xlabel(t1,'Error threshold, $\delta_r$','Interpreter','latex','FontSize',14)
ylabel(t1,'Noise level, $\sigma$','Interpreter','latex','FontSize',14)

cbh = colorbar(h(end)); 
% To position the colorbar as a global colorbar representing
% all tiles, 
cbh.Layout.Tile = 'east'; 

cbh = colorbar(h1(end)); 
% To position the colorbar as a global colorbar representing
% all tiles, 
cbh.Layout.Tile = 'east'; 

%% Plot the corr plots below: Is seems fitting?






%% look at proportion of predictions made. i.e top left should be blue blue. 

succ=zeros([sigl,gapl,etolsl]);

succs=tenderror;
succs(isnan(succs))=0;
succs=succs./succs;
succs(isnan(succs))=0;

succ=sum(succs,2)/numrels;

 figure(1092)
        clf
    t4=tiledlayout('flow');
    


    for gap=1:gapl
       h4= nexttile;
        succplot=reshape(succ(:,1,gap,:),[sigl],[etolsl]);
pcolor(etols,sigs,succplot);
set(gca, 'XScale', 'log');
  set(gca, 'YScale', 'log');
     %shading interp
title("$\Delta t=$"+timestepbetweendata(gap),'Interpreter','latex','FontSize',12);
 xticks([1e-2,1e1])
  yticks([1e-3,1e0])
caxis manual
caxis([0 1]);


hold on

    end
xlabel(t4,'Error threshold, $\delta_r$','Interpreter','latex','FontSize',14)
ylabel(t4,'Noise level, $\sigma$','Interpreter','latex','FontSize',14)
 cbh = colorbar(h4(end)); 
% % To position the colorbar as a global colorbar representing
% % all tiles, 
 cbh.Layout.Tile = 'east'; 
%%
su=succ<0.50;
rstest=rsq;
rstest(su)=0;
figure(107)
clf
 t4=tiledlayout('flow')


    for gap=1:gapl
        
        ebysig=reshape(rstest(:,gap,:),[sigl],[length(etols)]);

figure(107)
nexttile
 [M,I]=max(ebysig,[],"all");
 [X,Y]=ind2sub(size(ebysig),I);
plot(tranend(:,X),tenderror(X,:,gap,Y),'o','Markersize',4)
hold on
plot(t,t,'k')
 title("$\Delta t=$"+timestepbetweendata(gap)+",$R^2=$"+num2str(M, 3),'Interpreter','latex','FontSize',12);


    end
    
xlabel(t4,'Computed transient end time','Interpreter','latex','FontSize',14)
ylabel(t4,'Predicted transient end','Interpreter','latex','FontSize',14)






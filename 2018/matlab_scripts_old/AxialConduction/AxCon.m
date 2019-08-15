%%% 1D heat conduction in pipe + wall

% cells=5500;
mal=1;

cells=550*mal;
xxx=100*mal;
xx1=120*mal;
xx2=320*mal;
xx3=12/10*mal;
xx4=2/3;
%%xx4=0.;
ccc=2*sqrt(2);

CI1=0.11;

CI2=2/3;
%CI2=0.65;
%CI2=0.635;

CI3 = 12.7;
%%CI3 = 14;


L=5.5;
dz=L/cells;
z=linspace(0,L,cells);

Tend=0.001;               %sec
Tend=200;               %sec
dt=4E-2;            %sec
% dt=0.001;            %sec
dt=0.02;            %sec
di=12E-3;             %m
da=15E-3;             %m
V=(da^2-di^2)*pi/4*2; %m^3
CN=0.5;

TmS=TPT100(5);
%um=0.5230; %m/s
mdot; %kg/m3
% U=6.0969; %V
% I=355.8913; %A
% PP=U*I;
PP=Qin_out; %%% Qin_out
%PP=Q_Pt100; %%% Qin_out
% qv=0.12;
qv=0;
%kHL = dqhldT / (da*pi*5.45); %%% heat loss
%r_iso=0.29/(2*pi);
r_iso=0.05;
kHL0 = 1 ./ (ra/0.055 * log(r_iso/ra) );
dTmOH = kHL0*2*pi*ra*5.45*(124.5-T_inf)/(mdot*cpm);

PP=PP+kHL0*2*pi*ra*5.45*(Tw_SA-T_iso)-(lambda_Pet*(5.45)/di + zeta_mischer) * rho_m * Wmean^2/2 * Wmean * Ageo_quer;


% %DNS
% TmS=196.85;
% mdot=0.053811;
% PP=1.77924E4*2*pi*di*2;
% kHL0=0.;

T_inf %%% heat loss

qQ=ones(1,cells);
for i=1:cells
%    if(i>1200 && i<3200)
     if(i>xx1 && i<xx2)
        qQ(i)=PP / (V); %W/m^3
    else
        qQ(i)=0; %W/m^3
    end
end
        

lambdaW=15;  %W/mK
rhoW=7.93E3; %kg/m3
cpW=500;     %J/kgK

aW= lambdaW/(rhoW*cpW);

dAdVW = di*4/(da^2-di^2);
dAdVWA = da*4/(da^2-di^2);

dAdVF = 4/(di);



TmS=TmS+273.15;

Tw=ones(1,cells)*TmS;
Tm=ones(1,cells)*TmS;

lambdaF=calc_lambda(Tm-273.15);  %W/mK
rhoF=calc_rho(Tm-273.15);        %kg/m3
cpF=calc_cp(Tm-273.15)*1000;          %J/kgK
muF=calc_mu(Tm-273.15);           %Pas
PrF = muF.*cpF./lambdaF;    %-

lambdaFW=calc_lambda(Tw-273.15);  %W/mK
rhoFW=calc_rho(Tw-273.15);        %kg/m3
cpFW=calc_cp(Tw-273.15)*1000;          %J/kgK
muFW=calc_mu(Tw-273.15);           %Pas
PrFW = muFW.*cpFW./lambdaFW;    %-

aF= lambdaF./(rhoF.*cpF);

um=ones(1,cells);

NuGnie=ones(1,cells);
alpha=ones(1,cells);
qw=ones(1,cells);
qv=ones(1,cells);
khl=ones(1,cells);

qAw=ones(1,cells);
qAf=ones(1,cells);


qAw0=dt*aW/dz*1/dz*CN;
  aaa=1;  

        %%% ....
zeta_10000 = (0.79*log(10000)-1.64)^(-2);

t=0;
while t < Tend
    Tm0=Tm;
    Tw0=Tw;
    
    qAw(1) = 0;
    qAf(1) = 0;
    Tm(1) = 0;
    Tw(1) = 0;
    for i=2:(cells-1)
        %% RHS %% EXPLICIT PART
        um(i) = mdot / ( rhoF(i)* (di^2*pi/4) );
        %% additional flux term
        Re_m=um(i)*di*rhoF(i)/muF(i);
        %PrF(i)=Pr_m;
        
        zeta=(0.79*log(Re_m)-1.64)^(-2);
        %NuGnie(i)=((zeta/8)*(Re_m-1000)*PrF(i))/(1+12.7*sqrt(zeta/8)*(PrF(i)^(2/3)-1))*(PrF(i)/PrFW(i))^0.11;
        
        if(i>xx1 && i<xx2)
            NuGnie(i)=((zeta/8)*(Re_m-1000)*PrF(i))/(1+CI3*sqrt(zeta/8)*(PrF(i)^(CI2)-1))*(PrF(i)/PrFW(i))^CI1*(1+(xx3/(i-xxx))^(xx4));
            khl(i)=kHL0;
        else
            NuGnie(i)=((zeta/8)*(Re_m-1000)*PrF(i))/(1+CI3*sqrt(zeta/8)*(PrF(i)^(CI2)-1))*(PrF(i)/PrFW(i))^CI1;
            khl(i)=kHL0;
        end
        
%         if( (i>4 && i<6) || (i>315 && i<325) || (i>115 && i<125) )
%             khl(i)=kHL0*10;
%         end
        
        
        %%% ....
%         Nu_m_L_02300_Gnie = (83.326+(1.953*(2300*PrF(i)*di/L)^(1/3)-0.6)^3+(0.924*PrF(i)^(1/3)*(2300*di/L)^(1/2))^3)^(1/3);
%         Nu_m_L_10000_Gnie = ((zeta_10000/8)*10000*PrF(i))/(1+12.7*sqrt(zeta_10000/8)*(PrF(i)^(2/3)-1));%*(Pr_m/Pr_w)^0.11;
%         gam = (Re_m-2300)/(10000-2300);
%          NuGnie(i)= (1-gam) * Nu_m_L_02300_Gnie + gam * Nu_m_L_10000_Gnie;

        
        alpha(i)= NuGnie(i)*lambdaF(i)/di;
        
        qw(i) = alpha(i)*(Tw(i)-Tm(i));              %% Heatflux form the pipe wall to the fluid
        %qv(i) = kHL * (Tw(i) - (T_inf+273.15));    %% Heatloss over the isolation (T_iso from measurement)
        qv(i) = khl(i) * (Tw(i) - (T_inf+273.15));    %% Heatloss over the isolation (T_iso from measurement)
        
        %% 1. wall
        Tw(i) = Tw0(i)                                                      + ...
            dt * (                                                          + ...
                    (1-CN) * aW*1/dz*( (Tw0(i+1)-Tw0(i))/dz - (Tw0(i)-Tw0(i-1))/dz ) + ...
                    qQ(i)/(rhoW*cpW)                                        - ...
                    qw(i) * dAdVW/(rhoW*cpW)                                - ...
                    qv(i) * dAdVWA/(rhoW*cpW)... %0. ...%
                 );
        %% 2. fluid
        Tm(i) = Tm0(i)                                                        + ...
            dt * (                                                            - ...
                    (1-CN) * um(i)*1/dz*( (Tm0(i+1)+Tm0(i))/2  - (Tm0(i)+Tm0(i-1))/2  )   + ...
                    (1-CN) * aF(i)*1/dz*( (Tm0(i+1)-Tm0(i))/dz - (Tm0(i)-Tm0(i-1))/dz) + ...
                    qw(i) * dAdVF/(rhoF(i)*cpF(i)) ...
                 );
        %% IMPLICIT PART
        %% 1. wall
        if(i==2)
            Tw(i)=Tw(i)+qAw0*Tw0(1);
        end
        if(i==(cells-1))
            Tw(i)=Tw(i)+qAw0*Tw0(end);
        end
        
        qAw(i) = -qAw0 / ( (1+2*qAw0) + qAw(i-1)*qAw0 );
        Tw(i) = ( Tw(i) + Tw(i-1)*qAw0 ) / ( (1+2*qAw0) + qAw(i-1)*qAw0 );
        
        %% 2. fluid
        qAf0=dt*(aF(i)/dz)*1/dz*CN;
        qAf1=dt*(um(i)/2)*1/dz*CN;
        
        if(i==2)
            Tm(i)=Tm(i)+(qAf0+qAf1)*Tm0(1);
        end
        if(i==(cells-1))
            Tm(i)=Tm(i)+(qAf0-qAf1)*Tm0(end);
        end
    
        qAf(i) = -(qAf0-qAf1) / ( (1+2*qAf0) + qAf(i-1)*(qAf0+qAf1) );
        Tm(i) = ( Tm(i) + Tm(i-1)*(qAf0+qAf1) ) / ( (1+2*qAf0) + qAf(i-1)*(qAf0+qAf1) );
        
    end
    Tw(1) = Tw0(1);
    Tm(1) = Tm0(1);
    %% solving system
    for i=(cells-2):-1:2
         Tw(i) = Tw(i) - qAw(i)*Tw(i+1);
         Tm(i) = Tm(i) - qAf(i)*Tm(i+1);
    end
    
    Tw(end)=Tw(end-1);
    Tm(end)=Tm(end-1);
    
    lambdaF=calc_lambda(Tm-273.15);  %W/mK
    rhoF=calc_rho(Tm-273.15);        %kg/m3
    cpF=calc_cp(Tm-273.15)*1000;          %J/kgK
    muF=calc_mu(Tm-273.15);           %Pas
    PrF = muF.*cpF./lambdaF;    %-

    lambdaFW=calc_lambda(Tw-273.15);  %W/mK
    rhoFW=calc_rho(Tw-273.15);        %kg/m3
    cpFW=calc_cp(Tw-273.15)*1000;          %J/kgK
    muFW=calc_mu(Tw-273.15);           %Pas
    PrFW = muFW.*cpFW./lambdaFW;    %-

    aF= lambdaF./(rhoF.*cpF);
    
    cflU=um(i)*dt/dz;
    cflNu=muF(end)/rhoF(end)*dt/dz^2;

    if(cflU>ccc)
        dt=dt/1.1;
    elseif(cflU<ccc)
        dt=dt*1.1;
    end
    
    fprintf('t=%d cfl(U)=%d  cfl(nu)=%d \n',t,cflU,cflNu);
    t=t+dt;    
end

%%

%%
    %% %%% Bound it!
    i=1;
        um(i) = mdot / ( rhoF(i)* (di^2*pi/4) );
        Re_m=um(i)*di*rhoF(i)/muF(i);
        
        zeta=(0.79*log(Re_m)-1.64)^(-2);
        %NuGnie(i)=((zeta/8)*(Re_m-1000)*PrF(i))/(1+12.7*sqrt(zeta/8)*(PrF(i)^(2/3)-1))*(PrF(i)/PrFW(i))^0.11;
        
        if(i>xx1 && i<xx2)
            NuGnie(i)=((zeta/8)*(Re_m-1000)*PrF(i))/(1+12.7*sqrt(zeta/8)*(PrF(i)^(2/3)-1))*(PrF(i)/PrFW(i))^(CI2)*(1+(xx3/(i-xxx))^(xx4));
        else
            NuGnie(i)=((zeta/8)*(Re_m-1000)*PrF(i))/(1+12.7*sqrt(zeta/8)*(PrF(i)^(2/3)-1))*(PrF(i)/PrFW(i))^(CI2);
        end
        
%         Nu_m_L_02300_Gnie = (83.326+(1.953*(2300*PrF(i)*di/L)^(1/3)-0.6)^3+(0.924*PrF(i)^(1/3)*(2300*di/L)^(1/2))^3)^(1/3);
%         zeta_10000 = (0.79*log(10000)-1.64)^(-2);
%         Nu_m_L_10000_Gnie = ((zeta_10000/8)*10000*PrF(i))/(1+12.7*sqrt(zeta_10000/8)*(PrF(i)^(2/3)-1));%*(Pr_m/Pr_w)^0.11;
%         gamma = (Re_m-2300)/(10000-2300);
%         Nu_m_Gnie = (1-gamma) * Nu_m_L_02300_Gnie + gamma * Nu_m_L_10000_Gnie;
%         
%         NuGnie(i)=Nu_m_Gnie;
        
        alpha(i)= NuGnie(i)*lambdaF(i)/di;
        
        qw(i) = alpha(i)*(Tw(i)-Tm(i));
    i=cells;
        um(i) = mdot / ( rhoF(i)* (di^2*pi/4) );
        Re_m=um(i)*di*rhoF(i)/muF(i);
        zeta=(0.79*log(Re_m)-1.64)^(-2);
        %NuGnie(i)=((zeta/8)*(Re_m-1000)*PrF(i))/(1+12.7*sqrt(zeta/8)*(PrF(i)^(2/3)-1))*(PrF(i)/PrFW(i))^0.11;
        
        if(i>xx1 && i<xx2)
            NuGnie(i)=((zeta/8)*(Re_m-1000)*PrF(i))/(1+12.7*sqrt(zeta/8)*(PrF(i)^(2/3)-1))*(PrF(i)/PrFW(i))^(CI2)*(1+(xx3/(i-xxx))^(xx4));
        else
            NuGnie(i)=((zeta/8)*(Re_m-1000)*PrF(i))/(1+12.7*sqrt(zeta/8)*(PrF(i)^(2/3)-1))*(PrF(i)/PrFW(i))^(CI2);
        end
        
%         Nu_m_L_02300_Gnie = (83.326+(1.953*(2300*PrF(i)*di/L)^(1/3)-0.6)^3+(0.924*PrF(i)^(1/3)*(2300*di/L)^(1/2))^3)^(1/3);
%         zeta_10000 = (0.79*log(10000)-1.64)^(-2);
%         Nu_m_L_10000_Gnie = ((zeta_10000/8)*10000*PrF(i))/(1+12.7*sqrt(zeta_10000/8)*(PrF(i)^(2/3)-1));%*(Pr_m/Pr_w)^0.11;
%         gamma = (Re_m-2300)/(10000-2300);
%         Nu_m_Gnie = (1-gamma) * Nu_m_L_02300_Gnie + gamma * Nu_m_L_10000_Gnie;
%         
%         NuGnie(i)=Nu_m_Gnie;
        
        alpha(i)= NuGnie(i)*lambdaF(i)/di;
        
        qw(i) = alpha(i)*(Tw(i)-Tm(i));


%%  Temperatures
fig = figure;

zF = 1.2:0.1:3.2;
kF = (Tm(end)-Tm(1)) / (2);
%dF = Tw(2970) - kF*2.97;
dF = Tw(297*mal) - kF*2.97;

plot(z,Tm-273.15,'DisplayName','Tm');
hold all;
plot(z,Tw-273.15,'DisplayName','Tw');
plot(zF,kF*zF+dF-273.15,'-.k','DisplayName','fully developed');
plot(z(297),Tm(297)-273.15,'xk','LineWidth',5,'markerSize',20);
plot(z(297),Tw(297)-273.15,'xk','LineWidth',5,'markerSize',20);
plot(z(247),Tw(247)-273.15,'xk','LineWidth',5,'markerSize',20);

title('T / K')
set(gca,'FontSize',15);
set(gca,'YColor',[0 0 0],'YGrid','on');
set(gca,'XColor',[0 0 0],'XGrid','on');
set(gca,'YMinorGrid','on','YGrid','on',...
    'XMinorGrid','on',...
    'XGrid','on',...
    'YMinorTick','on',...
    'YColor',[0 0 0],...
    'XMinorTick','on',...
    'XColor',[0 0 0]);

ylabel('T / K')
xlabel('z / m')

%hleg=legend('show','Location','NorthEast');


%% dTemperatures
figure;

plot(z,Tw-Tm,'DisplayName','dT');

title('Temperature difference Tw-Tm')
set(gca,'FontSize',15);
set(gca,'YColor',[0 0 0],'YGrid','on');
set(gca,'XColor',[0 0 0],'XGrid','on');
set(gca,'YMinorGrid','on','YGrid','on',...
    'XMinorGrid','on',...
    'XGrid','on',...
    'YMinorTick','on',...
    'YColor',[0 0 0],...
    'XMinorTick','on',...
    'XColor',[0 0 0]);

ylabel('dT / K')
xlabel('z / m')
hleg=legend('show','Location','NorthEast');

% Nusselt number
figure;

plot(z,NuGnie,'DisplayName','Gnielinski');

title('Nusselt number')
set(gca,'FontSize',15);
set(gca,'YColor',[0 0 0],'YGrid','on');
set(gca,'XColor',[0 0 0],'XGrid','on');
set(gca,'YMinorGrid','on','YGrid','on',...
    'XMinorGrid','on',...
    'XGrid','on',...
    'YMinorTick','on',...
    'YColor',[0 0 0],...
    'XMinorTick','on',...
    'XColor',[0 0 0]);

ylabel('Nu / -')
xlabel('z / m')
hleg=legend('show','Location','SouthEast');

% alpha
figure;

plot(z,alpha,'DisplayName','alpha');

title('Heat transfer coefficient')
set(gca,'FontSize',15);
set(gca,'YColor',[0 0 0],'YGrid','on');
set(gca,'XColor',[0 0 0],'XGrid','on');
set(gca,'YMinorGrid','on','YGrid','on',...
    'XMinorGrid','on',...
    'XGrid','on',...
    'YMinorTick','on',...
    'YColor',[0 0 0],...
    'XMinorTick','on',...
    'XColor',[0 0 0]);

ylabel('alpha / W/m2K')
xlabel('z / m')
hleg=legend('show','Location','SouthEast');

% Heat Flux
figure;

plot(z,qw,'DisplayName','qw');

title('Wall heat flux')
set(gca,'FontSize',15);
set(gca,'YColor',[0 0 0],'YGrid','on');
set(gca,'XColor',[0 0 0],'XGrid','on');
set(gca,'YMinorGrid','on','YGrid','on',...
    'XMinorGrid','on',...
    'XGrid','on',...
    'YMinorTick','on',...
    'YColor',[0 0 0],...
    'XMinorTick','on',...
    'XColor',[0 0 0]);

ylabel('qw / W/m2')
xlabel('z / m')
hleg=legend('show','Location','NorthEast');

% Velocity
figure;

plot(z,um,'DisplayName','um');

title('Velocity')
set(gca,'FontSize',15);
set(gca,'YColor',[0 0 0],'YGrid','on');
set(gca,'XColor',[0 0 0],'XGrid','on');
set(gca,'YMinorGrid','on','YGrid','on',...
    'XMinorGrid','on',...
    'XGrid','on',...
    'YMinorTick','on',...
    'YColor',[0 0 0],...
    'XMinorTick','on',...
    'XColor',[0 0 0]);

ylabel('um / m/s')
xlabel('z / m')
hleg=legend('show','Location','NorthEast');

%%
% check mdot * cp * (Tm_in - Tm_out)
    Ac = 0.818;
    Bc = 3.664E-3;
    
T_in  = Tm(1);
T_out = Tm(end);
cpm=(Ac*(T_out-T_in)+Bc/2*(T_out^2-T_in^2))/(T_out-T_in)*1000;
%cpm=calc_cp(T_out)-273.15)*1000;

%mdot=m_dot_C1/3600; %in kg/s
mu_m = muF(end);
mdot=Re_m*pi*mu_m*di/4;
Q_into=mdot*cpm*(T_out-T_in)
dP=PP-Q_into

%%

if(0)
    %%
figure;

plot(z,Tm,'DisplayName','Tm');

title('T / K')
set(gca,'FontSize',15);
set(gca,'YColor',[0 0 0],'YGrid','on');
set(gca,'XColor',[0 0 0],'XGrid','on');
xlabel('z / m')
hleg=legend('show','Location','NorthEast');
end

%%
if(0)
    %%
figure;

plot(z,Tw,'DisplayName','Tw');

title('T / K')
set(gca,'FontSize',15);
set(gca,'YColor',[0 0 0],'YGrid','on');
set(gca,'XColor',[0 0 0],'XGrid','on');
xlabel('z / m')
hleg=legend('show','Location','NorthEast');
end
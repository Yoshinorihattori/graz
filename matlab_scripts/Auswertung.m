clc; clear all;

win=0;


%%%%%%%%% read measurement data
%HERE = 'Messreihen/'
HERE = './re270pr50/Daten_19.11.2018_15.42.44'

[stmfile, stmpath] =uigetfile(strcat(HERE,'*.txt'),'Prüfstandsdatei auswählen');
fid=fopen(fullfile(stmpath, stmfile),'rt');

date = stmfile(7:16);
time = strcat(stmfile(18:19) , ':' , stmfile(21:22) , ':' , stmfile(24:25));

indata=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',1);
fclose(fid);

% 1:  Tsa
% 2:  Tsi
% 3:  Ti
% 4:  m_dot C1
% 5:  m_dot C2
% 6:  Re C1
% 7:  Re C2
% 8:  v C1
% 9:  v Cerr2
% 10: U_MS
% 11: I_MS
% 12: P_MS
% 13: q
% 14: dp: m_dot
% 15: dt: Re
% 16: dp
% 17: nu Fluid
% 18: rho Fluid
% 19: lambda Fluid
% 20: cp Fluid
% 21: Pr Fluid
% 22: lambda Rohr
% 23: Nu
% 22: dp: T1
% 23: dp: T2
% 24: dp: T3

di=12E-3;
da=15E-3;
disoa = 90E-3;
L=2;
V=(da^2-di^2)*pi/4*L;
Ageo_MS=di*pi*L;
Ugeo_MS=di*pi;
Ageo_quer=di^2*pi/4;

Tsa(1,:)          = indata{1}(1:24);
Tsi(1,:)          = indata{2}(1:24);

Ti(1,:)           = indata{3}(1:8);
m_dot_C1(1,:)     = indata{4}(1);
m_dot_C2(1,:)     = indata{5}(1);
Re_C1(1,:)        = indata{6}(1);
Re_C2(1,:)        = indata{7}(1);
v_C1(1,:)         = indata{8}(1);
v_C2(1,:)         = indata{9}(1);
U_MS(1,:)         = indata{10}(1);
I_MS(1,:)         = indata{11}(1);
P_MS(1,:)         = indata{12}(1);
q(1,:)            = indata{13}(1);
P1(1,:)           = indata{14}(1);
P2(1,:)           = indata{15}(1);
dp(1,:)           = indata{16}(1);
nu_Fluid(1,:)     = indata{17}(1);
rho_Fluid(1,:)    = indata{18}(1);
lambda_Fluid(1,:) = indata{19}(1);
cp_Fluid(1,:)     = indata{20}(1);
Pr_Fluid(1,:)     = indata{21}(1);
lambda_Rohr(1,:)  = indata{22}(1);
Nu_Fluid(1,:)     = indata{23}(1);

I1(1,:)           = indata{24}(1);
U1(1,:)           = indata{25}(1);

dp_T1(1,:)        = indata{26}(1);
dp_T2(1,:)        = indata{27}(1);
dp_T3(1,:)        = indata{28}(1);
%Mittel (U*I) mean P / W
Pm(1,:)           = indata{29}(1);
TPT100_a(1,:)     = indata{30}(1:5);
TPT100_i(1,:)     = indata{31}(1:5);
TPT100_m(1,:)     = indata{32}(1:5); 

NuPt100(1,:)      = indata{33}(1);
NuQm(1,:)         = indata{34}(1);
%% zeta computation Einlesen
ZRem100(1,:)      = indata{35}(1);
Zzeta100(1,:)     = indata{36}(1);
Zqm100(1,:)       = indata{37}(1);

TPT100_T5(1,:)    = indata{38}(1);

%% T_aussen -> T_innen (check)
lambda_K = 0.16;
ri       = di/2;
ra       = da/2;
DeltaK   = 0.065E-3;
%lambda_N = @(T) 20;%20;%0.0135*T + 15.1622;%14;



%%% conduction equation for inner wall temperatuer
lambda_N = 20;

T1 = @(r,qv,qzu,TN,Ta)    qzu/(2*lambda_N)*ra^2.*(0.5-0.5.*(r/ra).^2+log(r/ra)) + Ta + ...
                      - qv/lambda_N*(ra+DeltaK).*(lambda_N/lambda_K*log(ra/(ra+DeltaK))+log(r/ra));
T2 = @(r,qv,Ta)      - qv/lambda_K*(ra+DeltaK).*log(r/(ra+DeltaK))+Ta;


if(0)
    %%
figure;
hold all;
r1=meshgrid(ri:0.005E-3:ra,1);
r2=meshgrid(ra:0.005E-3:(ra+DeltaK),1);
% plot(r1,T1(r1,0,3400/(Ageo_quer*L),100,100)-100,'-k','DisplayName','T1K qv=0');
% plot(r2,T2(r2,0,100)-100,'.k','DisplayName','T2K qv=0');
% 
% plot(r1,T1(r1,100,3400/(Ageo_quer*L),100,100)-100,'-r','DisplayName','T1K qv=100');
% plot(r2,T2(r2,100,100)-100,'.r','DisplayName','T2K qv=100');
% 
% plot(r1,T1(r1,1000,3400/(Ageo_quer*L),100,100)-100,'-b','DisplayName','T1K qv=1000');
% plot(r2,T2(r2,1000,100)-100,'.b','DisplayName','T2K qv=1000');


plot(r1,T1(r1,qvA,Qin_out/V,TPT100_a(3),TPT100_a(3)),'-b','DisplayName','T1K qv');
hold all;
plot(r2,T2(r2,qvA,TPT100_a(3)),'.b','DisplayName','T2K qv');


set(gca,'FontSize',15);
set(gca,'YColor',[0 0 0],'YGrid','on');
set(gca,'XColor',[0 0 0],'XGrid','on');
h=xlabel('r / m');
h1=ylabel('T - T_w / gradC');
end



%% Definition der Temperaturen und MesspositionenCoeff_lf
Tsi_up = [Tsi(1,1) Tsi(1,3) Tsi(1,5) Tsi(1,7) Tsi(1,10) Tsi(1,11) Tsi(1,13), ...
    Tsi(1,15) Tsi(1,17) Tsi(1,19) Tsi(1,20) Tsi(1,21) Tsi(1,23)];
Tsa_up = [Tsa(1,1) Tsa(1,3) Tsa(1,5) Tsa(1,7) Tsa(1,10) Tsa(1,11) Tsa(1,13), ...
    Tsa(1,15) Tsi(1,17) Tsi(1,19) Tsi(1,20) Tsi(1,21) Tsi(1,23)];

x_pos_Ts_up = [1.230 1.520 2.020 2.420 2.720 2.870 3.020 3.170 3.230 4.680 ...
    4.720 5.480 5.520];
x_pos_Ts_up_a = [1.230 1.520 2.020 2.420 2.720 2.870 3.020 3.170 3.230 4.680 ...
    4.720 5.480 5.520];

Tsi_down = [Tsi(1,2) Tsi(1,4) Tsi(1,6) Tsi(1,8) Tsi(1,12) Tsi(1,14) ...
    Tsi(1,16) Tsi(1,18) Tsi(1,22) Tsi(1,24)];
Tsa_down = [Tsa(1,2) Tsa(1,4) Tsa(1,6) Tsa(1,8) Tsa(1,12) Tsa(1,14) ...
    Tsa(1,16) ];
x_pos_Ts_down = [1.230 1.520 2.020 2.420 2.870 3.020 3.170 3.230 5.480 5.520];
x_pos_Ts_down_a = [1.230 1.520 2.020 2.420 2.870 3.020 3.170];

x_pos_t=[2.720 4.700];

T_Ms_fluid = [ Ti(3) Ti(1)];

x_pos_p=[3.3 4.3 4.4];

pos = 1:24;



x_pos_TPt100 = [2.920 0.05 2.720 5.500];
x_pos_TPt100 = [2.420 0.05 2.920 5.500];
x_pos_TPt100 = [1.970 2.470 2.970 5.500 0.05];
xTin  = 0.05;
xTout = 5.50;
%x_pos_TPt100 = [2.720 3.170 3.170 4.700];

x_MS_in = 1.2;
x_MS_out = 3.2;


%%



%%%%%%%%%% TPT100_a = TPT100_mean


TPT100=zeros(5,1);

TPT100(1:5)=TPT100_i';

% TPT100(2)=TPT100_a(2);
% TPT100(4)=TPT100_a(4);
% TPT100(5)=TPT100_T1;

Tm_aus=TPT100(4);
Tm=Tm_aus; %%%%%%%%%%%%%%% neu

%% iteration to update the temperatuer dependent material properties --> new wall temperature
for i=1:2

    Tw = TPT100(3);

    T_iso = Ti(2);
    T_inf = Ti(7);


%% Berechnung der Fluid Stoffwerte

    Pr_w     = calc_Pr(Tw);
    rho_w    = calc_rho(Tw);
    nu_w     = calc_nu(Tw);
    mu_w     = nu_w * rho_w;
    cp_w     = calc_cp(Tw);
    lambda_w = calc_lambda(Tw);

    Pr_m     = calc_Pr(Tm);
    rho_m    = calc_rho(Tm);
    nu_m     = calc_nu(Tm);
    mu_m     = nu_m * rho_m;
    cp_m     = calc_cp(Tm);
    lambda_m = calc_lambda(Tm);
    
    
    Re_C1;

%% Korrektur mittel Temperatur
    dqhldT = 1.;
    zeta_mischer = 32.46;

    Thl = (Tw+Tm)*0.5;
    qhl= dqhldT* (Thl-T_iso);

    qhlm2 = qhl / (da*pi*5.45);

    Re_m = Re_C1;
    lambda_Pet = (1.8*log10(Re_m)-1.5)^(-2.);

%%% cp mean
    T_in=TPT100(5)+273.15;
    T_out=TPT100(4)+273.15;
        Ac = 0.818;
        Bc = 3.664E-3;
    cpm=(Ac*(T_out-T_in)+Bc/2*(T_out^2-T_in^2))/(T_out-T_in)*1000;

    Wmean=Re_m*nu_m/di;

    mdot = m_dot_C1/3600;


%%% Q= DT_pt100 / Dx
    
    fit=polyfit(x_pos_TPt100(1:3),TPT100(1:3)',1);
    gradT_PT100=fit(1);
    
    Q_Pt100 = gradT_PT100*mdot*cpm*L;

    qUI=P_MS/Ageo_MS;

    dqw2 = (lambda_Pet*(5.45-2.92)/di + zeta_mischer) * rho_m * Wmean^2/2 * Wmean * Ageo_quer - qhl*(5.45-2.53)/5.45;

    %Qin_out=max(mdot*cpm*(T_out-T_in)-(lambda_Pet*(5.45-2.92)/di + zeta_mischer) * rho_m * Wmean^2/2 * Wmean * Ageo_quer + qhl,0);
    Qin_out=mdot*cpm*(T_out-T_in);
    dQ=Qin_out-P_MS;

%%% Q=m_dot * cp_m * DT;
    gradT = (Qin_out/L) / ( (mdot) * cpm );

    dTm = dqw2/(mdot*cpm);


    SA_Tm=gradT*(x_MS_out-x_pos_TPt100(3));

%%%%%%

    Tm_aus = TPT100(4) - dTm;
    Tm = Tm_aus - SA_Tm;


    qw = Qin_out/ Ageo_MS;
    if(Qin_out>0)
        qvA=qhlm2;
    else
        qvA=0;
    end


    Pr_w     = calc_Pr(Tw);
    rho_w    = calc_rho(Tw);
    nu_w     = calc_nu(Tw);
    mu_w     = nu_w * rho_w;
    cp_w     = calc_cp(Tw);
    lambda_w = calc_lambda(Tw);

    Pr_m     = calc_Pr(Tm);
    rho_m    = calc_rho(Tm);
    nu_m     = calc_nu(Tm);
    mu_m     = nu_m * rho_m;
    cp_m     = calc_cp(Tm);
    lambda_m = calc_lambda(Tm);

%P_spez_m2 = P_MS/A;
    P=Pm;
    P=P_MS;

    P_spez_m2 = P/Ageo_MS;
    P_spez_m3 = P/V;

%% Berechnung von Re
    Re_m = Re_C1;
    Re_w = Re_C2;

    Re_m = m_dot_C1/3600*4/(nu_m*rho_m*di*pi);
    Re_w = m_dot_C1/3600*4/(nu_w*rho_m*di*pi);           %%% check why there is a differnece!!!!



%% Korrektur Tw @ Ende Messtrecke
%%% Korrektur
    qVol = Qin_out/V;

    Tw=T1(ri,qvA,qVol, (TPT100_a(3)), (TPT100_a(3)));
    Tw_SA=T1(ri,qvA,qVol, (TPT100_a(1)), (TPT100_a(1)));
    Tw_SA2=T1(ri,qvA,qVol, (TPT100_a(2)), (TPT100_a(2)));

    TPT100(3) = Tw;
    TPT100(1) = Tw_SA;
    TPT100(2) = Tw_SA2;

    Tsi_up(1:8)=T1(ri,qvA,Qin_out/V, Tsa_up(1:8) ,Tsa_up(1:8));
    Tsi_down(1:7)=T1(ri,qvA,Qin_out/V, Tsa_down(1:7) ,Tsa_down(1:7));

end

%%% Q= gradtent of Wall Temp sensors
xFit=x_pos_Ts_down(2:end-4);
TFit=Tsi_down(2:end-4);
TFit(1)=( TFit(1)+Tsi_up(2)  )/2;
TFit(2)=( TFit(2)+Tsi_up(3)  )/2;
TFit(3)=( TFit(3)+Tsi_up(4)  )/2;
TFit(4)=( TFit(4)+Tsi_up(6)  )/2;
TFit(5)=( TFit(5)+Tsi_up(7)  )/2;
fit=polyfit(xFit,TFit,1);
gradTfit=fit(1);

off = Tw - (gradTfit*x_pos_TPt100(3) + fit(2));

TFit     = TFit     + off;
Tsi_up   = Tsi_up   + off;
Tsi_down = Tsi_down + off;
Tsa_up   = Tsa_up   + off;
Tsa_down = Tsa_down + off;
    
%% Update Fluidproperties
%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Pr_m     = calc_Pr(Tm);
rho_m    = calc_rho(Tm);
nu_m     = calc_nu(Tm);
mu_m     = nu_m * rho_m;
cp_m     = calc_cp(Tm);
lambda_m = calc_lambda(Tm);
Re_m = m_dot_C1/3600*4/(nu_m*rho_m*di*pi);

Pr_w     = calc_Pr(Tw);
rho_w    = calc_rho(Tw);
nu_w     = calc_nu(Tw);
mu_w     = nu_w * rho_w;
cp_w     = calc_cp(Tw);
lambda_w = calc_lambda(Tw);

%% qscale


qScale     = ( gradTfit/gradT );
qScalePT100= ( gradT_PT100/gradT );




%% plot
if(1)
figure;
hold all;
plot(x_pos_TPt100,TPT100,'xr','DisplayName','TPt','MarkerSize',10,'LineWidth',4);
line( [x_pos_TPt100(5) x_MS_in],[TPT100(5) TPT100(5)]);
line( [x_MS_in x_MS_out],[TPT100(5) Tm_aus]);
line( [x_MS_out x_pos_TPt100(4)],[Tm_aus TPT100(4)]);

plot(x_pos_Ts_up,Tsi_up,'xg-','DisplayName','Tsi_o','LineWidth',2);
plot(x_pos_Ts_down,Tsi_down,'xm-','DisplayName','Tsi_u','LineWidth',2);

plot(x_pos_Ts_up_a,Tsa_up,'xg-.','DisplayName','Tsa_o');
plot(x_pos_Ts_down_a,Tsa_down,'xm-.','DisplayName','Tsa_u');

offset = Tw - x_pos_TPt100(3) * gradT;
plot(x_pos_Ts_up,gradT*x_pos_Ts_up+offset,'-.k','DisplayName','fully developed flow Qin_out');

offset = Tw -  x_pos_TPt100(3)* gradTfit;
plot(xFit,gradTfit*xFit+offset,'-*k','DisplayName','fully developed flow fit');

offset = Tw -  x_pos_TPt100(3) * gradT_PT100;
plot(x_pos_TPt100(1:3),gradT_PT100*x_pos_TPt100(1:3)+offset,'-*c','DisplayName','Pt100 fit');

plot(x_pos_TPt100(3),Tm,'kx','DisplayName','SA Tm','MarkerSize',15,'LineWidth',4);



set(gca,'FontSize',15);
set(gca,'YColor',[0 0 0],'YGrid','on');
set(gca,'XColor',[0 0 0],'XGrid','on');
h=xlabel('Messposition x / m');
h1=ylabel('Temperatur T / °C');

saveas(gcf,'temperatures','epsc');

end

%% nusselt
%Messung
Nu_ideal =  qUI* di / (lambda_m * ( Tw-Tm ) );
Nu_QmdotCpDt_m = qw * di / (lambda_m * ( Tw-Tm ) );
Nu_QmdotCpDt_w = qw * di / (lambda_w * ( Tw-Tm ) );
qw_mdotCpDt = qw;

Nu_fit   = gradTfit   *mdot*cpm/(di*pi) * di / ( lambda_m * ( Tw-(Tm_aus - gradTfit   *(x_MS_out-x_pos_TPt100(3)) ) ) );
qw_fit = gradTfit   *mdot*cpm/(di*pi);
Nu_Pt100 = gradT_PT100*mdot*cpm/(di*pi) * di / ( lambda_m * ( Tw-(Tm_aus - gradT_PT100*(x_MS_out-x_pos_TPt100(3)) ) ) );
qw_Pt100 = gradT_PT100*mdot*cpm/(di*pi);

NuMessung=Nu_QmdotCpDt_w;

% Gnielinski
zeta_Kon=(0.79*log(Re_m)-1.64)^(-2);
zeta_Pet=(1.8*log10(Re_m)-1.5)^(-2.);
%%% Colebrook (1939)
% syms zeta;
% rau=0.015E-3;
% rau=0.015E-3;
% Cole = 1/sqrt(zeta)==-2.0*log10((rau/di)/3.7+2.51/(Re_m*sqrt(zeta)));
% zeta_Cole = double(vpasolve(Cole,zeta));

Wmean=Re_m*nu_m/di;
zeta_m = dp / ( (1.0)/di * rho_m * Wmean^2/2);

zeta= zeta_Pet;
% zeta= zeta_Cole;

NuGnie=((zeta/8.)*Re_m*Pr_m)/(1.+12.7*sqrt(zeta/8.)*(Pr_m^(2./3.)-1.))*(Pr_m/Pr_w)^0.11;
NuGnie1000=((zeta/8)*(Re_m-1000)*Pr_m)/(1+12.7*sqrt(zeta/8)*(Pr_m^(2/3)-1))*(Pr_m/Pr_w)^0.11;

% Mean
Nu_m_L_02300_Gnie = (83.326+(1.953*(2300*Pr_m*di/L)^(1/3)-0.6)^3+(0.924*Pr_m^(1/3)*(2300*di/L)^(1/2))^3)^(1/3);
zeta_Kon_10000=(0.79*log(10000)-1.64)^(-2);
zeta_Pet_10000=(1.8*log10(10000)-1.5)^(-2);
zeta_10000 = zeta_Pet_10000;
Nu_m_L_10000_Gnie = ((zeta_10000/8)*10000*Pr_m)/(1+12.7*sqrt(zeta_10000/8)*(Pr_m^(2/3)-1))*(Pr_m/Pr_w)^0.11;
gamma = (Re_m-2300)/(10000-2300);
Nu_m_Gnie = (1-gamma) * Nu_m_L_02300_Gnie + gamma * Nu_m_L_10000_Gnie;


% Wall
Nu_w_L_02300_Gnie = (83.326+(1.953*(2300*Pr_w*di/L)^(1/3)-0.6)^3+(0.924*Pr_w^(1/3)*(2300*di/L)^(1/2))^3)^(1/3);
Nu_w_L_10000_Gnie = ((zeta_10000/8)*10000*Pr_w)/(1+12.7*sqrt(zeta_10000/8)*(Pr_w^(2/3)-1))*(Pr_w/Pr_w)^0.11;
gamma = (Re_w-2300)/(10000-2300);
Nu_w_Gnie = (1-gamma) * Nu_w_L_02300_Gnie + gamma * Nu_w_L_10000_Gnie;


% DB
NuDB=0.023*Re_m^(4/5)*Pr_m^0.4;


fprintf('\nmdot\n');
fprintf('Messung qw=(mDotCpDt) ... ReM: %10.4f\n',Re_m);
%fprintf('qw = %10.4f %10.4f %10.4f\n',qw_mdotCpDt,qw_fit,qw_Pt100);
fprintf('qw = %10.4f\n',qw_mdotCpDt);
%fprintf('Messung       : %10.4f ; %10.4f ; %10.4fn',Nu_QmdotCpDt_m,Nu_Pt100,Nu_fit);
fprintf('Messung       : %10.4f\n',NuMessung);
fprintf('NuGnie vt     : %10.4f\n',NuGnie);
fprintf('NuGnie -1000  : %10.4f\n',NuGnie1000);
fprintf('Nu     mean   : %10.4f\n',Nu_m_Gnie);
fprintf('Nu     DB     : %10.4f\n',NuDB);
fprintf('\n');


%% Messfehler Nu
DT=0.2;
Dm=0.002*mdot;
Dlambda=7.328E-5*DT;
Dcp=3.664E-3*DT;

DNu_Nu=sqrt((Dm/mdot)^2+(Dcp/cpm)^2+(Dlambda/lambda_w)^2 + DT^2*((1/(TPT100(5)-Tm))^2+(TPT100(5)-Tw)^2/((-Tw+Tm)*(TPT100(5)-Tm))^2+(1/(Tw-Tm))^2));
Nu_QmdotCpDt_w*DNu_Nu;
%% Neues Tm

Pr_m     = calc_Pr(Tm);
rho_m    = calc_rho(Tm);
nu_m     = calc_nu(Tm);
mu_m     = nu_m * rho_m;
cp_m     = calc_cp(Tm);
lambda_m = calc_lambda(Tm);
Re_m = m_dot_C1/3600*4/(nu_m*rho_m*di*pi);



Wmean=Re_m*nu_m/di;




zeta_Kon=(0.79*log(Re_m)-1.64)^(-2);
zeta_Pet=(1.8*log10(Re_m)-1.5)^(-2.);

cf_Pet = zeta_Pet/4;
% %%% Colebrook (1939)
% syms zeta;
% rau=0.015E-3;
% rau=0.015E-3;
% Cole = 1/sqrt(zeta)==-2.0*log10((rau/di)/3.7+2.51/(Re_m*sqrt(zeta)));
% zeta_Cole = double(vpasolve(Cole,zeta));

zeta= zeta_Pet;

% zeta=(0.79*log(Re_m)-1.64)^(-2);
dp_pet=zeta*1/di*rho_m*Wmean^2/2;
tau_w_pet=dp_pet*di/(4);
w_tau_pet=sqrt(tau_w_pet/rho_w);

zeta_messung=dp*di/rho_m/Wmean^2*2;
tau_w=dp*di/(1*4);

cf_M=tau_w/(rho_m*Wmean^2/2);

%% cf Messfehler
DDp=0.0035*dp;
Drho=0.616*DT;

Dcf_cf = sqrt((DDp/dp)^2+(Drho/rho_m)^2+(2*Dm/mdot)^2);
cf_M*Dcf_cf;
%%

w_tau=sqrt(tau_w/rho_m);

ReTau=w_tau*di/nu_w;
ReTau_pet=w_tau_pet*di/nu_w;

TTau=qw_mdotCpDt/(rho_w*cp_w*1000*w_tau);
TTau_pet=qw/(rho_w*cp_w*1000*w_tau_pet);
hTau=cp_w*1000*TTau;
hTau_pet=cp_w*1000*TTau_pet;



zeta_Kon=(0.79*log(Re_m)-1.64)^(-2);
dp_kon=zeta_Kon*1/di*rho_m*Wmean^2/2;

%% 

%%%% Output for further postProcessing
if(0)
  fprintf('%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f,%16.3f\n',Tm,TPT100(5),Tw,mdot,cpm,Re_m,Re_w,ReTau,Pr_m,Pr_w,P_MS,Qin_out,dQ,qw,qvA,dp,Nu_QmdotCpDt_w,Nu_Pt100,Nu_fit,NuGnie,NuGnie1000,Nu_m_Gnie,DNu_Nu,cf_M,cf_Pet,Dcf_cf);

  fprintf('%26.10f,%26.10f\n',Dcf_cf,DNu_Nu);
  fprintf('%26.10f,%26.10f,%26.10f\n',cf_M,cf_Pet,Dcf_cf);
endif

%%
if(win)
%% Schreiben LATEX
%save_dataInport.tex


fileID = fopen('../Auswertung/Report_neu/dataInport.tex','w');

fprintf(fileID,'\\newcommand{\\Date}{%s}\n',date);
fprintf(fileID,'\\newcommand{\\Time}{%s}\n',time);

fprintf(fileID,'\\newcommand{\\ReNumberTm}{%10.2f}\n',Re_m);
fprintf(fileID,'\\newcommand{\\ReNumberTw}{%10.2f}\n',Re_w);
fprintf(fileID,'\\newcommand{\\PrNumberTm}{%10.2f}\n',Pr_m);
fprintf(fileID,'\\newcommand{\\PrNumberTw}{%10.2f}\n',Pr_w);
fprintf(fileID,'\\newcommand{\\massflow}{%10.2f}\n',mdot);

fprintf(fileID,'\\newcommand{\\Frhom}{%10.2f}\n',rho_m);
fprintf(fileID,'\\newcommand{\\Flambdam}{%10.2f}\n',lambda_m);
fprintf(fileID,'\\newcommand{\\Fcpm}{%10.2f}\n',cp_m);
fprintf(fileID,'\\newcommand{\\Fnum}{%10.2E}\n',nu_m);

fprintf(fileID,'\\newcommand{\\Frhow}{%10.2f}\n',rho_w);
fprintf(fileID,'\\newcommand{\\Flambdaw}{%10.2f}\n',lambda_w);
fprintf(fileID,'\\newcommand{\\Fcpw}{%10.2f}\n',cp_w);
fprintf(fileID,'\\newcommand{\\Fnuw}{%10.2E}\n',nu_w);


fprintf(fileID,'\\newcommand{\\Wmean}{%10.4f}\n',Wmean);

fprintf(fileID,'\\newcommand{\\zetaM}{%10.4f}\n',zeta_messung);
fprintf(fileID,'\\newcommand{\\dpM}{%10.2f}\n',dp);
fprintf(fileID,'\\newcommand{\\cfM}{%14.8f}\n',zeta_messung/4);

fprintf(fileID,'\\newcommand{\\zetaPet}{%10.4f}\n',zeta);
fprintf(fileID,'\\newcommand{\\dpPet}{%10.2f}\n',dp_pet);
fprintf(fileID,'\\newcommand{\\cfPet}{%14.8f}\n',zeta/4);

fprintf(fileID,'\\newcommand{\\PMS}{%10.2f}\n',P_MS);
fprintf(fileID,'\\newcommand{\\PMSme}{%10.2f}\n',P_spez_m2);
fprintf(fileID,'\\newcommand{\\PMSma}{%10.2f}\n',P_spez_m3);
fprintf(fileID,'\\newcommand{\\UMS}{%10.2f}\n',U_MS);
fprintf(fileID,'\\newcommand{\\IMS}{%10.2f}\n',I_MS);
fprintf(fileID,'\\newcommand{\\qscale}{%10.2f}\n',(1-qScalePT100)*100);

fprintf(fileID,'\\newcommand{\\WallTemp}{%10.2f}\n',Tw);
fprintf(fileID,'\\newcommand{\\MeanTemp}{%10.2f}\n',Tm);

fprintf(fileID,'\\newcommand{\\WallTempK}{%10.2f}\n',Tw+273.15);
fprintf(fileID,'\\newcommand{\\MeanTempK}{%10.2f}\n',Tm+273.15);

fprintf(fileID,'\\newcommand{\\NusseltID}{%10.2f}\n',Nu_ideal);
fprintf(fileID,'\\newcommand{\\NusseltM}{%10.2f}\n',NuMessung);

fprintf(fileID,'\\newcommand{\\NusseltGnieSimpM}{%10.2f}\n',NuGnie);
fprintf(fileID,'\\newcommand{\\NusseltGnieSimpMT}{%10.2f}\n',NuGnie1000);
fprintf(fileID,'\\newcommand{\\NuGnieTm}{%10.2f}\n',Nu_m_Gnie);
fprintf(fileID,'\\newcommand{\\NuGnieTw}{%10.2f}\n',Nu_w_Gnie);
fprintf(fileID,'\\newcommand{\\NuDB}{%10.2f}\n',NuDB);

fprintf(fileID,'\\newcommand{\\Tmin}{%10.2f}\n',T_in-273.15);
fprintf(fileID,'\\newcommand{\\Tmout}{%10.2f}\n',T_out-273.15);


fprintf(fileID,'\\newcommand{\\DNSRetau}{%10.2f}\n',ReTau_pet);
fprintf(fileID,'\\newcommand{\\DNStauw}{%10.4f}\n',tau_w_pet);
fprintf(fileID,'\\newcommand{\\DNSwtau}{%10.4f}\n',w_tau_pet);
fprintf(fileID,'\\newcommand{\\DNSTfric}{%10.4f}\n',TTau_pet);
fprintf(fileID,'\\newcommand{\\DNSHfric}{%10.4f}\n',hTau_pet);



fprintf(fileID,'\\newcommand{\\DNSrhoN}{%10.2f}\n',rho_w);
fprintf(fileID,'\\newcommand{\\DNScpN}{%10.2f}\n',cp_w);
fprintf(fileID,'\\newcommand{\\DNSnuN}{%10.2E}\n',nu_w);
fprintf(fileID,'\\newcommand{\\DNSmuN}{%10.2E}\n',mu_w);
fprintf(fileID,'\\newcommand{\\DNSlambdaN}{%10.2f}\n',lambda_w);




fprintf(fileID,'\\newcommand{\\Qinout}{%10.2f}\n',Qin_out);
fprintf(fileID,'\\newcommand{\\Qinoutqw}{%10.2f}\n',Qin_out/Ageo_MS);
fprintf(fileID,'\\newcommand{\\QPTH}{%10.2f}\n',gradT_PT100*mdot*cpm*L);
fprintf(fileID,'\\newcommand{\\QPTHqw}{%10.2f}\n',gradT_PT100*mdot*cpm*L/Ageo_MS);
fprintf(fileID,'\\newcommand{\\qw}{%10.2f}\n',qw);
fprintf(fileID,'\\newcommand{\\qv}{%10.2f}\n',qvA);




fclose(fileID);

end 
end 
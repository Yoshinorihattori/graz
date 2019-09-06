
# coding: utf-8

# In[6]:


#calcuration libraries
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import os
import codecs
import time
#from natsort import natsorted
#get_ipython().run_line_magic('matplotlib', 'inline')


# In[7]:

def calculate_mean(data):
    s = sum(data)
    N = len(data)
    mean =s/N
    return mean

#平均からの偏差を求める
def find_difference(data):
    mean = calculate_mean(data)
    diff = []
    for num in data:
        diff.append(num-mean)
    return diff

def calculate_variance(data):
    diff = find_difference(data)
    #差の２乗を求める
    squared_diff = []
    for d in diff:
        squared_diff.append(d**2.)
    #分散を求める
    sum_squared_diff = sum(squared_diff)
    variance = sum_squared_diff/len(data)
    return variance


# Material properties for Ethylen Glycol
A = 2.0148;
B = 4.50E-3;
def cp(T):
    cp = (A + B * (T+273.15));
    return cp
def cp_K(T):
    cp_K = A + B * (T);
    cp_K =  cp_K *1000;##???
    return cp_K
C = 0.2134;
D = 6.071E-4;
def Lambda(T):
    Lambda = C + D * (T+273.15);
    return Lambda
def Lambda_K(T):
    Lambda_K = C + D * (T);
    return Lambda_K
E = 1.1001E-4;
F = 325.85;
G = -207.30;
def mu(T):
    mu = E * np.exp( F / ( (T+273.15) + G) );
    return mu
def mu_K(T):
    mu_K = E * np.exp( F / ( (T) + G) );
    return mu_K
H = 1268.28;
I = -0.66;
def rho(T):
    rho = H + I * (T+273.15);
    return rho
def rho_K(T):
    rho_K = H + I * (T);
    return rho_K
def nu(T):
    nu = mu(T)/rho(T);
    return nu
def nu_K(T):
    nu_K = mu_K(T)/rho_K(T);
    return nu
def Pr(T):
    Pr = ( mu(T) * cp(T) * 1000 ) / Lambda(T);
    return Pr
def Pr_K(T):
    Pr_K = ( mu_K(T) * cp_K(T)) / Lambda_K(T);
    return Pr

# conduction equation for inner wall temperatuer 110-115
#Lambda_K = 0.16
DeltaK = 0.065E-3
Lambda_N = 20#!!!!!!!!!!!!!!!!!!!!!!!!
Lambda_C = 0.16 #thermal conductivity[W/mK] of capton tape
Lambda_iso = 0.055 #thermal conductivity[W/mK] of glass wool, outer isolation
riso = 0.05 #thickness of isolation
#temperature distribution in pipe(ri << r << ra)
def T1(r,qv,qzu,TN,Ta):
    T1 = (qzu/(2*Lambda_N)*ra**2)*(0.5-0.5*(r/ra)**2+np.log(r/ra)) + Ta - qv/Lambda_N*(ra+DeltaK)*(Lambda_N/Lambda_C*np.log(ra/(ra+DeltaK))+np.log(r/ra))
    return T1
#temperature distribution in capton tape(ra << r << ra+DeltaK)
def T2(r,qv,Ta):
    T2 = - qv/Lambda_C*(ra+DeltaK)*np.log(r/(ra+DeltaK)) + Ta
    return T2

#Experimental Facilities
di=12E-3
da=15E-3
ri = di/2
ra = da/2
disoa = 90E-3
L=2
V=(da**2-di**2)*np.pi/4*L
Ageo_MS = di*np.pi*L #円菅内部の試験部表面積
Ugeo_MS = di*np.pi #円菅内面積
Ageo_quer = di**2*np.pi/4

#PT100 position
#x_pos_TPt100_old = 
x_pos_TPt100_new =      [0.030, 2.400, 2.600, 2.800, 3.000, 3.000, 3.140, 5.500]
x_pos_TPt100_new_tmp =  [2.400, 2.600, 2.800, 3.000, 3.140]
TPT100_new_tmp =        [0,0,0,0,0]
xTin  = 0.03
xTout = 5.50
x_MS_in = 1.2
x_MS_out = 3.2


# In[8]:


filename_list = []
T_in_list = []
T_out_list = []
Tw_list = []
Tm_list = []
mdot_list = []
qw_list = []
Re_m_list = []
Re_w_list = []
ReTau_list = []
Pr_m_list = []
Pr_w_list = []
NuMessung_list = []
cf_M_list = []
delta_NuMessung_list = []
delta_cf_M_list = []
df_list_all = []

mu_m_list = []
mu_w_list = []

I_MS_list = []

Nu_turb_Gni4Tau_list = []


# ### read measurement data 6-115

# In[9]:


where = './data_re14000_pr10_new/'

n_samples  = 0
for f_name in sorted(os.listdir(where)):
 n_samples = n_samples + 1

std_cf       = np.zeros(n_samples)
std_Re        = np.zeros(n_samples)
std_Nu        = np.zeros(n_samples)

mean_cf       = np.zeros(n_samples)
mean_Re        = np.zeros(n_samples)
mean_Nu        = np.zeros(n_samples)
nsamples   = np.zeros(n_samples)

n_samples  = 0

for fname in sorted(os.listdir(where)):
    filename = where + fname
    if filename == where + 'rubbish':
        continue
    elif filename == where + '.DS_Store':
        continue
    elif filename == where + '.DS_Store.txt':
        continue
    else:
        txt = codecs.open(filename, encoding ='cp1252')
        data = np.loadtxt(txt, skiprows = 1)
        Tsa = data[:,0] #Temperatur aussen: Ts_aussen / °C
        Tsi = data[:,1] #Temperatur innen: Ts innen / °C
        Ti = data[:,2] #Fluidtemperatur: T / °C
        m_dot_C1 = data[0,3] #Massenstrom: m_dot C1 / kg/h
        m_dot_C2 = data[0,4] #Massenstrom: m_dot C2 / kg/h
        Re_C1 = data[0,5] #Re C1 / -
        Re_C2 = data[0,6] #Re C2 / -
        v_C1 = data[0,7] #Geschwindigkeit: v C1 / m/s
        v_C2 = data[0,8] #Geschwindigkeit: v C2 / m/s
        U_MS = data[0,9] #Spannungsabfall Messstrecke: U_MS / V
        I_MS = data[0,10] #Strom Messstrecke: I_MS / A
        P_MS = data[0,11] #Leistung Messstrecke: P_MS / W
        q = data[0,12] #spzifische Leistung: q / W/m3
        P1 = data[0,13] #Druckmessung: P1 / Pa
        P2 = data[0,14] #Druckmessung: P2 / Pa
        dp = data[0,15] #Differenzdruck: dp / bar
        nu_Fluid = data[0,16] #kin. Viscositaet Fluid: nu / m2/s
        rho_Fluid = data[0,17] #Dichte Fluid: rho / kg/m3
        Lambda_Fluid = data[0,18] #Waermeleitfaehigkeit Fluid: lambda / W/mK
        cp_Fluid = data[0,19] #spez. Waermekapazitaet Fluid: cp /J/kgK
        Pr_Fluid = data[0,20] #Pr / -
        Lambda_Rohr = data[0,21] #Waermeleitfaehigkeit Rohr: lambda / W/mK
        Nu_Fluid = data[0,22] #Nu / -
        I1 = data[0,23] #Strom I / A
        U1 = data[0,24] #Spannung U / V
        dp_T1 = data[0,25] #Re_tau / - 関数名とtxtが合っていない？
        dp_T2 = data[0,26] #Re_tau_Pet / -　関数名とtxtが合っていない？
        dp_T3 = data[0,27] #dp Pet / Pa
        #Mittel (U*I) mean P / W
        Pm = data[0,28] #Mittel (U*I) mean P / W
        TPT100_a_old = data[:5,29]#PT100 T aussen(5 value)
        TPT100_i_old = data[:5,30]#PT100 T innen(5 value)
        TPT100_m_old = data[:5,31]#PT100 T mittel(5 value)
        TPT100_a_new = data[:8,38]#!!!!!!!!!!!!!
        TPT100_i_new = data[:8,39]
        TPT100_m_new = data[:8,40]
        NuPt100 = data[0,32] #Nu Pt100
        NuQm = data[0,33] #Nu Qm
        # zeta computation Einlesen
        ZRem100 = data[0,34] #zeta Re
        Zzeta100 = data[0,35] #zeta
        Zqm100 = data[0,36] #zeta qw
        TPT100_T5 = data[0,37] #PT100 T5
        
        TPT100_out = TPT100_a_new
        TPT100_in  = TPT100_i_new
        TPT100_mean= TPT100_m_new
        TPT100_new = TPT100_i_new
        Tm_aus = TPT100_new[7]
        Tm = Tm_aus

        n_samples  = n_samples  + 1
        # iteration to update the temperatuer dependent material properties --> new wall temperature 205-339
        # first loop doesnt include dq heat loss.
        for i in range(2):
            Tw = TPT100_new[6]
            T_iso = Ti[1]
            T_inf = Ti[6]

            Pr_w = Pr(Tw)
            Pr_w     = Pr(Tw)
            rho_w    = rho(Tw)
            nu_w     = nu(Tw)
            mu_w     = nu_w * rho_w
            cp_w     = cp(Tw)
            Lambda_w = Lambda(Tw)
            Pr_m     = Pr(Tm)
            rho_m    = rho(Tm)
            nu_m     = nu(Tm)
            mu_m     = nu_m * rho_m
            cp_m     = cp(Tm)
            Lambda_m = Lambda(Tm)
            
            dqhldT = 1
            zeta_mischer = 32.46
            Thl = (Tw + Tm)*0.5
            qhl = dqhldT * (Thl - T_iso)
            qhlm2 = qhl / (da*np.pi*5.45)
            Re_m = Re_C1
            Lambda_Pet = (1.8*np.log10(Re_m)-1.5)**(-2)
            
            #Specific heat capasity, cp mean
            T_in = TPT100_new[0] + 273.15
            T_out = TPT100_new[7] + 273.15
            Ac = 2.0148
            Bc = 4.50E-3
            cpm = (Ac*(T_out-T_in)+Bc/2*(T_out**2-T_in**2))/(T_out-T_in)*1000 #Eq(2.108)
            Wmean = Re_m*nu_m/di #velocisty
            mdot = m_dot_C1/3600
            
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            #x_pos_TPt100_new_tmp[0] = x_pos_TPt100_new[1]
            #x_pos_TPt100_new_tmp[1] = x_pos_TPt100_new[2]
            #x_pos_TPt100_new_tmp[2] = x_pos_TPt100_new[3]
            #x_pos_TPt100_new_tmp[3] = 0.5 * (x_pos_TPt100_new[4] + x_pos_TPt100_new[5])
            #x_pos_TPt100_new_tmp[4] = x_pos_TPt100_new[6]
            
            TPT100_new_tmp[0] = TPT100_new[1]
            TPT100_new_tmp[1] = TPT100_new[2]
            TPT100_new_tmp[2] = TPT100_new[3]
            TPT100_new_tmp[3] = 0.5 * (TPT100_new[4] + TPT100_new[5])
            TPT100_new_tmp[4] = TPT100_new[6]
            
            fit = np.polyfit(x_pos_TPt100_new_tmp, TPT100_new_tmp,1)
            gradT_PT100 = fit[0]
            Q_Pt100 = gradT_PT100*mdot*cpm*L
            qUI = P_MS/Ageo_MS
            dqw2 = (Lambda_Pet*(5.45-2.92)/di + zeta_mischer) * rho_m * Wmean**2 /2 * Wmean * Ageo_quer - qhl*(5.45-2.53)/5.45
            
            Qin_out = mdot*cpm*(T_out-T_in)
            dQ = Qin_out-P_MS
            gradT = (Qin_out/L) / (mdot * cpm )
            dTm = dqw2/(mdot*cpm)
            SA_Tm = gradT*(x_MS_out-x_pos_TPt100_new[6])
            Tm_aus = TPT100_new[7] - dTm
            Tm = Tm_aus - SA_Tm
            
            qw = Qin_out/ Ageo_MS
            if(Qin_out>0):
                qvA=qhlm2
                #print('qvA=qhlm2')
            else:
                qvA=0
                #print('qvA=0')
            
            Pr_w     = Pr(Tw)
            rho_w    = rho(Tw)
            nu_w     = nu(Tw)
            mu_w     = nu_w * rho_w
            cp_w     = cp(Tw)
            Lambda_w = Lambda(Tw)
            Pr_m     = Pr(Tm)
            rho_m    = rho(Tm)
            nu_m     = nu(Tm)
            mu_m     = nu_m * rho_m
            cp_m     = cp(Tm)
            Lambda_m = Lambda(Tm)
            
            #P_spez_m2 = P_MS/A;
            P=Pm
            P=P_MS
            #Reynolds number[-]
            Re_m = m_dot_C1/3600*4./(nu_m*rho_m*di*np.pi)
            Re_w = m_dot_C1/3600*4./(nu_w*rho_m*di*np.pi)
            
            #Wall temperature calibration
            qVol = Qin_out / V
            Tw = T1(ri,qvA,qVol,TPT100_a_new[6], TPT100_a_new[6])
            Tw_SA6 = T1(ri,qvA,qVol, TPT100_a_new[5], TPT100_a_new[5])
            Tw_SA5 = T1(ri,qvA,qVol, TPT100_a_new[4], TPT100_a_new[4])
            Tw_SA4 = T1(ri,qvA,qVol, TPT100_a_new[3], TPT100_a_new[3])
            Tw_SA3 = T1(ri,qvA,qVol, TPT100_a_new[2], TPT100_a_new[2])
            Tw_SA2 = T1(ri,qvA,qVol, TPT100_a_new[1], TPT100_a_new[1])
            Tw_SA1 = T1(ri,qvA,qVol, TPT100_a_new[0], TPT100_a_new[0])
            Tw_SA8 = T1(ri,qvA,qVol, TPT100_a_new[7], TPT100_a_new[7])
            TPT100_new[6] = Tw
            TPT100_new[0] = Tw_SA1
            TPT100_new[1] = Tw_SA2
            TPT100_new[2] = Tw_SA3
            TPT100_new[3] = Tw_SA4
            TPT100_new[4] = Tw_SA5
            TPT100_new[5] = Tw_SA6
            TPT100_new[7] = Tw_SA8
            
            Nu_QmdotCpDt_w = qw * di / (Lambda_w * (Tw-Tm ))
            NuMessung = Nu_QmdotCpDt_w
            tau_w = dp*di/(1*4)
            
            ###ADDDDDDDDDDD!!!!!!!!!
            Pr_w     = Pr(Tw)
            rho_w    = rho(Tw)
            nu_w     = nu(Tw)
            mu_w     = nu_w * rho_w
            cp_w     = cp(Tw)
            Lambda_w = Lambda(Tw)
            Pr_m     = Pr(Tm)
            rho_m    = rho(Tm)
            nu_m     = nu(Tm)
            mu_m     = nu_m * rho_m
            cp_m     = cp(Tm)
            Lambda_m = Lambda(Tm)
            
            w_tau=(tau_w/rho_m)**(1./2.) #share velocity
            ReTau = w_tau*di/nu_w
            ###ADDDDDDDDDDD!!!!!!!!!
            Wmean = Re_m*nu_m/di
            cf_M = tau_w/(rho_m*Wmean**2/2.) #Friction coefficient Eq(4.8)
        
#        print(filename)
        #桁数は四捨五入ではなく、「丸め」であることに注意
        #print('Tin','{:.4f}'.format(T_in), 'Tout','{:.4f}'.format(T_out),'Tw', '{:.4f}'.format(Tw), 'Tm', '{:.4f}'.format(Tm))
        #print('mdot', '{:.5f}'.format(mdot))
        #print('qw', '{:.5e}'.format(qw))
#        print("Re[-]",'{:.5e}'.format(Re_m),'{:.4f}'.format(ReTau))
        #print('Pr[-]','{:.2f}'.format(Pr_m),'{:.2f}'.format(Pr_w))
        #print('Nu','{:.4f}'.format(NuMessung))
#        print('cf_M','{:.10f}'.format(cf_M))
        
        # Measurement uncertainty
        #### absolute error
        T_e = 0.04#New PT100??????????????
        #mass flow rate [kg / s]
        mdot_e = 0.20E-3 * mdot
        #pressure[Pa]
        p_e = 0.35E-3 * dp #dpでいいの???????????
        #density [kg / m^3]
        rho_e = 0.66 * T_e
        #specific heat capasity [J / kg K]
        cp_e = 4.5 * T_e
        #thermal conductivity of fluid [W / m K]
        Lambda_e = 6.071e-4 * T_e
        
        #Uncertainty in each measurement influencing
        delta_mdot = (mdot_e / mdot)**2
        delta_cp = (cp_e/cp_m)**2 #cpm、cp_mどちら？
        delta_Lambda = (Lambda_e/Lambda_m)**2
        delta_T = (T_e/(T_in-Tm))**2 + ((T_e * (T_in-Tw)) / ((T_in-Tm)*(Tm-Tw)))**2 + (T_e/(Tm-Tw) )**2#T_outは影響しない？
        #Measurement uncertainty for NuMessung
        delta_NuMessung = (delta_mdot + delta_cp + delta_Lambda + delta_T)**(1./2.) * NuMessung
        #Reduction of the uncertainty of the cp calibration is the most effective, in order to reduce the uncertainty of Nu.
        
        #Uncertainty in each measurement influencing
        delta_p = (p_e / dp)**2
        delta_rho = (rho_e / rho_m)**2
        #Measurement uncertainty for cf_M
        delta_cf_M = (delta_p + delta_rho + (2.*delta_mdot))**(1./2.) * cf_M
        #Reduction of the uncertainty of the p calibration is the most effective, in order to reduce the uncertainty of cf.
        
        xi_Pet4Tau = (1.8 * np.log10(Re_m) - 1.5)**(-2)#Petukhov
        Nu_turb_Gni4Tau = ((xi_Pet4Tau/8. * Re_m * Pr_m) / (1. + 12.7 * (xi_Pet4Tau/8)**0.5 * (Pr_m**(2/3) - 1))) *(Pr_m/Pr_w)**0.11
        Nu_turb_Gni4Tau_list.append(Nu_turb_Gni4Tau)

        
        #Nu_turb_Gni4Tau_list.append(Nu_turb_Gni4Tau)
        
        filename_list.append(filename)
        T_in_list.append(T_in)
        T_out_list.append(T_out)
        Tw_list.append(Tw)
        Tm_list.append(Tm)
        mdot_list.append(mdot)
        qw_list.append(qw)
        Re_m_list.append(Re_m)
        Re_w_list.append(Re_w)
        ReTau_list.append(ReTau)
        Pr_m_list.append(Pr_m)
        Pr_w_list.append(Pr_w)
        NuMessung_list.append(NuMessung)
        cf_M_list.append(cf_M)
        delta_NuMessung_list.append(delta_NuMessung)
        delta_cf_M_list.append(delta_cf_M) 
        
        mu_m_list.append(mu_m)
        mu_w_list.append(mu_w)
        
        
        I_MS_list.append(I_MS)
        #--------------------------------------------
        
        epoch = os.path.getmtime(filename)
        filename_day = time.strftime('%d%m%Y', time.localtime(epoch))
        filename_time = time.strftime('%H%M%S', time.localtime(epoch))
        #https://tonari-it.com/python-file-get-time-epoch/#toc5
        #print(filename_day)
        #print(filename_time)

#        variance = calculate_variance(cf_M_list)
        variance = calculate_variance(NuMessung_list)
        mean     = calculate_mean    (NuMessung_list)
        std_Nu [n_samples-1] = variance**0.5
        mean_Nu[n_samples-1] = mean 
#        print(Re_m_list, std[n_samples-1])

        variance = calculate_variance(cf_M_list)
        mean     = calculate_mean    (cf_M_list)
        std_cf[n_samples-1] = variance**0.5
        mean_cf[n_samples-1] = mean 

        variance = calculate_variance(Re_m_list)
        mean     = calculate_mean    (Re_m_list)
        std_Re[n_samples-1] = variance**0.5
        mean_Re[n_samples-1] = mean 




#        print(np.size(Re_m_list),n_samples)
        nsamples[n_samples-1] = n_samples    

        
# In[45]:

fig = plt.figure(1)
ax1 = plt.subplot(131) 
plt.plot(nsamples, std_Re, color='blue',linestyle="dotted", label=r'$std^{Re}$')
plt.grid(True)
plt.legend()
plt.xlabel(r'# samples')


ax1 = plt.subplot(132) 
plt.plot(nsamples, std_cf, color='blue',linestyle="dotted", label=r'$std^{cf}$')
plt.grid(True)
plt.legend()
plt.xlabel(r'# samples')


ax1 = plt.subplot(133) 
plt.plot(nsamples, std_Nu, color='blue',linestyle="dotted", label=r'$std^{Nu}$')
plt.grid(True)
plt.legend()
plt.xlabel(r'# samples')








fig = plt.figure(3)
ax1 = plt.subplot(131) 
plt.plot(nsamples, mean_Re, color='blue',linestyle="dotted", label=r'$\bar{Re}$')
plt.grid(True)
plt.legend()
plt.xlabel(r'# samples')


ax1 = plt.subplot(132) 
plt.plot(nsamples, mean_cf, color='blue',linestyle="dotted", label=r'$\bar{cf}$')
plt.grid(True)
plt.legend()
plt.xlabel(r'# samples')


ax1 = plt.subplot(133) 
plt.plot(nsamples, mean_Nu, color='blue',linestyle="dotted", label=r'$\bar{Nu}$')
plt.grid(True)
plt.legend()
plt.xlabel(r'# samples')
#plt.title(r'$std$',fontdict=24)
#OutPutFig = './recf_pr15.pdf'
#plt.tight_layout()
#plt.savefig(OutPutFig)
#print('File name:', OutPutFig)


# In[35]:


plt.show()
# In[ ]:





# In[18]:




# In[19]:




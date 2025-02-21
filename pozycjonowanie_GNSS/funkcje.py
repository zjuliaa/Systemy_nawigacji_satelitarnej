#Julia Zapała (325710)
import numpy as np
from datetime import date
import math 

def readrnxnav(file):
    m=1
    nav=np.zeros((2000,37))
    inav=np.zeros((2000))
    n=-1
    with open(file, "r") as f:
        for s in f:                
            answer = s.find('END OF HEADER') 
            if answer != -1:
                break
        for s in f:
            s = s.replace('D', 'E')
            if m==1:
                prn=int(s2n(s,1,2))
                a = np.empty((1,6))
                a[:] = np.NaN
                a[0,0:6]=np.array(s2e(s,4,23))
            else:
                a=np.append(a,s2n(s,4,19))
            for x in range (3):
                p=23+x*19
                a=np.append(a,s2n(s,p,19))
            if m<8:
                m+=1
            else:
                n+=1
                nav[n,:]=a
                inav[n]=prn
                m=1
        nav=nav[0:n+1,:]
        inav=inav[0:n+1]
        inav = inav.astype(int)
    f.close()
    return nav, inav

def readrnxobs(file, time_start, time_end, GNSS = 'G'):
    with open(file, "r") as f: 
        for s in f:
            label = s[59:]
            if label.find('SYS / # / OBS TYPES') == 1:
                if s[0] == GNSS:
                    p = 7
                    types_header = []
                    for i in range(int(s[4:4+2])):
                        if p > 58:
                            p = 7
                            s = next(f)
                        types_header.append(s[p:p+3])
                        p += 4
                
            elif label.find('APPROX POSITION XYZ')!= -1:
                xr = np.array(([float(s[1:1+14]), float(s[15:15+14]), float(s[29:29+14])]))
            
            elif label.find('END OF HEADER') == 1:
                break
            types_of_obs = ['C1C', 'C2W']
        ind = np.zeros((len(types_header)))
        for n in range(len(types_of_obs)):
            i=(types_header.index(types_of_obs[n])) if types_of_obs[n] in types_header else -1
            if i>-1:
                ind[i]=n+1
        
        obs = np.zeros((150000, len(types_of_obs)))*np.nan
        iobs = np.zeros((150000, 3))
        n = 0
        for s in f:
            label = s[0]
            if label == '>':
                epoch = s2e(s,2,29)
                y = epoch[0]
                tt = date2tow(epoch)[1] - date2tow(epoch)[2] * 86400
                if tt > (date2tow(time_end)[1] - date2tow(epoch)[2] * 86400):
                    break
                else:
                    flag = int(round(tt))>=(date2tow(time_start)[1] - date2tow(time_start)[2] * 86400)
                if flag:
                    number_of_all_sats = int(s[33:33+2])
                    iobs[n+np.arange(0,number_of_all_sats),1] = tt
                    iobs[n+np.arange(0,number_of_all_sats),2] = date2tow(epoch)[1]
                    for sat in range(number_of_all_sats):
                        s = next(f)
                        p = 3
                        if s[0] == GNSS:
                            for i in range(len(types_header)):
                                if ind[i] != 0:
                                    obs[n+sat, int(ind[i] - 1)] = s2n(s,p,15)
                                    iobs[n+sat,0] = s2n(s,1,2)
                                p+=16
                    n += number_of_all_sats
        obs = obs[0:n, :]
        iobs = iobs[0:n,:]
        obs = np.delete(obs,iobs[:,0]==0, axis=0)
        iobs = np.delete(iobs,iobs[:,0]==0, axis=0)
        f.close()
        iobs = iobs.astype(int)
    return obs, iobs

def s2e(s,p,n):
    epoch = [int(s[p:p+4]), int(s[p+5:p+5+2]), int(s[p+8:p+8+2]), int(s[p+11:p+11+2]), int(s[p+14:p+14+2]), float(s[p+17:n])]
    return epoch 

def date2tow(data):    
    dday=date.toordinal(date(data[0],data[1],data[2])) - (date.toordinal(date(1980,1,6)))
    week = dday//7
    dow = dday%7
    tow = dow * 86400 + data[3] * 3600 + data[4] * 60 + data[5]
    return week, tow, dow

def s2n(s,p,n):
    a = s[p:p+n]
    if (not (a and not a.isspace())):
        a = np.nan
    else:
        a = float(a)        
    return a

nav_file = 'C:\Sem4\SNS\projekt2\projekt2\BRDC00WRD_R_20240650000_01D_GN.rnx'
nav, inav= readrnxnav(nav_file)

zdrowe = nav[:, 30] == 0
inav = inav[zdrowe]
nav = nav[zdrowe, :]
sat = 2

# data_obliczen = [2024, 3, 5, 11, 15, 0]
# week, tow, _ = date2tow(data_obliczen)


time_start =  [2024, 3, 5, 11, 15, 0]  
time_end =    [2024, 3, 5, 11, 15, 0] 
week, tow = date2tow(time_start)[0:2]
week_end, tow_end = date2tow(time_end)[0:2]

def xyz2blh(X, Y, Z):
    a = 6378137
    e2 = 0.00669438002290
    p = np.sqrt(X**2 + Y**2)
    phi = np.arctan2(Z, p * (1 - e2))
    phi_prev = 2 * np.pi
    N = a / (np.sqrt(1 - e2 * np.sin(phi) * np.sin(phi)))
    while abs(phi - phi_prev) > 1e-10:
        N = a / (np.sqrt(1 - e2 * np.sin(phi) * np.sin(phi)))
        h = p / np.cos(phi) - N
        phi_prev = phi
        phi = np.arctan2(Z, p * (1 - e2 * N / (N + h)))
    lam = np.arctan2(Y, X)
    h = p / np.cos(phi) - N
    phi = np.rad2deg(phi)
    lam = np.rad2deg(lam)
    return phi, lam, h

def Np(B):
    a = 6378137
    e2 = 0.00669438002290
    N = a/(1-e2*(np.sin(B)**2))**0.5
    return N

def hirvonen(X,Y,Z):
    e2 = 0.00669438002290
    r = (X**2 + Y**2)**0.5
    B = math.atan(Z/(r*(1-e2)))
    
    while 1:
        N = Np(B)
        H = r/np.cos(B) - N
        Bst = B
        B = math.atan(Z/(r*(1-(e2*(N/(N+H))))))    
        if abs(Bst-B)<(0.00001/206265):
            break
    L = math.atan2(Y,X)
    N = Np(B)
    H = r/np.cos(B) - N
    return B, L, H 


def obl_wspolrzednych_sat(sow, week, sat):
    indeks_satelity = inav==sat
    nav_satelity = nav[indeks_satelity, :]
    toe = nav_satelity[:, 17]
    roznica = np.abs(sow - toe)
    index_najmniejszej_roznicy = np.argmin(roznica)
    wiersz_nav = nav_satelity[index_najmniejszej_roznicy, :]
    e = wiersz_nav[14] 
    pierwiastek_a = wiersz_nav[16] 
    Omega_0 = wiersz_nav[19] 
    perigee = wiersz_nav[23] 
    M0 = wiersz_nav[12] 
    Toe = wiersz_nav[17] 
    i0 = wiersz_nav[21] 
    dt = (wiersz_nav[9]/1000)*np.pi/180
    af0 = wiersz_nav [6] 
    af1 = wiersz_nav[7] 
    af2 = wiersz_nav[8] 
    gps_week = wiersz_nav[27] 
    t_razem_z_tygodiami = week * 7 * 86400 + sow
    toe_razem_z_tygodniami = gps_week * 7 * 86400 + Toe
    tk = t_razem_z_tygodiami - toe_razem_z_tygodniami
    dn = wiersz_nav[11] 
    Omega_kropka = wiersz_nav[24] 
    IDOT = wiersz_nav[25] 
    Cuc = wiersz_nav[13] 
    Cus = wiersz_nav[15] 
    Cic = wiersz_nav[18] 
    Cis = wiersz_nav[20] 
    Crc = wiersz_nav[22] 
    Crs = wiersz_nav[10] 
    ni = 3.986005 * 10**14 
    omegaE = 7.2921151467 * 10**-5
    a = pierwiastek_a**2
    n0 = np.sqrt(ni/(a**3))
    n = n0 + dn
    Mk = M0 + n * tk
    Ek = Mk
    Ei=0 
    while abs(Ek-Ei)>10**-12:
        Ei = Ek
        Ek = Mk + e * math.sin(Ei)
    vk = np.arctan2(np.sqrt(1-e**2)*np.sin(Ek), np.cos(Ek) - e)
    PHI_K = vk + perigee
    duk = Cus* np.sin(2 * PHI_K) + Cuc* np.cos(2 * PHI_K)
    drk = Crs* np.sin(2 * PHI_K) + Crc* np.cos(2 * PHI_K)
    dik = Cis* np.sin(2 * PHI_K) + Cic* np.cos(2 * PHI_K) 
    uk = PHI_K + duk
    rk = a*(1 - e * np.cos(Ek)) + drk
    ik = i0 + IDOT * tk + dik
    xk = rk * np.cos(uk)
    yk = rk * np.sin(uk)
    Omega_K = Omega_0 + (Omega_kropka - omegaE) * tk - omegaE * Toe
    Xk = xk * np.cos(Omega_K) - yk*np.cos(ik) * np.sin(Omega_K)
    Yk = xk * np.sin(Omega_K) + yk*np.cos(ik) * np.cos(Omega_K)
    Zk = yk * np.sin(ik)
    c = 299792458.0 
    dts = af0 + af1 * tk + af2 * (tk**2)
    dt_rel = (-2 * np.sqrt(ni))/(c**2) * e * pierwiastek_a * np.sin(Ek)
    dt_s_rel = dts + dt_rel
    XYZ = [Xk, Yk, Zk]
    return Xk, Yk, Zk, XYZ, dt_s_rel
    


def iono (t, phi_r, lam_r, el, az):
    alfaa = [2.4214E-08,  0.0000E+00, -1.1921E-07,  5.9605E-08]
    betaa = [1.2902E+05,  0.0000E+00, -1.9661E+05, -6.5536E+04]
    c = 299792458.0
    el_sem = el / 180
    psi_sem = 0.0137/(el_sem+0.11)-0.022
    alfa_sem = [a / 180 for a in alfaa]
    phi_iip_sem = (phi_r / 180) + psi_sem * np.cos(np.deg2rad(az))
    if phi_iip_sem >  0.416666:
        phi_iip_sem = 0.416
    elif phi_iip_sem < -0.416666:
        phi_iip_sem = -0.416
    phi_ipp = phi_iip_sem * 180
    lam_iip_sem = (lam_r * 1/180) + psi_sem * np.sin(np.deg2rad(az)) / np.cos(np.deg2rad(phi_ipp))
    lam_iip = lam_iip_sem * 180
    phi_m_sem = phi_iip_sem + 0.064 * np.cos(np.deg2rad(lam_iip) - 1.617)
    t_new = 43200 * lam_iip_sem + t
    t_new = t_new % 86400
    AION = alfaa[0] + alfaa[1] * phi_m_sem + alfaa[2] * phi_m_sem**2 + alfaa[3] * phi_m_sem**3
    if AION < 0:
        AION = 0
    PION = betaa[0] + betaa[1] * phi_m_sem + betaa[2] * phi_m_sem**2 + betaa[3] * phi_m_sem**3
    if PION < 72000:
        PION = 72000
    ion_rad = (2 * np.pi * (t_new - 50400)) / PION
    mf = 1 + 16 * ((0.53 - el_sem) ** 3)
    if abs(ion_rad) <= np.pi/2:
        dIL1 = c * mf * (5 * 10 ** (-9) + AION * (1 - ion_rad**2 / 2 + ion_rad**4 / 24))
    else:
        dIL1 = c * mf * (5 * 10 ** (-9))
    return dIL1


    
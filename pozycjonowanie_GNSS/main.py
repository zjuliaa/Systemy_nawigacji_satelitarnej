import numpy as np
from funkcje import *
import math as mat
import matplotlib.pyplot as plt

nav_file = r'C:\Sem4\SNS\projekt2\projekt2\BRDC00WRD_R_20240650000_01D_GN.rnx'
obs_file = r'C:\Sem4\SNS\projekt2\projekt2\JOZ200POL_R_20240650000_01D_30S_MO.rnx'
time_start =  [2024, 3, 5, 0, 0, 0]  
time_end =    [2024, 3, 5, 23, 59, 59] 
obs, iobs = readrnxobs(obs_file, time_start, time_end, 'G')
nav, inav = readrnxnav(nav_file)
zdrowe = nav[:, 30] == 0
inav = inav[zdrowe]
nav = nav[zdrowe, :]
xr0 = [3660000.,  1400000.,  5000000.]
xyz_r_ref = [3664880.9100,  1409190.3850,  5009618.2850]
el_mask = 10 
week, tow = date2tow(time_start)[0:2]
week_end, tow_end = date2tow(time_end)[0:2]
t = tow
index_t = iobs[:, 2] == t
Pobs = obs[index_t, 0]
sats = obs[index_t, 0]
dtr = 0
dt=30
rotab= []
t0=291.15
p0=1013.25
Rh0=0.5
c1=77.64
c2=-12.96
c3=3.718*10**5
h=152
c = 299792458.0
alfaa = [2.4214E-08,  0.0000E+00, -1.1921E-07,  5.9605E-08]
betaa = [1.2902E+05,  0.0000E+00, -1.9661E+05, -6.5536E+04]
dx_val = []
dy_val = []
dz_val = []
xr0_list = []
for t in range(tow, tow_end+1, 30):
    dx_values = []
    dy_values=[]
    dz_values=[]
    index_t=iobs[:,2]==t
    Pobs= obs[index_t,0]
    sats=iobs[index_t,0]
    dtr=0
    tau=0.07
    for i in range(5):
        A = []
        y = []
        if i != 0:
            x = np.linalg.inv(A_matrix.T @ A_matrix) @ A_matrix.T @ y_wektor
            xr0 = [xr0[0] + x[0], xr0[1] + x[1], xr0[2] + x[2]]
            dtr = dtr + x[3]/c 
        for n, sat in enumerate(sats):
            if i == 0:
                tau = 0.07
            else:
                tau = rotab[n]/c
            ts=t-tau+dtr
            c = 299792458.0
            Xk, Yk, Zk, XYZ, dt_s_rel =obl_wspolrzednych_sat(ts,week,sat)
            we=7.2921151467*10**-5
            alfa=we*tau
            Rz=np.array([[np.cos(alfa), np.sin(alfa), 0], [-np.sin(alfa), np.cos(alfa), 0], [0, 0, 1]])
            xstsr= Rz@np.array([Xk, Yk, Zk])
            ro=np.sqrt((xstsr[0] - xr0[0])**2+(xstsr[1] - xr0[1])**2+(xstsr[2] - xr0[2])**2)
            if i == 0:
                rotab.append(ro)
            else: 
                rotab[n] = ro 
            phi, lam, h = xyz2blh(xr0[0], xr0[1], xr0[2])
            phi_rad, lam_rad = np.deg2rad(phi), np.deg2rad(lam)
            R = np.array([[-np.sin(phi_rad) * np.cos(lam_rad), -np.sin(lam_rad), np.cos(phi_rad) * np.cos(lam_rad)],
                                [-np.sin(phi_rad) * np.sin(lam_rad), np.cos(lam_rad), np.cos(phi_rad) * np.sin(lam_rad)],
                                [np.cos(phi_rad), 0, np.sin(phi_rad)]])
            xyz_satelity = np.array([xstsr[0], xstsr[1], xstsr[2]])
            xyz_odbiornika = np.array([xr0[0], xr0[1], xr0[2]])
            wektor_s_o = xyz_satelity - xyz_odbiornika
            xrneu= R.T@wektor_s_o
            neu= xrneu
            az= np.rad2deg(np.arctan2(neu[1], neu[0]))
            el= np.rad2deg(np.arcsin(neu[2]/np.sqrt(neu[0]**2+neu[1]**2+neu[2]**2)))
            if az < 0:
                az = az + 360
            if el > el_mask:
                p0 = 1013.25
                t0 = 291.15
                Rh0 = 0.5
                if n==0:
                    tropo=0
                    dIL1 = 0
                    dTw0=0
                    dTd0=0
                    IONO = 0
                else:
                    hort=h-31
                    p=p0*(1-0.0000226*hort)**5.225
                    temp=t0-0.0065*hort
                    Rh=Rh0*mat.exp(-0.0006396*hort)
                    Nd0=c1*(p/temp)
                    e=6.11*Rh*10**((7.5*(temp-273.15))/(temp-35.85))
                    Nw0 = c2 * (e/temp) + c3 * (e/(temp**2))
                    N0=Nd0+Nw0
                    hd= 40136+148.72*(temp-273.15)
                    dTdo=0.002277*p
                    dTwo=0.002277*((1255/temp)+0.05)*e
                    tropo=(1/np.sin(np.deg2rad(el)))*(dTdo+dTwo)
                    IONO = iono(t, xr0[1], xr0[0], el, az)
                psr= np.sqrt(((xyz_satelity[0]-xyz_odbiornika[0])**2)+((xyz_satelity[1]-xyz_odbiornika[1])**2)+((xyz_satelity[2]-xyz_odbiornika[2])**2))
                wiersz_macierzy_A= np.array([-(xyz_satelity[0]-xyz_odbiornika[0])/psr,-(xyz_satelity[1]-xyz_odbiornika[1])/psr,-(xyz_satelity[2]-xyz_odbiornika[2])/psr, 1]) 
                A.append(wiersz_macierzy_A)
                Pcal = ro - c*dt_s_rel + c*dtr
                el_wekt_wolnych = Pobs[n] - Pcal
                y.append(el_wekt_wolnych)  
        xr0_list.append(xr0)          
        A_matrix = np.array(A)
        y_wektor = np.array(y)
        dx_values.append(xr0[0] - xyz_r_ref[0])
        dy_values.append(xr0[1] - xyz_r_ref[1])
        dz_values.append(xr0[2] - xyz_r_ref[2])
        dx=dx_values[-1]
        dy=dy_values[-1]
        dz=dz_values[-1]
    dx_val.append(dx)
    dy_val.append(dy)
    dz_val.append(dz)

time = np.arange(0, len(dx_val)*30, 30) / 3600  

sigma_dx = np.std(dx_val)
sigma_dy = np.std(dy_val)
sigma_dz = np.std(dz_val)

rms_dx = np.sqrt(np.mean(np.square(dx_val)))
rms_dy = np.sqrt(np.mean(np.square(dy_val)))
rms_dz = np.sqrt(np.mean(np.square(dz_val)))
plt.figure(figsize=(12, 8))

plt.subplot(311)
plt.plot(time, dx_val, label='dx')
plt.xlabel('Time [hours]')
plt.ylabel('dx [m]')
plt.text(0.95, 0.95, f'$\\sigma$ = {sigma_dx:.2f}m\nRMS = {rms_dx:.2f}m', 
            verticalalignment='top', horizontalalignment='right', transform=plt.gca().transAxes, 
            bbox = dict(facecolor = 'white', alpha = 0.5))
plt.grid()

plt.subplot(312)
plt.plot(time, dy_val, label='dy')
plt.xlabel('Time [hours]')
plt.ylabel('dy [m]')
plt.text(0.95, 0.95, f'$\\sigma$ = {sigma_dy:.2f}m\nRMS = {rms_dy:.2f}m', 
            verticalalignment='top', horizontalalignment='right', transform=plt.gca().transAxes, 
            bbox = dict(facecolor = 'white', alpha = 0.5))
plt.grid()

plt.subplot(313)
plt.plot(time, dz_val, label='dz')
plt.xlabel('Time [hours]')
plt.ylabel('dz [m]')
plt.text(0.95, 0.95, f'$\\sigma$ = {sigma_dz:.2f}m\nRMS = {rms_dz:.2f}m', 
            verticalalignment='top', horizontalalignment='right', transform=plt.gca().transAxes, 
            bbox = dict(facecolor = 'white', alpha = 0.5))
plt.grid()

plt.tight_layout()
plt.show()

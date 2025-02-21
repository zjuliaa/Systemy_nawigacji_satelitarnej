import numpy as np
import math
from pylab import *
import matplotlib.pyplot as plt 
from matplotlib.pyplot import rc, rcParams, grid 
import math 
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pylab import *
import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk
from tkcalendar import DateEntry
import cartopy.crs as ccrs
from matplotlib import rc
from matplotlib.figure import Figure
from tkinter import messagebox

canvas1 = None
canvas2 = None
canvas3 = None
canvas4 = None
canvas6 = None
canvas7 = None

def oblicz():
    global canvas1, canvas2, canvas3, canvas4, canvas6, canvas7
    for canvas in [canvas1, canvas2, canvas3, canvas4, canvas6, canvas7]:
        if canvas is not None:
            canvas.get_tk_widget().destroy()

    def read_yuma(almanac_file):
        if almanac_file:
            alm = open(almanac_file)
            alm_lines = alm.readlines()
            all_sat = []
            for idx, value in enumerate(alm_lines):
                if value[0:3]=='ID:':
                    one_sat_block = alm_lines[idx:idx+13]
                    one_sat = []
                    for line in one_sat_block:
                        data = line.split(':')
                        one_sat.append(float(data[1].strip()))
                    all_sat.append(one_sat)
            alm.close()
            all_sat = np.array(all_sat)
            return (all_sat)
        
    def read_alm(file):
        m = 0
        with open(file, "r") as f:
            block = []
            nav_data = []
            for s in f:                
                if m<13:
                    m+=1
                    block.append(s)
                else:
                    block_array = np.genfromtxt(block,delimiter=10).T
                    nav_data.extend(block_array)
                    m = 0
                    block = []
        nav_data = np.array(nav_data)        
        return nav_data

    def get_prn_number(nav_data):
        prns = []
        for nav in nav_data:
            nsat = nav[0]
            if 0<nsat<=37:
                prn = int(nsat)
                prns.append(prn)
            elif 38<=nsat<=64:
                prn = 100 + int(nsat-37)
                prns.append(prn)
            elif 111<=nsat<=118:
                prn = 400 + int(nsat-110)
                prns.append(prn)
            elif 201<=nsat<=263:
                prn = 200 + int(nsat-200)
                prns.append(prn)    
            elif 264<=nsat<=310:
                prn = 300 + int(nsat-263)
                prns.append(prn)
            elif 311<=nsat:
                prn = 300 + int(nsat-328)
                prns.append(prn)           
            else: 
                prn = 500 + int(nsat)
                prns.append(prn)
        return prns

    def get_alm_data(file):
        nav_data = read_alm(file)
        prns = get_prn_number(nav_data)
        nav_data[:,0] = prns
        return nav_data

    def create_prn_alm(nav_data):
        prns = []
        for nav in nav_data:
            nsat = nav[0]
            if 0<nsat<=37:
                prn = 'G'+str(int(nsat)).zfill(2)
                prns.append(prn)
            elif 38<=nsat<=64:
                prn = 'R'+str(int(nsat-37)).zfill(2)
                prns.append(prn)
            elif 111<=nsat<=118:
                prn = 'Q'+str(int(nsat-110)).zfill(2)
                prns.append(prn)
            elif 201<=nsat<=263:
                prn = 'E'+str(int(nsat-200)).zfill(2)
                prns.append(prn)    
            elif 264<=nsat<=310:
                prn = 'C'+str(int(nsat-263)).zfill(2)
                prns.append(prn)
            elif 311<=nsat:
                prn = 'C'+str(int(nsat-328)).zfill(2)
                prns.append(prn)           
            else: 
                prn = 'S'+str(int(nsat)).zfill(2)
                prns.append(prn)
        return prns

    def get_alm_data_str(alm_file):
        alm_data = read_alm(alm_file)
        selected_system = selected_system_var.get()
        if selected_system == "GPS":
            alm_data = alm_data[alm_data[:, 0] <= 32, :]
            title = 'GPS'
        elif selected_system == "Galileo":
            alm_data = alm_data[(alm_data[:, 0] >= 201) & (alm_data[:, 0] <= 236), :]
            title = 'Galileo'
        elif selected_system == "GLONASS":
            alm_data = alm_data[(alm_data[:, 0] >= 38) & (alm_data[:, 0] <= 61), :]
            title = 'GLONASS'
        prns = create_prn_alm(alm_data) 
        return alm_data, prns, title   

    nav, prns, title = get_alm_data_str('Almanac2024053.alm')
    niesatelity = nav[:, 0] >399
    prns = create_prn_alm(nav)
    satelity = nav[:, 0] < 400
    nav = nav[satelity, :]

    def julday(y,m,d,h = 0):
        if m <= 2:
            y = y - 1
            m = m + 12
        jd = math.floor(365.25*(y+4716))+math.floor(30.6001*(m+1))+d+h/24-1537.5
        return jd

    def get_gps_time(y, m, d, h = 0, mnt = 0, s = 0):
        days = julday(y, m, d) - julday(1980, 1, 6)
        week = days // 7
        day = days%7
        sow = day * 86400 + h * 3600 + mnt * 60 + s
        return int(week), sow

    wiersz_nav = nav[0,:]

    def obliczenie_wspolrzednych_satelity(wiersz_nav, y, m, d, h = 0, mnt = 0, s=0):
        week, sec_of_week = get_gps_time(y,m,d,h,mnt,s)
        actual_sec = week * 7 * 86400 + sec_of_week
        mi = 3.986005 * 10**14
        OMGE = 7.2921151467 * 10**-5
        svprn = wiersz_nav[0]
        health = wiersz_nav[1]
        e = wiersz_nav[2]
        toa = wiersz_nav[7]
        i = (54 + wiersz_nav[8]) * np.pi/180
        Omega_dot = (wiersz_nav[9]/1000) * np.pi/180
        sqrtA = wiersz_nav[3]
        Omega = wiersz_nav[4] * np.pi/180
        omega = wiersz_nav[5] * np.pi/180
        M0 = wiersz_nav[6] * np.pi/180
        af0 = wiersz_nav[10]
        af1 = wiersz_nav[11]
        gps_week = wiersz_nav[12]
        actual_nav = gps_week * 7 * 86400 + toa
        tk = actual_sec - actual_nav
        a = sqrtA**2
        n = math.sqrt(mi/a**3)
        M_k = M0 + n * tk
        E_k = M_k
        E_i = 0
        while abs(E_k - E_i) > 10**-12:
            E_i = E_k 
            E_k = M_k + e * np.sin(E_i)
        vk = np.arctan2(np.sqrt(1 - e**2) * np.sin(E_k), np.cos(E_k) - e)
        uk = vk + omega
        rk = a * (1 - e * math.cos(E_k))
        Omega_k = Omega + (Omega_dot - OMGE) * tk - OMGE* toa
        xk = rk * np.cos(uk)
        yk = rk * np.sin(uk)
        X = xk * np.cos(Omega_k) - yk * np.cos(i) * np.sin(Omega_k)
        Y = xk * np.sin(Omega_k) + yk * np.cos(i) * np.cos(Omega_k)
        Z = yk * np.sin(i)
        return np.array([X, Y, Z])

    def blh2xyz(phi, lam, h):
        phi_rad = np.deg2rad(phi)
        lam_rad = np.deg2rad(lam)
        a = 6378137
        e2 = 0.00669438002290
        N = a / (np.sqrt(1 - e2 * np.sin(phi_rad) * np.sin(phi_rad)))
        X= (N + h) * np.cos(phi_rad) * np.cos(lam_rad)
        Y = (N + h) * np.cos(phi_rad) * np.sin(lam_rad)
        Z = (N * (1 - e2) + h) * np.sin(phi_rad)
        return np.array([X, Y, Z])

    maska = maska_obserwacji_combobox.get()
    maska = float(maska)
    year, m, d = date_entry.get_date().year, date_entry.get_date().month, date_entry.get_date().day
    
    a = 6378137
    e2 = 0.00669438002290
    phi = phi_entry.get()
    lam = lambda_entry.get()
    try:
        phi_val = float(phi_entry.get())
        lam_val = float(lambda_entry.get())
        if not -90 <= phi_val <= 90:
            messagebox.showerror("Błąd", "Wartość phi poza zakresem! Wprowadź wartość od -90 do 90.")
            return
        if not -180 <= lam_val <= 180:
            messagebox.showerror("Błąd", "Wartość lambda poza zakresem! Wprowadź wartość od -180 do 180.")
            return
    except ValueError:
        messagebox.showerror("Błąd", "Proszę wprowadzić poprawne wartości numeryczne.")

    h= h_entry.get()
    phi = float(phi)
    lam = float(lam)
    h = float(h)
    xyz_odbiornika = blh2xyz(phi, lam, h) 
    phi_rad = np.deg2rad(phi)
    lam_rad = np.deg2rad(lam)
    hours_in_day = 24
    minutes_in_hour = 60
    seconds_in_minute = 60
    R = np.array([[-np.sin(phi_rad) * np.cos(lam_rad), -np.sin(lam_rad), np.cos(phi_rad) * np.cos(lam_rad)],
                                [-np.sin(phi_rad) * np.sin(lam_rad), np.cos(lam_rad), np.cos(phi_rad) * np.sin(lam_rad)],
                                [np.cos(phi_rad), 0, np.sin(phi_rad)]])
    satelity_wspolrzedne = []
    A_macierz= []
    elewacje = [[] for _ in range(len(nav))] 
    azymuty = [[] for _ in range(len(nav))]
    xyz_sat = [[] for _ in range(len(nav))]

    for hour in range(hours_in_day):
        for minute in range(minutes_in_hour):
            A_minuta = []
            for index, sat in enumerate(nav):
                xyz_satelity = obliczenie_wspolrzednych_satelity(sat, year,m,d,hour,minute,0)
                xyz_sat[index].append(xyz_satelity)
                wektor_s_o = xyz_satelity - xyz_odbiornika
                xrneu= R.T@wektor_s_o
                dlugosc_s_o = np.linalg.norm(wektor_s_o)
                dlugosc = np.linalg.norm(xrneu)
                neu= xrneu
                az= np.rad2deg(np.arctan2(neu[1], neu[0]))
                el= np.rad2deg(np.arcsin(neu[2]/np.sqrt(neu[0]**2+neu[1]**2+neu[2]**2)))
                elewacje[index].append(el)
                if az < 0:
                    az = az + 360
                azymuty[index].append(az)
                if el > maska:
                    psr= np.sqrt(((xyz_satelity[0]-xyz_odbiornika[0])**2)+((xyz_satelity[1]-xyz_odbiornika[1])**2)+((xyz_satelity[2]-xyz_odbiornika[2])**2))
                    wiersz_macierzy_A= np.array([-(xyz_satelity[0]-xyz_odbiornika[0])/psr,-(xyz_satelity[1]-xyz_odbiornika[1])/psr,-(xyz_satelity[2]-xyz_odbiornika[2])/psr, 1]) 
                    A_minuta.append(wiersz_macierzy_A)
                satelity_wspolrzedne.append(xyz_satelity)
            A_macierz.append(A_minuta)

    prns_with_index = [prns[index] for index in range(len(prns))]
    elewacje = np.array(elewacje)
    satelity_wspolrzedne = np.array(satelity_wspolrzedne)
    xyz_sat = np.array(xyz_sat)

    # wykres elewacji
    fig1, ax1 = plt.subplots(figsize=(12, 6))
    for idx, sat_elevations in enumerate(elewacje):
        filtered_elevations = [el if el > maska else None for el in sat_elevations]
        plt.plot(filtered_elevations, label = prns_with_index[idx])
    plt.xticks(ticks=[i*60 for i in range(24)], labels=[f'{i}:00' for i in range(24)], rotation=45)  
    plt.xlabel('Godzina dnia')
    plt.ylabel('Elewacja [°]')
    plt.title(f'Elewacja Satelitów {title}')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='x-small') 
    plt.subplots_adjust(right=0.75)
    plt.grid(True)
    plt.tight_layout()
    canvas1 = FigureCanvasTkAgg(fig1, master=tab1)
    canvas1.draw()
    canvas1.get_tk_widget().pack(fill='both', expand=True)

    # wykres widoczności satelitów
    visible_satellites_per_minute = [np.sum(np.array(elewacje)[:, minute] > maska) for minute in range(1440)]
    minutes = np.arange(1440)
    hours = minutes / 60  
    fig2, ax2 = plt.subplots(figsize=(15, 7))  
    ax2.plot(minutes, visible_satellites_per_minute)
    ax2.fill_between(minutes, visible_satellites_per_minute, step="pre", color='green', alpha=0.3)
    ax2.set_xlabel('Godzina dnia')
    xticks = np.arange(0, 1441, 60) 
    xticklabels = [f'{int(xtick/60)}:00' for xtick in xticks]  
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xticklabels, rotation=45) 
    ax2.set_ylabel('Liczba widocznych satelitów')
    ax2.set_ylim(bottom=0)  
    ax2.set_title(f'Liczba widocznych satelitów {title}')
    ax2.grid(True)
    plt.tight_layout()
    canvas2 = FigureCanvasTkAgg(fig2, master=tab2)
    canvas2.draw()
    canvas2.get_tk_widget().pack(fill='both', expand=True)

    # wykres widoczności satelitów w czasie
    visibility = np.array(elewacje) > maska  
    fig3, ax3 = plt.subplots(figsize=(15, 7)) 
    for sat_idx, sat_visibility in enumerate(visibility):
        for minute in range(1440):
            if sat_visibility[minute]:  
                ax3.fill_between([minute, minute+1], [sat_idx - 0.4, sat_idx - 0.4], [sat_idx + 0.4, sat_idx + 0.4],  color=plt.cm.tab20(sat_idx % 20), step="pre", linewidth=0)
    ax3.set_yticks(range(len(visibility)))
    ax3.set_yticklabels(prns) 
    ax3.set_ylim(-1, len(visibility))
    ax3.set_xlim(0, 1440)
    xticks = list(range(0, 1440, 60)) 
    xticklabels = [f"{i // 60:02d}:{i % 60:02d}" for i in xticks]  
    ax3.set_xticks(xticks) 
    ax3.set_xticklabels(xticklabels, rotation=45)  
    ax3.set_xlabel('Godzina dnia')
    ax3.set_ylabel('Satelita')
    ax3.set_title(f'Widoczności satelitów {title}')
    ax3.grid(True)
    plt.tight_layout()
    canvas3 = FigureCanvasTkAgg(fig3, master=tab3)
    canvas3.draw()
    canvas3.get_tk_widget().pack(fill='both', expand=True)

    # wykres DOP
    GDOPs = []
    PDOPs = []
    TDOPs = []
    HDOPs = []
    VDOPs = []
    for A in A_macierz:
        A_transpose = np.transpose(A)
        ATA = np.dot(A_transpose, A)
        Q = np.linalg.inv(ATA)
        GDOP = np.sqrt(Q[0, 0] + Q[1, 1] + Q[2, 2] + Q[3, 3])
        PDOP = np.sqrt(Q[0, 0] + Q[1, 1] + Q[2, 2])
        TDOP = np.sqrt(Q[3, 3])
        HDOP = np.sqrt(Q[0, 0] + Q[1, 1])  
        VDOP = np.sqrt(Q[2, 2])  
        GDOPs.append(GDOP)
        PDOPs.append(PDOP)
        TDOPs.append(TDOP)
        HDOPs.append(HDOP)
        VDOPs.append(VDOP)
    times = np.linspace(0, 24, 1440) 
    fig4, ax4 = plt.subplots(figsize=(15, 7)) 
    plt.plot(times, GDOPs, label='GDOP', linewidth=2)
    plt.plot(times, PDOPs, label='PDOP', linewidth=2)
    plt.plot(times, TDOPs, label='TDOP', linewidth=2)
    plt.plot(times, HDOPs, label='HDOP', linewidth=2)
    plt.plot(times, VDOPs, label='VDOP', linewidth=2)
    plt.legend()
    plt.xlabel('Godzina dnia')
    plt.ylabel('Wartość DOP')
    plt.xticks(np.arange(0, 25, 1), [f'{i}:00' for i in range(25)], rotation=45)
    ax4.set_title(f'Wartości DOP {title}')
    plt.grid(True)
    canvas4 = FigureCanvasTkAgg(fig4, master=tab4)
    canvas4.draw()
    canvas4.get_tk_widget().pack(fill='both', expand=True)

    elewacje_rad = [np.radians(90 - np.array(el)) for el in elewacje]
    azymuty_rad = [np.radians(az) for az in azymuty]
    wybrana_godzina = hour_cb.get()
    wybrana_minuta = minute_cb.get()
    godzina = int(wybrana_godzina)
    minuta = int(wybrana_minuta)
    aktualny_moment = godzina * 60 + minuta

    # wykres skyplot dla okreslonej godziny
    elewacje_radiany = [np.radians(el) for el in elewacje]
    elewacje_radiany = np.array(elewacje_radiany)
    elewacje_t = elewacje_radiany.T
    fig6 = plt.figure()
    ax6 = fig6.add_subplot(111, polar=True)
    maska_rad = np.radians(maska)
    widoczne_satelity = []
    for idx, elewacja_satelity in enumerate(elewacje_t[aktualny_moment]):
        if elewacja_satelity >= maska_rad:
            widoczne_satelity.append(idx)  
    for idx in widoczne_satelity:
        ax6.plot(azymuty_rad[idx], elewacje_rad[idx], label=prns_with_index[idx], alpha=0.5)
        ax6.plot(azymuty_rad[idx][aktualny_moment], elewacje_rad[idx][aktualny_moment], 'o')
        ax6.annotate(prns_with_index[idx], xy=(azymuty_rad[idx][aktualny_moment], elewacje_rad[idx][aktualny_moment]), xytext=(5, 5), textcoords='offset points', ha='right', va='bottom')
    ax6.set_theta_zero_location('N')  
    ax6.set_theta_direction(-1)  
    ax6.set_ylim(0, np.radians(90))  
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    ax6.set_title(f'Skyplot {title} o {wybrana_godzina}:{wybrana_minuta}')
    ax6.set_xticks(np.deg2rad([0, 45, 90, 135, 180, 225, 270, 315]))
    ax6.set_xticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'])
    ax6.set_yticks(np.radians([0, 15, 30, 45, 60, 75, 90]))
    ax6.set_yticklabels(['90°', '75°', '60°', '45°', '30°', '15°', 'Horyzont'])
    canvas6 = FigureCanvasTkAgg(fig6, master=tab6)
    canvas6.draw()
    canvas6.get_tk_widget().pack(fill='both', expand=True)
    
    #Wykres Groundtrack
    def latlon(XYZ):
        r_delta = np.linalg.norm(XYZ[:2])
        sinA = XYZ[1] / r_delta
        cosA = XYZ[0] / r_delta
        Lon = math.atan2(sinA, cosA)
        if Lon < -math.pi:
            Lon = Lon + 2 * math.pi
        Lat = math.asin(XYZ[2] / np.linalg.norm(XYZ))
        return math.degrees(Lat), math.degrees(Lon)

    rc('grid', color='gray', linewidth=0.5, linestyle='--')
    fontsize = 20
    rc('xtick', labelsize=fontsize)
    rc('ytick', labelsize=fontsize)
    rc('font', size=fontsize)
    fig7 = plt.figure(figsize=(14, 7))
    ax7 = fig7.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    plt.subplots_adjust(right=0.75)  
    ax7.coastlines()
    ax7.set_global()
    ax7.grid(True)
    ax7.set_title(f'Współrzędne satelitów {title} o {wybrana_godzina}:{wybrana_minuta}')
    plt.yticks(range(-90, 91, 30))
    plt.xticks(range(-180, 181, 30))
    for idx, (elewacja_trajektoria, xyz_trajektoria) in enumerate(zip(elewacje, xyz_sat)):
        latitudes, longitudes = [], []
        segment_started = False
        for moment, elewacja in enumerate(elewacja_trajektoria):
            if np.radians(elewacja) >= maska_rad:
                if not segment_started:
                    segment_started = True
                    segment_start_idx = moment
                lat, lon = latlon(xyz_trajektoria[moment])
                latitudes.append(lat)
                longitudes.append(lon)
            else:
                if segment_started:
                    ax7.plot(longitudes[segment_start_idx:moment], latitudes[segment_start_idx:moment], '-', transform=ccrs.Geodetic(), label=prns_with_index[idx] if segment_start_idx == 0 else "", alpha=0.5)
                    latitudes, longitudes = [], []
                    segment_started = False
        if segment_started:
            ax7.plot(longitudes, latitudes, '-', transform=ccrs.Geodetic(), label=prns_with_index[idx] if segment_start_idx == 0 else "", alpha=0.5)
        if np.radians(elewacja_trajektoria[aktualny_moment]) >= maska_rad:
            lat, lon = latlon(xyz_trajektoria[aktualny_moment])
            ax7.plot(lon, lat, 'o', color='red', transform=ccrs.Geodetic())
            ax7.text(lon, lat, prns_with_index[idx], fontsize=12, color='black', ha='center', va='center', transform=ccrs.Geodetic())
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    canvas7 = FigureCanvasTkAgg(fig7, master=tab7)
    canvas7.draw()
    canvas7.get_tk_widget().pack(fill='both', expand=True)
    
#Interfejs graficzny
def create_datetime_entry(frame, label, row):
    tk.Label(frame, text=label).grid(row=row, column=0, sticky=tk.W, padx=10, pady=5)
    date_entry = DateEntry(frame, width=12, background='darkblue', foreground='white', borderwidth=2)
    date_entry.grid(row=row, column=1, padx=10, pady=5, sticky=tk.W)
    return date_entry

root = tk.Tk()
root.title("Zadanie 1")
notebook = ttk.Notebook(root)
notebook.pack(fill=tk.BOTH, expand=True)
frame_inputs = tk.Frame(notebook)
notebook.add(frame_inputs, text="Dane")

default_phi = "52"  
default_lambda = "21"  
default_h = "100"
phi_entry = tk.Entry(frame_inputs)
phi_entry.grid(row=2, column=1, padx=10, pady=5)
phi_entry.insert(0, default_phi)  
lambda_entry = tk.Entry(frame_inputs)
lambda_entry.grid(row=3, column=1, padx=10, pady=5)
lambda_entry.insert(0, default_lambda)  
h_entry = tk.Entry(frame_inputs)
h_entry.grid(row=4, column=1, padx=10, pady=5)
h_entry.insert(0, default_h)
label3 = tk.Label(frame_inputs, text="Współrzędne miejsca obserwacji (phi, lam, h):", font=("Arial", 12))
label3.grid(row=2, column=0, padx=10, pady=10, sticky=tk.W)
date_entry = create_datetime_entry(frame_inputs, "Data:", 0)
label_maska_obserwacji = tk.Label(frame_inputs, text="Maska obserwacji [°]:")
label_maska_obserwacji.grid(row=6, column=0, sticky=tk.W, padx=10, pady=5)
mask_values = ['5', '10', '15'] 
maska_obserwacji_combobox = ttk.Combobox(frame_inputs, values=mask_values, state='readonly', width=15)
maska_obserwacji_combobox.grid(row=6, column=1, padx=10, pady=5)
maska_obserwacji_combobox.set('5')
selected_system_var = tk.StringVar(value="GPS")
gps_radiobutton = tk.Radiobutton(frame_inputs, text="GPS", variable=selected_system_var, value="GPS")
galileo_radiobutton = tk.Radiobutton(frame_inputs, text="Galileo", variable=selected_system_var, value="Galileo")
glonass_radiobutton = tk.Radiobutton(frame_inputs, text="GLONASS", variable=selected_system_var, value="GLONASS")
gps_radiobutton.grid(row=8, column=1, sticky=tk.W, padx=10, pady=5)
galileo_radiobutton.grid(row=8, column=2, sticky=tk.W, padx=10, pady=5)
glonass_radiobutton.grid(row=8, column=3, sticky=tk.W, padx=10, pady=5)
label4 = tk.Label(frame_inputs, text="Wybierz epokę dla wykresu Skyplot i Groundtrack:", font=("Arial", 12))
label4.grid(row=9, column=0, padx=10, pady=10, sticky=tk.W)
hour_cb = ttk.Combobox(frame_inputs, values=[f"{i:02d}" for i in range(24)], width=5)
hour_cb.grid(row=9, column=1, padx=10, pady=10)
hour_cb.set("00")
minute_cb = ttk.Combobox(frame_inputs, values=[f"{i:02d}" for i in range(0, 60, 5)], width=5)  
minute_cb.grid(row=9, column=2, padx=10, pady=10)
minute_cb.set("00")
button1 = tk.Button(frame_inputs, text="Zamknij program", width=30, command=root.destroy)
button1.grid(row=10, column=0, columnspan=2, padx=10, pady=10)
button2 = tk.Button(frame_inputs, text="Oblicz", width=30, command=oblicz)
button2.grid(row=11, column=0, columnspan=2, padx=10, pady=10)

tab1 = ttk.Frame(notebook)
notebook.add(tab1, text="Wykres elewacji")
tab2 = ttk.Frame(notebook)
notebook.add(tab2, text="Wykres widoczności satelitów")
tab3 = ttk.Frame(notebook)
notebook.add(tab3, text="Wykres widoczności satelitów w czasie")
tab4 = ttk.Frame(notebook)
notebook.add(tab4, text="Wykres DOP")
tab6 = ttk.Frame(notebook)
notebook.add(tab6, text="Wykres skyplot")
tab7 = ttk.Frame(notebook)
notebook.add(tab7, text="Wykres Groundtrack")
notebook.pack(expand=True, fill="both")
root.mainloop()

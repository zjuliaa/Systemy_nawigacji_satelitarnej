import numpy as np
import matplotlib.pyplot as plt
import math

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

plik = r'C:\Sem4\SNS\projekt3\MIMA.pos'
wyniki= np.genfromtxt(plik, comments= '%')
xyzref = [3655333.847, 1403901.067, 5018038.047]
dxyz = wyniki[:,2:2+3] - xyzref
ratio = wyniki[:, -1]
time = np.arange(0, len(dxyz)*30, 30)/3600

phi, lam, h = hirvonen(xyzref[0], xyzref[1], xyzref[2])
Rneu = np.array([[-np.sin(phi) * np.cos(lam), -np.sin(lam), np.cos(phi) * np.cos(lam)],
            [-np.sin(phi) * np.sin(lam), np.cos(lam), np.cos(phi) * np.sin(lam)],
            [np.cos(phi), 0, np.sin(phi)]])
dneu = []
for dx in dxyz:
    dn = Rneu.T@dx
    dneu.append(dn)
dneu = np.array(dneu)

fig, ax = plt.subplots(3, 1, figsize=(10, 8))
for i in range(3):
    ax[i].plot(time, dneu[:, i])
    ax[i].set_ylabel(['N (m)', 'E (m)', 'U (m)'][i])

ax[2].set_xlabel('Time (hours)')
plt.suptitle('Różnica współrzędnych z referencyjnymi NEU')
plt.tight_layout(rect=[0, 0, 1, 0.95])

fix = wyniki[:, 5] ==1
floats = wyniki[:,5]==2
tfix= time[fix]
tfloat = time[floats]
time_fix = time[fix]
time_float = time[floats]
fig, ax = plt.subplots(3, 1, figsize=(10, 8))
for i in range(3):
    ax[i].plot(time_fix, dneu[fix, i], color='red', label='Fixed')
    ax[i].plot(time_float, dneu[floats, i], color='green', label='Float')
    ax[i].set_ylabel(['N (m)', 'E (m)', 'U (m)'][i])
    ax[i].legend()

ax[2].set_xlabel('Time (hours)')
plt.suptitle('Liczba poszczególnych typów roziązań wektora')
plt.tight_layout(rect=[0, 0, 1, 0.95])

fig, ax = plt.subplots()
ax.scatter(dneu[floats, 1], dneu[floats, 0], color='green', label='Float')
ax.scatter(dneu[fix, 1], dneu[fix, 0], color='red', label='Fixed')
ax.axhline(0)
ax.axvline(0)
ax.axis('equal')
ax.legend()
plt.title('Rozkład punktów')

std = np.std(dneu[fix,:], axis=0)
fig, ax = plt.subplots()
ax.bar(['N', 'E', 'U'], std)
plt.xlabel('Współrzędne')
plt.title('Odchylenie standardowe')

rms_dx = np.sqrt(np.mean(np.square(dneu[fix,0])))
rms_dy = np.sqrt(np.mean(np.square(dneu[fix,1])))
rms_dz = np.sqrt(np.mean(np.square(dneu[fix,2])))
rms_values = [rms_dx, rms_dy, rms_dz]
coordinates = ['N', 'E', 'U']
fig, ax = plt.subplots()
plt.bar(coordinates, rms_values)
plt.title('RMS dla poszczególnych współrzędnych')
plt.xlabel('Współrzędne')

fig, ax = plt.subplots(1, 1)
plt.title('Wartość testu ratio')
ax.plot(time, ratio)
plt.xlabel('Czas [h]')
plt.grid()

max_errors = np.max(np.abs(dneu), axis=0)
fig, ax = plt.subplots()
ax.bar(['N', 'E', 'U'], max_errors)
plt.title('Maksymalne wartości błędów')
plt.xlabel('Współrzędne')

plt.show()

# Legge le coordinate di un profilo alare (file .dat o .csv) — tipicamente coppie (x,y) ordinate lungo il bordo superiore e inferiore.
# Calcola proprietà geometriche fondamentali:
# area (superficie del profilo proiettata),
# baricentro (centroide) rispetto al bordo d’attacco,
# spessore massimo e posizione,
# corda (chord length) e Mean Aerodynamic Chord (MAC, semplice approssimazione),
# linea di curvatura / camber line e camber massimo,
# momento d’inerzia (secondo momento d’area) attorno all’asse x (utile per stime strutturali).
# Visualizza il profilo, la camber line, il punto di spessore massimo e il centroide con grafici.


import numpy as np
import matplotlib.pyplot as plt
import os as os
import tkinter as tk
from tkinter import filedialog, messagebox


# scelgo il file da aprire
def selezione_file_finestra():
    """Apri una finestra per scegliere il file. Ritorna percorso (stringa) o '' se annullato."""
    root = tk.Tk()
    root.withdraw()  # nasconde la finestra principale
    percorso = filedialog.askopenfilename(
        title="Seleziona il file con le coordinate (x,y) del profilo NACA",
        filetypes=[("Text/CSV", "*.txt *.dat *.csv"), ("Tutti i file", "*.*")]
    )
    root.destroy()
    return percorso or ""


#leggo coordinate da file 
def leggi_dati(percorso):
    coordinates = []
    with open(percorso,'r') as f:
        for a in f:
            a = a.strip() #strip mi toglie spazi bianchi da inizio e fine di una stringa
            if not a or a.startswith('#'): continue
            parts = a.replace(',',' ').split()
            if len(parts)>=2:
                try:
                    x = float(parts[0])
                    y = float(parts[1])
                    coordinates.append((x,y))
                except: pass
    return np.array(coordinates)


def Gauss_Green_area_MomentiInerzia(coordinates):
    x = coordinates[:,0]
    y = coordinates[:,1]
    # if not (np.allclose(x[0],x[-1]) and np.allclose(y[0],y[-1])):
    #     x = np.append(x,x[0])
    #     y = np.append(y,y[0])
    cross = x[:-1]*y[1:] - x[1:]*y[:-1] 
    A = 0.5*np.sum(cross)
    Cx = (1/(6*A))*np.sum((x[:-1]+x[1:])*cross)
    Cy = (1/(6*A))*np.sum((y[:-1]+y[1:])*cross)
    Baricentro = np.array([Cx,Cy])
    I_xx = (1.0/12.0)*np.sum((y[:-1]**2 + y[:-1]*y[1:] + y[1:]**2)*cross)
    I_yy = (1.0/12.0)*np.sum((x[:-1]**2 + x[:-1]*y[1:] + y[1:]**2)*cross)
    return abs(A), x, Baricentro, I_xx, I_yy

#trovo TE e LE oltre a spessori: 

# trovo LE (leading edge) e TE (trailing edge) 
    # coordinates è un array (N,2) con colonne [x, y]
def TE_LE_spessore_camber (coordinates):
    if coordinates.size == 0:
        TE = LE = None
    else:
        x_vals = coordinates[:, 0]
        y_vals = coordinates[:, 1]
    # Leading edge: minimo di x (punto più vicino al bordo d'attacco)
        le_idx = int(np.argmin(x_vals))
        LE = coordinates[le_idx]
    #Trailing edge: massimo di x (punto più lontano, tipicamente x~=1 per profili NACA)
        te_idx = int(np.argmax(x_vals))
        TE = coordinates[te_idx]
    #Semi spessore positivo e negativo: escludendo il punto di LE (x min)
    x_vals_sup = x_vals[:le_idx]
    y_vals_sup = y_vals[:le_idx]
    x_vals_inf = x_vals[le_idx+1:]
    y_vals_inf = y_vals[le_idx+1:]

    #aumento sensibilità di calcolo con griglia nell'intervallo di x e identifico gli estremi
    x_min = max(x_vals_sup.min(),x_vals_inf.min())
    x_max = min(x_vals_sup.max(),x_vals_inf.max())
    x_griglia = np.linspace(x_min,x_max,1200)

    ordine_sup = np.argsort(x_vals_sup); 
    x_vals_sup_ord, y_vals_sup_ord = x_vals_sup[ordine_sup], y_vals_sup[ordine_sup]
    ordine_inf = np.argsort(x_vals_inf); 
    x_vals_inf_ord, y_vals_inf_ord = x_vals_inf[ordine_inf], y_vals_inf[ordine_inf]

    #interpolo le curve superiori e inferiori sulla griglia di punti
    y_vals_sup_i = np.interp(x_griglia, x_vals_sup_ord, y_vals_sup_ord)
    y_vals_inf_i = np.interp(x_griglia, x_vals_inf_ord, y_vals_inf_ord)
    #spessore punto per punto
    spessore = y_vals_sup_i - y_vals_inf_i
    # print(spessore)
    spessore = np.abs(spessore)
    idx_max = int(np.argmax(spessore))
    print("Max thickness:", format(spessore[idx_max],".5f"), "at x =", format(x_griglia[idx_max],".5f"))

    #calcolo il camber del profilo:
    c_x = 0.5*(y_vals_sup_i + y_vals_inf_i)
    camber_max = np.max(np.abs(c_x))
    idx_camber_max = int(np.argmax(camber_max))

    plt.figure(figsize=(12,3))
    plt.plot(x_vals_sup, y_vals_sup, '-o', label='upper')
    plt.plot(x_vals_inf, y_vals_inf, '-o', label='lower')
    plt.axvline(coordinates[le_idx,0], color='k', linestyle='--')
    plt.legend()
    plt.show(block = False)
    return {
        "TE": TE,
        "LE": LE,
        "spessore": spessore,
        "idx_max": idx_max,
        "x_griglia": x_griglia,
        "y_vals_sup_i": y_vals_sup_i,
        "y_vals_inf_i": y_vals_inf_i,
        "camber_max":camber_max,
        "idx_camber_max":idx_camber_max,        
    }


percorso = selezione_file_finestra()
coordinates = leggi_dati(percorso)
if not np.allclose(coordinates[0], coordinates[-1]):
                coordinates = np.vstack([coordinates, coordinates[0]])
Area, x, Baricentro, I_xx, I_yy = Gauss_Green_area_MomentiInerzia(coordinates)
result = TE_LE_spessore_camber(coordinates)
calcoli = dict(Area = Area, Baricentro = Baricentro, coordinata_x = x.reshape(-1,2), spessore = result["spessore"],
                x_spessore_max = result["x_griglia"][result["idx_max"]], camber_max = result["camber_max"],
                x_camber_max  = result["x_griglia"][result["idx_camber_max"]], I_xx = I_xx, I_yy = I_yy)
# con reshape(-1,2) calcola righe automatiche

#visualizzo i risultati 
for q,w in calcoli.items():
     print(q,w,sep=" = ")

#plotting del profilo con dati in input:
plt.figure(figsize=[12,3])
plt.plot(coordinates[:,0], coordinates[:,1], '-r')
plt.grid(True)
plt.show()
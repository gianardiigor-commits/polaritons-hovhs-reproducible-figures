# Adapted by Igor Gianardi (first author) from code originally written by
# Michele Pini (second author; ORCID: https://orcid.org/0000-0001-5522-5109).

#%% Initialization cell
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib
import numpy as np
import os

dir_path = os.path.abspath(os.path.dirname(__file__))
os.chdir(dir_path)

#Set fonts
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'

#font di "Latex"  Prova a vedere quali sono quelli disponibili


#Definining function to import float numbers from txt files
def get_float(txt):    #get_float è una funzione di Python        
    try:
        return float(txt) #prova a convertire il txt in un float
    except Exception:
        if (len(txt) > 0): print('note: could not parse',txt,'as float')
        return np.NaN
    
#Definining function to import float numbers from txt files and converts D->E (needed for fortran double numbers)
def get_float_DtoE(txt): #
    try:
        return float(txt.decode().replace('D', 'E'))  #Prima di convertire in float rimpiazza D con E per leggere bene i dati txt generati da Fortran
    except Exception:
        if (len(txt) > 0): print('note: could not parse',txt,'as float')
        return np.NaN
     
    
        if (len(txt) > 0): print('note: could not parse',txt,'as string')
        return np.NaN    

#Definining function to import integer numbers from txt files
def get_int(txt):
    try:
        return int(txt)
    except Exception:
        if (len(txt) > 0): print('note: could not parse',txt,'as integer')
        return np.NAN 


# In the following the letter q does not denote the phtonic momentum but the light-matter coupling g
#%% Input read cell - A(q,omega)
#It reads input data from txt files

nq = 300
nOmega = 400

#Path of the file to read

path1 = 'DplotCheq1.txt'
print(path1)

#Empty arrays needed afterwards
Aqw_flip = []
q_a = []

n_col=3 #number of columns to read
n_skip_rows=3 #number of rows to skip at the beginning of the txt file 

#Array of converters from txt to float
cnv = dict(zip(range(n_col), [get_float_DtoE] * n_col))
#this is equivalent to cnv={0:get_float_DtoE,1:get_float_DtoE,2:get_float_DtoE}
#It defines the function to convert from the strings in the txt file to floats numbers


for i in range(0, nq, 1):
    file1_in = open(path1, 'rt')
    q_tmp_a, Omega_a, Aqw_a = \
       np.loadtxt(file1_in, unpack=True, skiprows=n_skip_rows + i * nOmega,
                  usecols=range(0, n_col), max_rows=nOmega, converters=cnv)
   # converte da g/W → g/Eg (W/Eg = 0.3)
    q_tmp_a *= 0.3
    Aqw_flip_tmp = np.flip(Aqw_a, axis=0)
    Aqw_flip.append(Aqw_flip_tmp)
    q_a.append(q_tmp_a[0])

Aqw_hm = np.array(Aqw_flip).transpose()

qmin_tot = q_a[0]
qmax_tot = q_a[nq - 1]
Omega_tot_min = Omega_a[0]
Omega_tot_max = Omega_a[nOmega - 1]

print(qmin_tot)
print(qmax_tot)
print(Omega_tot_min)
print(Omega_tot_max)

#%% Renormalizing data

Omega_min = 0.75
Omega_max = 1.2

min_Aqw = 0
max_Aqw = 20

for i in range(0, nOmega):
    if Omega_a[i] >= Omega_min:
        imin = max(i - 1, 0)
        print(imin)
        print(Omega_a[imin])
        break

for i in range(0, nOmega):
    if Omega_a[i] >= Omega_max or i == nOmega - 1:
        imax = i
        print(imax)
        print(Omega_a[imax])
        break

Aqw_hm_rnm = Aqw_hm.copy()
for i in range(0, nq):
    for j in range(0, nOmega):
        if Aqw_hm[j, i] > max_Aqw:
            Aqw_hm_rnm[j, i] = max_Aqw
        if Aqw_hm[j, i] < min_Aqw:
            Aqw_hm_rnm[j, i] = min_Aqw

Aqw_hm_rnm_cut = Aqw_hm_rnm[nOmega - imax - 1:nOmega - imin][:]

dq = (q_a[1] - q_a[0])
qmin = q_a[0] - dq / 2
qmax = q_a[nq - 1] + dq / 2

print(qmin)
print(qmax)

#%% heatmap plot cell
AxesLabelSize = 18
TitleSize = 18
TicksSize = 18
LegendSize = 16
LabelSize = 18

xtick_major_spacing1 = 0.05     # passo più adatto al nuovo xmax
xtick_minor_spacing1 = 0.01
ytick_major_spacing_rep = 0.05
ytick_minor_spacing_rep = 0.01
ticks_pad = 4
title_pad = 10

xmin = 0
xmax = 0.15
ymin = 0.9
ymax = 1.1

ratio = 0.15/ (ymax - ymin)

fig = plt.figure()
fig, ax = plt.subplots(1, 1, sharex=True, sharey=False, figsize=(6, 6))

ax.set_ylim([ymin, ymax])
ax.yaxis.set_tick_params(which='both', direction='out', left='on', pad=ticks_pad, labelsize=TicksSize)
ax.yaxis.set_major_locator(ticker.MultipleLocator(ytick_major_spacing_rep))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(ytick_minor_spacing_rep))
ax.set_xlim([xmin, xmax])
ax.xaxis.set_major_locator(ticker.MultipleLocator(xtick_major_spacing1))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(xtick_minor_spacing1))
ax.xaxis.set_tick_params(which='both', direction='out', bottom='on', pad=ticks_pad, labelsize=TicksSize)

ax.set_xlabel(r'$g/E_g$', fontsize=AxesLabelSize, labelpad=1)
ax.set_ylabel(r'$\omega/E_g$', fontsize=AxesLabelSize, labelpad=0)

cax = ax.imshow(Aqw_hm_rnm_cut, interpolation='none', cmap='inferno',
                extent=[qmin, qmax, Omega_min, Omega_max], aspect=ratio)

cbar = plt.colorbar(cax, fraction=0.03, pad=0.04, ax=ax)
cbar.ax.tick_params(labelsize=LegendSize)

ax.set_title(r'$A(\omega) E_g$', fontsize=TitleSize, pad=title_pad)


# Colore associato a intensità 20 nella mappa 'inferno'
color_at_20 = matplotlib.colormaps['inferno'](20 / max_Aqw)


# Etichetta omega_P(g) in colore mappa, posizione assoluta
ax.text(0.045, 0.950, r'$\omega_{\mathrm{P}}(g)$',
        fontsize=LabelSize, color=color_at_20, ha='left', va='top')

#Etichetta bianca su tre righe, centrata in g/W=0.25 e appena sopra ω/Eg = 1
ax.text(0.12, 1.025, 'Interband\np-h continuum',
        fontsize=LegendSize, color='white', ha='center', va='bottom')


#Etichetta bianca su tre righe, centrata in g/W=0.25 e appena sopra ω/Eg = 1
ax.text(0.12, 0.960, 'Zero\nabsorption',
      fontsize=LegendSize, color='white', ha='center', va='bottom')

ax.axhline(y=1.0008, color='white', linestyle='--', linewidth=1.4)

#ax.text(0.02, ymax - 0.01, r'$\omega_0 = E_g$',
#       fontsize=LabelSize, color='white', ha='left', va='top')

ax.text(
    0.01, ymax - 0.01,
    r'$\quad \,\,\,\,\omega_0 = E_g$'      # 1st math fragment
    '\n'                     # real newline (not raw)
    r'$W/E_g = 0.3$',        # 2nd math fragment
    fontsize=LabelSize,
    color='white',
    ha='left',
    va='top',
    bbox=dict(
        boxstyle='round,pad=0.2',
        edgecolor='white',
        facecolor='none',
        linewidth=1.2
    )
)

                              

# Lettura del file (2 colonne: x = g/W, y = DelR)
#file_path = 'DelRdA_P0.4350_data.txt'
#g_data, DelR_data = np.loadtxt(file_path, unpack=True)

# Applica le trasformazioni: -y + 1
#DelR_transformed = -DelR_data + 1

# Filtro dei dati fino a g/W <= 0.5
#mask = g_data <= 0.5
#g_filtered = g_data[mask]
#DelR_filtered = DelR_transformed[mask]

# Plot sullo stesso ax
#ax.plot(g_filtered, DelR_filtered, color=color_at_20, linestyle='--', linewidth=1.4,
#       label=r'Modified $\Delta_R$: $1 - \Delta_R$')

plt.savefig('Aqw_vs_gtoGap_final.eps', format='eps', bbox_inches='tight', pad_inches=0.05)
plt.savefig('Aqw_vs_gtoGap_final.pdf', format='pdf', bbox_inches='tight', pad_inches=0.05)


#%%
a=[1,1]
b=a
b[1]=2
print(a)

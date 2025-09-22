# -*- coding: utf-8 -*-
#%% Initialization cell
import matplotlib.pyplot as plt # plot library
import matplotlib.ticker as ticker # ticks library
import numpy as np # numerical library
import os

dir_path=os.path.abspath(os.path.dirname(__file__))
os.chdir(dir_path)
path = os.getcwd()
print(path)

# Set fonts
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'

# Definining function to import float numbers from txt files
def get_float(txt):
    try:
        return float(txt)
    except Exception:
        if len(txt) > 0:
            print('note: could not parse', txt, 'as float')
        return np.NaN

#%% Input read cell

# Array of polarizations 
Pa = [0.435, 0.50, 0.70]
nP = len(Pa)

# Empty arrays to fill with data for both polarizations
omega_a_P = []
DOS_up_a_P = []
DOS_down_a_P = []

# Number of columns to read
n_col = 3

# Array of converters from txt to float
cnv = dict(zip(range(n_col), [get_float] * n_col))  # Usare get_float per tutti i dati

# Importing data
for i in range(0, nP):
    P = Pa[i]
    # Setting path to txt file for polarization P (P is written as a 4 decimal float)
    path = f'Dos_P{P:.4f}_data.txt'
    print(path)
    # Opening file
    file_in = open(path, 'rt')
    # Importing the columns as temporary arrays
    omega_a_tmp, DOS_down_a_tmp, DOS_up_a_tmp = \
    np.loadtxt(file_in, unpack=True, skiprows=3, usecols=range(n_col), converters=cnv)
    # Appending the temporary arrays to the arrays that contain all the polarizations
    omega_a_P.append(omega_a_tmp)
    DOS_up_a_P.append(DOS_up_a_tmp)
    DOS_down_a_P.append(DOS_down_a_tmp)
    # Closing file
    file_in.close()

#%% plot cell

# Font sizes
AxesLabelSize=20
TitleSize=20
TicksSize=19
LegendSize=18
LabelSize=20

# List of colors (for more colors, check https://matplotlib.org/stable/gallery/color/named_colors.html)
color_list = ['Blue', 'BlueViolet', 'Red']

# List of dash styles (linea continua e tratteggiata)
dash_list = [[1, 0], [4, 2], [1, 1.5]]  # [1,0] rappresenta una linea continua, [4,2] una linea tratteggiata

# Linewidth
lw = 1.6

# Generating figure
fig = plt.figure()

# Generating subplots (2,1 per due subplot verticali)
fig, axs = plt.subplots(2, 1, sharex=True, sharey=False, figsize=(6, 10))

# Adjust the space between the two vertical plots
fig.subplots_adjust(hspace=0.1)

# Setting boundaries of x and y axes
xmin = 0.81
xmax = 1.19
ymin1 = -1.5
ymax1 = 0.2
ymin2 = -1.5
ymax2 = 1.5

# Setting spacing between major and minor ticks on both x and y axes
xtick_major_spacing = 0.05
xtick_minor_spacing = 0.01
ytick_major_spacing = 0.5
ytick_minor_spacing = 0.1

ytick_major_spacing_lower = 0.5
ytick_minor_spacing_lower = 0.1


# Setting the distance of numbers from the plot
ticks_pad = 7

# Coordinates of textbox for subplot label (a, b)
x_textbox = 0.90
y_textbox = 0.94

# Coordinates of legends
lgloc1 = (0.06, 0.35)  # upper plot
lgloc2 = (0.06, 0.6)   # lower plot

# Legends settings
lg_handle_size = 28 / LegendSize  # Size of the handles 
lg_label_sp = 0.3  # Vertical spacing
lg_borderpad = 0.2  # Border padding
lg_handletextpad = 0.5  # Distance between handles and text

# Setting different limits for the y-axes
axs[0].set_ylim(ymin1, ymax1)  # First subplot: ymin1=-8, ymax1=1
axs[1].set_ylim(ymin2, ymax2)  # Second subplot: ymin2=-5, ymax2=5

# Setting common x-axis ticks and subplot-specific y-axis ticks
for i, ax in enumerate(axs.flat):
    ax.set_xlim(xmin, xmax)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(xtick_major_spacing))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(xtick_minor_spacing))

    if i == 0:  # upper subplot
        ax.yaxis.set_major_locator(ticker.MultipleLocator(ytick_major_spacing))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(ytick_minor_spacing))
    else:  # lower subplot
        ax.yaxis.set_major_locator(ticker.MultipleLocator(ytick_major_spacing_lower))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(ytick_minor_spacing_lower))

    ax.xaxis.set_tick_params(which='both', direction='in', top='on', pad=ticks_pad, labelsize=TicksSize)
    ax.yaxis.set_tick_params(which='both', direction='in', right='on', pad=ticks_pad, labelsize=TicksSize)


# Upper subplot
ax = axs[0]
# Plotting data
labels = ["2D Check.", "1D Parab.", "2D Parab." ]
for i in range(0, nP):
    P = Pa[i]
    ax.plot(omega_a_P[i], DOS_up_a_P[i], dashes=dash_list[i], color=color_list[i], linewidth=lw, label=labels[i])
    
# Setting axis labels
ax.set_ylabel(
    r'$\operatorname{Im}\Sigma^{R}_{\mathrm{ph}}'
    r'\left/\left(g/W\right)^{2}E_g \right.$',
    fontsize=AxesLabelSize
)

# Adding textbox with subplot label (a)
ax.text(x_textbox, y_textbox, r'$(a)$', transform=ax.transAxes, fontsize=LabelSize,
        horizontalalignment='left', verticalalignment='center')

# Plotting legend for the upper subplot
ax.legend(loc=lgloc1, fontsize=LegendSize, frameon=0, handlelength=lg_handle_size,
          labelspacing=lg_label_sp, borderpad=lg_borderpad, title=r'', title_fontsize=LegendSize,
          handletextpad=lg_handletextpad)

# Lower subplot
ax = axs[1]
# Plotting data
labels = ["2D Check.", "1D Parab.", "2D Parab." ]
for i in range(0, nP):
    P = Pa[i]
    ax.plot(omega_a_P[i], DOS_down_a_P[i], dashes=dash_list[i], color=color_list[i], linewidth=lw, label=labels[i])

# Setting axis labels
ax.set_ylabel(
    r'$\operatorname{Re}\Sigma^{R}_{\mathrm{ph}}'
    r'\left/\left(g/W\right)^{2}E_g \right.$',
    fontsize=AxesLabelSize
)

ax.set_xlabel(r'$\omega/E_g$', fontsize=AxesLabelSize)




# Adding textbox with subplot label (b)
ax.text(x_textbox, y_textbox, r'$(b)$', transform=ax.transAxes, fontsize=LabelSize,
        horizontalalignment='left', verticalalignment='center')

# Plotting legend for the lower subplot
ax.legend(loc=lgloc2, fontsize=LegendSize, frameon=0, handlelength=lg_handle_size,
          labelspacing=lg_label_sp, borderpad=lg_borderpad, title=r'', title_fontsize=LegendSize,
          handletextpad=lg_handletextpad)
          
axs[1].text(0.956, 0.05, '$W/E_g= 0.3$', 
            transform=axs[1].transAxes,
            fontsize=LabelSize,
            color='black',
            ha='right',
            va='bottom',
            bbox=dict(boxstyle='round,pad=0.3',
                      edgecolor='black',
                      facecolor='none',
                      linewidth=1.2))


# Exporting eps file
plt.savefig('Figure-Self_Final.eps', format='eps', bbox_inches='tight', pad_inches=0.05)

# Exporting PDF file
plt.savefig('Figure-Self_Final.pdf', format='pdf', bbox_inches='tight', pad_inches=0.05)



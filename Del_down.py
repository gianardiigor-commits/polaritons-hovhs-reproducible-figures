#%% Initialization and input read cell
import matplotlib.pyplot as plt #plot library
import matplotlib.ticker as ticker #ticks library
import numpy as np #numerical library
import os

# Set the working directory to the current file's location
dir_path = os.path.abspath(os.path.dirname(__file__))
os.chdir(dir_path)
path = os.getcwd()
print(path)

# Set fonts
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'

# Defining function to import float numbers from txt files
def get_float(txt):
    try:
        return float(txt)
    except Exception:
        if (len(txt) > 0): print('note: could not parse', txt, 'as float')
        return np.NaN

# Array of polarizations (reduced to 4)
Pa = [0.435, 0.50, 0.70]
nP = len(Pa)

# Empty arrays to fill with data for all the polarizations
omega_a_P = []
DelR_a_P = []

# Number of columns to read (2: one for frequency and one for DelR)
n_col = 2
# Array of converters from txt to float
cnv = dict(zip(range(n_col), [get_float]*n_col))

# Importing data
for i in range(0, nP):
    P = Pa[i]
    # Setting path to txt file for polarization P (P is written as a 4 decimal float)
    path = f'DelRdA_gtoGap_P{P:.4f}_data.txt'
    print(path)
    # Opening file
    with open(path, 'rt') as file_in:
        # Importing the columns as temporary arrays (frequency and DelR)
        omega_a_tmp, DelR_a_tmp = np.loadtxt(file_in, unpack=True, skiprows=3, usecols=range(n_col), converters=cnv)
        # Appending the temporary arrays to the arrays that contain all the polarizations
        omega_a_P.append(omega_a_tmp)
        DelR_a_P.append(DelR_a_tmp)

#%% Plot cell

# Font sizes
AxesLabelSize = 18
TitleSize = 20
TicksSize = 19
LegendSize = 16
LabelSize = 16

# List of colors (for more colors, check https://matplotlib.org/stable/gallery/color/named_colors.html)
color_list = ['Blue', 'BlueViolet', 'Red']
# List of dash styles (defined as [filled space,empty space,filled space,...])
dash_list = [[1, 0], [4, 2], [1, 1.5]]
# Linewidth
lw = 1.6

# Generating figure
fig = plt.figure()
# Creating a single plot (fig, ax) instead of subplots
fig, ax = plt.subplots(figsize=(6, 6))

# Setting boundaries of x and y axes
xmin = 0
xmax = 0.5
ymin = 0
ymax = 0.5

# Setting spacing between major and minor ticks on both x and y axes
xtick_major_spacing = 0.1
xtick_minor_spacing = 0.05
ytick_major_spacing = 0.1
ytick_minor_spacing = 0.02

# Setting the distance of numbers from the plot (needed to avoid superpositions on the bottom left)
ticks_pad = 7

# This sets the limits for the plot
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

# Setting the ticks
ax.xaxis.set_major_locator(ticker.MultipleLocator(xtick_major_spacing))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(xtick_minor_spacing))
ax.yaxis.set_major_locator(ticker.MultipleLocator(ytick_major_spacing))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(ytick_minor_spacing))
ax.xaxis.set_tick_params(which='both', direction='in', top='on', pad=ticks_pad, labelsize=TicksSize)
ax.yaxis.set_tick_params(which='both', direction='in', right='on', pad=ticks_pad, labelsize=TicksSize)

# Plotting data for the 3 polarizations with custom labels
labels = ["2D Check. (HOVHS)", "1D Parab.", "2D Parab."]
for i in range(0, nP):
    P = Pa[i]
    ax.plot(omega_a_P[i], DelR_a_P[i], dashes=dash_list[i], color=color_list[i], label=labels[i], linewidth=lw)

# Adding mathematical expressions as annotations
# ax.text(0.6, 0.019, r'$\sim\left(\frac{g}{W}\right)^{2} \log{\left(\frac{W}{g}\right)}$', 
#        fontsize=LabelSize, va='center', ha='left', color='black', fontweight='bold')
# ax.text(0.6, 0.055, r'$\sim\left(\frac{g}{W}\right)^{\frac{4}{3}} $', 
#        fontsize=LabelSize, va='center', ha='left', color='black', fontweight='bold')
# ax.text(0.6, 0.139, r'$\sim\left(\frac{g}{W}\right)^{\frac{4}{3}} \log^{\frac{2}{3}}\!\left(\frac{W}{g}\right)$', 
#        fontsize=LabelSize, va='center', ha='left', color='black', fontweight='bold')


# Setting axis labels
ax.set_ylabel(r'$\delta \omega_{\mathrm{P}} / \omega_0$', fontsize=AxesLabelSize)
ax.set_xlabel(r'$g / E_g$', fontsize=AxesLabelSize)


# Coordinates of legend
lgloc = (0.07, 0.45)

# Plotting legend
lg_handle_size = 28 / LegendSize  # Size of the handles
lg_label_sp = 0.3  # Vertical spacing
lg_borderpad = 0.2  # Border padding
lg_handletextpad = 0.5  # Distance between handles and text

ax.legend(loc=lgloc, fontsize=LegendSize, frameon=0, handlelength=lg_handle_size,
          labelspacing=lg_label_sp, borderpad=lg_borderpad, handletextpad=lg_handletextpad)

ax.text(
    0.04, ymax - 0.02,
    r'$\quad \,\,\,\,\,\omega_0 = E_g$'      # 1st math fragment
    '\n'                     # real newline (not raw)
    r'$W/E_g = 0.3$',        # 2nd math fragment
    fontsize=LabelSize,
    color='black',
    ha='left',
    va='top',
    bbox=dict(
        boxstyle='round,pad=0.2',
        edgecolor='black',
        facecolor='none',
        linewidth=1.2
    )
)
          
#Legend con Frame:

#legend = ax.legend(
#    loc=lgloc,
#    fontsize=LegendSize,
#  frameon=True,
#   handlelength=lg_handle_size,
#  labelspacing=lg_label_sp,
#   borderpad=lg_borderpad,
#    handletextpad=lg_handletextpad
#)

# Setting black border with custom thickness
#legend.get_frame().set_edgecolor('black')
#legend.get_frame().set_linewidth(1.2)  
        
# Exporting EPS file
plt.savefig('Figure-Del_Down_Final_gtoGap_NoLeg_4.09.25.eps', format='eps', bbox_inches='tight', pad_inches=0.05)

# Exporting PDF file
plt.savefig('Figure-Del_Down_Final_gtoGap_NoLeg_4.09.25.pdf', format='pdf', bbox_inches='tight', pad_inches=0.05)

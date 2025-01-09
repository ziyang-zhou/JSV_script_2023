import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.signal import find_peaks
import os
import pandas as pd
import pdb
import sys
sys.path.append('/home/ziyz1701/storage/CD_airfoil/tool/analysis_folder')  # Add the parent directory to the path
from functions import analysis

#############################################################################################
#                                        Load Inputs                                        #
#############################################################################################
# Load the settings dataframe:
settings = pd.read_csv("../setting.csv", index_col= 0)
# Unpack Relevant Settings:
project_folder = settings.at["project_folder", settings.columns[0]]
save_folder = settings.at["save_folder", settings.columns[0]]
folder_list = eval(settings.at["folder_list", settings.columns[0]])
name_list = eval(settings.at["name_list", settings.columns[0]])
profile_color = eval(settings.at["profile_color", settings.columns[0]])
linestyle_list = eval(settings.at["linestyle_list", settings.columns[0]])
variable = 'mean_cf'
# Physical parameters:
M_exp = eval(settings.at["M_exp", settings.columns[0]])
M_PF = eval(settings.at["M_PF", settings.columns[0]])
U_0 = eval(settings.at["U_0", settings.columns[0]])
#signal analysis
chord = eval(settings.at["chord", settings.columns[0]])

linestyle_list_b = eval(settings.at["linestyle_list_b", settings.columns[0]])

#plot setting
linewidth_3D_result = eval(settings.at["linewidth_3D_result", settings.columns[0]])
linewidth_2D_result = eval(settings.at["linewidth_2D_result", settings.columns[0]])
linewidth_list = eval(settings.at["linewidth_list", settings.columns[0]])
marker_list = eval(settings.at["marker_list", settings.columns[0]])
marker_size = eval(settings.at["marker_size", settings.columns[0]])
fontsize_var = eval(settings.at["fontsize_var", settings.columns[0]])
probe_number_list = eval(settings.at["probe_number_list", settings.columns[0]])
plot_size_square = eval(settings.at["plot_size_square", settings.columns[0]])
plot_size_rectangle = eval(settings.at["plot_size_rectangle", settings.columns[0]])
name_list = eval(settings.at["name_list", settings.columns[0]])
profile_color = eval(settings.at["profile_color", settings.columns[0]])

#############################################################################################
#                                        Calculate profile                      #
#############################################################################################

#This script correlates the signal in mic_data to the variable found in the probes in probe_path

save_path = 'aoa_correlation/'
read_path = 'aoa_correlation/'

nb_chunk = 1 #was 12
tstart = 0.2 # initial t to start from
tend = 0.9
def next_greater_power_of_2(x):  
    return 2**(x-1).bit_length()

#Load signal 2 #########################################################################################################
# Path to the folder containing CSV files
folder_path = project_folder + read_path + 'temporal_plane_jet_020_010.csv'

data = np.loadtxt(folder_path, delimiter=',', skiprows=1)
time = data[:,0]
x_velocity = data[:,3]
y_velocity = data[:,4]

keep_time = (time>tstart) & (time<tend)
y_velocity = y_velocity[keep_time] #construct the second signal which is to be correlated against
x_velocity = x_velocity[keep_time]
sig_2 = analysis.get_aoa(x_velocity,y_velocity)#reconstruct the aoa

peaks,_ = find_peaks(sig_2, height=0.75, width = 50)
troughs,_ = find_peaks(-sig_2, height=-0.8, width = 50)

figure_scale_factor = 1.3
plt.figure(figsize=(plot_size_rectangle[0]//figure_scale_factor,plot_size_rectangle[1]//figure_scale_factor))
# Enable LaTeX rendering
plt.rc('text', usetex=True)
# Use plain TeX font (Computer Modern) for legend
plt.rc('legend', fontsize=fontsize_var, fancybox=True, framealpha=0.7, edgecolor='black', facecolor='white', borderaxespad=0.5, labelspacing=0.5)
plt.rc('font', family='serif', weight='normal')

plt.tick_params(axis='both', which='major', labelsize=fontsize_var)
ax = plt.gca()
ax.yaxis.offsetText.set(size=fontsize_var)
plt.yticks(fontsize=fontsize_var)
plt.xticks(fontsize=fontsize_var)



plt.plot(time[keep_time],sig_2,label=r'$\alpha_i$',linewidth=linewidth_list[0])
plt.plot(time[keep_time][peaks], sig_2[peaks], 'rx', label='Peak',markersize=marker_size)
plt.plot(time[keep_time][troughs], sig_2[troughs], 'bo', label='Trough',markersize=marker_size)
plt.tight_layout()
plt.xlabel('Time [s]', fontsize=fontsize_var, fontweight='bold')
plt.ylabel(r'$\alpha_i$', fontsize=fontsize_var, fontweight='bold')
plt.ylim(0.0, 4.5)
#plt.xlim(0.05, 0.25)
plt.grid(alpha=0.2, color='black')
plt.legend(prop = { "size": fontsize_var*0.7})

# Fill the area between x-axis values 0.38 and 0.51 with a shade of red
fill_start = 0.38
fill_end = 0.51
plt.fill_betweenx(y=np.linspace(0, 4.5, 100), x1=fill_start, x2=fill_end, color='red', alpha=0.3, label='Shaded Area')
plt.tight_layout()
#plt.savefig(project_folder + save_path + 'angle_of_attack.jpg', format='jpeg', dpi=300)
plt.savefig(project_folder + save_path + 'angle_of_attack.pdf', dpi=300)

plt.close()

print('Peaks at {}'.format(time[keep_time][peaks]))
print('Troughs at {}'.format(time[keep_time][troughs]))

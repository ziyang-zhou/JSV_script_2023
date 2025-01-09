import pandas as pd
import numpy as np
from scipy import signal as sig
import matplotlib.pyplot as plt
import sys
import pdb
# Import the ScalarFormatter class
from matplotlib.ticker import ScalarFormatter
sys.path.append('/home/ziyz1701/storage/CD_airfoil/tool/analysis_folder')  # Add the parent directory to the path
from functions import analysis, futils
#############################################################################################
#                                        Load Inputs                                        #
#############################################################################################
# Load the settings dataframe:
settings = pd.read_csv("setting.csv", index_col= 0)
# Unpack Relevant Settings:
project_folder = settings.at["project_folder", settings.columns[0]]
save_folder = settings.at["save_folder", settings.columns[0]]
folder_list = eval(settings.at["folder_list", settings.columns[0]])
name_list = eval(settings.at["name_list", settings.columns[0]])
profile_color = eval(settings.at["profile_color", settings.columns[0]])
linestyle_list = eval(settings.at["linestyle_list", settings.columns[0]])

# Physical parameters:
M_exp = eval(settings.at["M_exp", settings.columns[0]])
M_PF = eval(settings.at["M_PF", settings.columns[0]])
U_0 = eval(settings.at["U_0", settings.columns[0]])
#signal analysis
chord = eval(settings.at["chord", settings.columns[0]])

#plot setting
linewidth_3D_result = eval(settings.at["linewidth_3D_result", settings.columns[0]])
linewidth_2D_result = eval(settings.at["linewidth_2D_result", settings.columns[0]])
linewidth_list = eval(settings.at["linewidth_list", settings.columns[0]])
marker_list = eval(settings.at["marker_list", settings.columns[0]])
marker_size = eval(settings.at["marker_size", settings.columns[0]])
fontsize_var = eval(settings.at["fontsize_var", settings.columns[0]])
plot_size_square = eval(settings.at["plot_size_square", settings.columns[0]])
plot_size_rectangle = eval(settings.at["plot_size_rectangle", settings.columns[0]])

# Load data
minus_cp_data = np.genfromtxt(project_folder + 'mean_cp/' + 'mean_cp_shear_layer_refined.csv', delimiter=',')
minus_cp_data_2 = np.genfromtxt(project_folder + 'mean_cp/'+'mean_cp_shear_layer_refined_with_trip.csv', delimiter=',')
minus_cp_2011_data = np.genfromtxt(project_folder + 'mean_cp/'+'2011-mean-cp-extracted.csv', delimiter=',')
minus_cp_2d_ff_data = np.genfromtxt(project_folder + 'mean_cp/'+'mean_cp_2d_freefield.csv', delimiter=',')
minus_cp_2d_lips_data = np.genfromtxt(project_folder + 'mean_cp/'+'mean_cp_2D_lips.csv', delimiter=',')
minus_cp_expt = np.genfromtxt(project_folder + 'mean_cp/'+'Cp_8deg_ECL.txt', delimiter='\t\t')
minus_cp_expt2 = np.genfromtxt(project_folder + 'mean_cp/'+'Cp_8deg_UdS.txt', delimiter='\t\t\t\t\t\t')

minus_cp_data_3 = np.genfromtxt(project_folder + 'mean_cp/'+'mean_cp_reference_DNS_damping_zone_corrected.csv', delimiter=',')
#############################################################################################
#                                        Calculate profile                      #
#############################################################################################

# Extract data columns
minus_cp = minus_cp_data[:, 0:]
minus_cp_2 = minus_cp_data_2[:, 0:]
minus_cp_2011 = minus_cp_2011_data[:, 0:]
minus_cp_2d_lips = minus_cp_2d_lips_data[:, 0:]
minus_cp_2d_ff = minus_cp_2d_ff_data[:, 0:]

minus_cp_3 = minus_cp_data_3[:, 0:]

# Coordinate adjustment
minus_cp[:, 0] -= 0.15
minus_cp_2[:, 0] -= 0.15
minus_cp_2011[:, 0] -= 0.15
minus_cp_2d_lips[:, 0] -= 0.15
minus_cp_2d_ff[:, 0] -= 0.15

minus_cp_3[:,0] -= 0.15

# Non-dimensionalize
minus_cp[:, 1:] *= -1
minus_cp_2[:, 1:] *= -1
minus_cp_3[:, 1:] *= -1
minus_cp_2011[:, 1:] *= -1
minus_cp_2d_lips[:, 1:] *= -1
minus_cp_2d_ff[:, 1:] *= -1

minus_cp[:, 0] /= 0.1356
minus_cp_2[:, 0] /= 0.1356
minus_cp_3[:, 0] /= 0.1356
minus_cp_2011[:, 0] /= 0.1356
minus_cp_2d_lips[:, 0] /= 0.1356
minus_cp_2d_ff[:, 0] /= 0.1356

plt.figure(figsize=plot_size_rectangle)
# Enable LaTeX rendering
plt.rc('text', usetex=True)
# Use plain TeX font (Computer Modern) for legend
plt.rc('legend', fontsize=fontsize_var, fancybox=True, framealpha=0.7, edgecolor='black', facecolor='white', borderaxespad=0.5, labelspacing=0.5)
plt.rc('font', family='serif', weight='normal')
# Plotting
plt.plot(minus_cp_2011[:, 0], minus_cp_2011[:, 1], label='3D DNS', color=profile_color[0],linewidth=linewidth_list[0])
plt.plot(minus_cp_2011[:, 0], minus_cp_2011[:, 2], label='', color=profile_color[0],linewidth=linewidth_list[0])
plt.plot(minus_cp[:, 0], minus_cp[:, 1], label='3D DNS-SLR', color=profile_color[1], linestyle='--',linewidth=linewidth_list[0])
plt.plot(minus_cp[:, 0], minus_cp[:, 2], label='', color=profile_color[1], linestyle='--',linewidth=linewidth_list[0])
plt.plot(minus_cp_2[:, 0], minus_cp_2[:, 1], label='3D DNS-SLRT', color=profile_color[2],linestyle='-.',linewidth=linewidth_list[0])
plt.plot(minus_cp_2[:, 0], minus_cp_2[:, 2], label='', color=profile_color[2],linestyle='-.',linewidth=linewidth_list[0])
plt.plot(minus_cp_2d_ff[:, 0], minus_cp_2d_ff[:, 1], label='2D DNS free field', color=profile_color[3],linewidth=linewidth_list[3], linestyle='-')
plt.plot(minus_cp_2d_ff[:, 0], minus_cp_2d_ff[:, 2], label='', color=profile_color[3],linewidth=linewidth_list[3], linestyle='-')
plt.plot(minus_cp_2d_lips[:, 0], minus_cp_2d_lips[:, 1], label='2D DNS lips', color=profile_color[4],linewidth=linewidth_list[4], linestyle=':')
plt.plot(minus_cp_2d_lips[:, 0], minus_cp_2d_lips[:, 2], label='', color=profile_color[4],linewidth=linewidth_list[4], linestyle=':')
plt.plot(minus_cp_expt[:, 0], minus_cp_expt[:, 1], '-s', label='Experiment (ECL)', linestyle='none', markersize=marker_size, markeredgecolor='blue', markerfacecolor='none')
plt.plot(minus_cp_expt2[:, 0], minus_cp_expt2[:, 1], '-^', label='Experiment (UdS)', linestyle='none', markersize=marker_size, markeredgecolor='black', markerfacecolor='none')

# Legend and labels
plt.legend(prop = { "size": fontsize_var*0.75 })
plt.xlabel('$x/c$', fontsize=fontsize_var, fontweight='bold')
plt.ylabel('$-\overline{C_p}$', fontsize=fontsize_var, fontweight='bold')
plt.xlim([-1, 0])
plt.ylim([-1.5, 3.0])
plt.grid(alpha=0.2, color='black')

plt.tick_params(axis='both', which='major', labelsize=fontsize_var)
ax = plt.gca()
ax.yaxis.offsetText.set(size=fontsize_var)
plt.yticks(fontsize=fontsize_var)
plt.xticks(fontsize=fontsize_var)

#plt.savefig(project_folder + 'mean_cp.jpg', format='jpeg', dpi=300)
plt.savefig(project_folder + 'mean_cp.pdf', dpi=300)

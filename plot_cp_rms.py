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
density = eval(settings.at["density", settings.columns[0]])
#signal analysis
chord = eval(settings.at["chord", settings.columns[0]])
U_ref_list = eval(settings.at["U_ref_list", settings.columns[0]])

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
minus_cp_data = np.genfromtxt(project_folder + 'mean_cp/' + 'prms_DNS_SLR_filtered.csv', delimiter=',')
minus_cp_data_2 = np.genfromtxt(project_folder + 'mean_cp/'+'prms_DNS_SLRT_filtered.csv', delimiter=',')
minus_cp_2011_data = np.genfromtxt(project_folder + 'mean_cp/'+'prms_DNS_filtered.csv', delimiter=',')
minus_cp_2d_ff_data = np.genfromtxt(project_folder + 'mean_cp/'+'prms_2d_ff_filtered.csv', delimiter=',')
minus_cp_2d_lips_data = np.genfromtxt(project_folder + 'mean_cp/'+'prms_2d_lips_filtered.csv', delimiter=',')
minus_cp_expt = np.genfromtxt(project_folder + 'mean_cp/'+'Cp_8deg_ECL.txt', delimiter='\t\t')
minus_cp_expt2 = np.genfromtxt(project_folder + 'mean_cp/'+'Cp_8deg_UdS.txt', delimiter='\t\t\t\t\t\t')

#############################################################################################
#                                        Calculate profile                      #
#############################################################################################

# Extract data columns
minus_cp = minus_cp_data[:, 0:]
minus_cp_2 = minus_cp_data_2[:, 0:]
minus_cp_2011 = minus_cp_2011_data[:, 0:]
minus_cp_2d_lips = minus_cp_2d_lips_data[:, 0:]
minus_cp_2d_ff = minus_cp_2d_ff_data[:, 0:]

# Coordinate adjustment
minus_cp[:, 0] -= minus_cp[1, 0]
minus_cp_2[:, 0] -= minus_cp_2[1, 0]
minus_cp_2011[:, 0] -= minus_cp_2011[1, 0]
minus_cp_2d_lips[:, 0] -= minus_cp_2d_lips[1, 0]
minus_cp_2d_ff[:, 0] -= minus_cp_2d_ff[1, 0]

# Non-dimensionalize
minus_cp[:, 0] /= 0.1356
minus_cp_2[:, 0] /= 0.1356
minus_cp_2011[:, 0] /= 0.1356
minus_cp_2d_lips[:, 0] /= 0.1356
minus_cp_2d_ff[:, 0] /= 0.1356

#Normalize to get the Cprms
minus_cp[:, 1:] = minus_cp[:, 1:]/(0.5*density*U_ref_list[1]**2) #DNS-SLR
minus_cp_2[:, 1:] = minus_cp_2[:, 1:]/(0.5*density*U_ref_list[2]**2) #DNS-SLRT
minus_cp_2011[:, 1:] = minus_cp_2011[:, 1:]/(0.5*density*U_ref_list[0]**2) #DNS
minus_cp_2d_lips[:, 1:] = minus_cp_2d_lips[:, 1:]/(0.5*density*U_ref_list[3]**2) #2D lips
minus_cp_2d_ff[:, 1:] = minus_cp_2d_ff[:, 1:]/(0.5*density*U_ref_list[4]**2) #2D freefield

plt.figure(figsize=plot_size_rectangle)
# Enable LaTeX rendering
plt.rc('text', usetex=True)
# Use plain TeX font (Computer Modern) for legend
plt.rc('legend', fontsize=fontsize_var, fancybox=True, framealpha=0.7, edgecolor='black', facecolor='white', borderaxespad=0.5, labelspacing=0.5)
plt.rc('font', family='serif', weight='normal')
# Plotting
plt.plot(minus_cp_2011[:, 0], minus_cp_2011[:, 1], label='3D DNS', color=profile_color[0],linewidth=linewidth_list[0])
#plt.plot(minus_cp_2011[:, 0], minus_cp_2011[:, 2], label='', color=profile_color[0],linewidth=linewidth_list[0])
plt.plot(minus_cp[:, 0], minus_cp[:, 1], label='3D DNS-SLR', color=profile_color[1], linestyle='--',linewidth=linewidth_list[0])
#plt.plot(minus_cp[:, 0], minus_cp[:, 2], label='', color=profile_color[1], linestyle='--',linewidth=linewidth_list[0])
plt.plot(minus_cp_2[:, 0], minus_cp_2[:, 1], label='3D DNS-SLRT', color=profile_color[2],linestyle='-.',linewidth=linewidth_list[0])
#plt.plot(minus_cp_2[:, 0], minus_cp_2[:, 2], label='', color=profile_color[2],linestyle='-.',linewidth=linewidth_list[0])
plt.plot(minus_cp_2d_ff[:, 0], minus_cp_2d_ff[:, 1], label='2D DNS free field', color=profile_color[3],linewidth=linewidth_list[3], linestyle='-')
#plt.plot(minus_cp_2d_ff[:, 0], minus_cp_2d_ff[:, 2], label='', color=profile_color[3],linewidth=linewidth_list[3], linestyle='-')
plt.plot(minus_cp_2d_lips[:, 0], minus_cp_2d_lips[:, 1], label='2D DNS lips', color=profile_color[4],linewidth=linewidth_list[4], linestyle=':')
#plt.plot(minus_cp_2d_lips[:, 0], minus_cp_2d_lips[:, 2], label='', color=profile_color[4],linewidth=linewidth_list[4], linestyle=':')

#Add probe 3 and probe 5 location
p3_x = -0.128516/0.1356+1
p5_x = -0.12376/0.1356+1

plt.axvline(x=p3_x, color='black', linestyle=':',linewidth=1)
plt.axvline(x=p5_x, color='black', linestyle=':',linewidth=1)

# Legend and labels
plt.legend(prop = { "size": fontsize_var*0.75 })
plt.xlabel('$x/c$', fontsize=fontsize_var, fontweight='bold')
plt.ylabel('${C_p}_{rms}$', fontsize=fontsize_var, fontweight='bold')
plt.xlim([0, 1.0])
plt.ylim([0, 0.025])
plt.grid(alpha=0.2, color='black')

plt.tick_params(axis='both', which='major', labelsize=fontsize_var)
ax = plt.gca()
ax.yaxis.offsetText.set(size=fontsize_var)
plt.yticks(fontsize=fontsize_var)
plt.xticks(fontsize=fontsize_var)

#plt.savefig(project_folder + 'mean_cp.jpg', format='jpeg', dpi=300)
plt.tight_layout()
plt.savefig(project_folder + 'mean_cprms_filtered.pdf', dpi=300)

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
linestyle_list_b = eval(settings.at["linestyle_list_b", settings.columns[0]])
variable = 'mean_cf'
# Physical parameters:
M_exp = eval(settings.at["M_exp", settings.columns[0]])
M_PF = eval(settings.at["M_PF", settings.columns[0]])
U_0 = eval(settings.at["U_0", settings.columns[0]])
#signal analysis
chord = eval(settings.at["chord", settings.columns[0]])
U_ref_list = eval(settings.at["U_ref_list", settings.columns[0]])
U_inf_list = eval(settings.at["U_inf_list", settings.columns[0]])

#plot setting
linewidth_3D_result = eval(settings.at["linewidth_3D_result", settings.columns[0]])
linewidth_2D_result = eval(settings.at["linewidth_2D_result", settings.columns[0]])
linewidth_list = eval(settings.at["linewidth_list", settings.columns[0]])
marker_list = eval(settings.at["marker_list", settings.columns[0]])
marker_size = eval(settings.at["marker_size", settings.columns[0]])
fontsize_var = eval(settings.at["fontsize_var", settings.columns[0]])
plot_size_square = eval(settings.at["plot_size_square", settings.columns[0]])
plot_size_rectangle = eval(settings.at["plot_size_rectangle", settings.columns[0]])

U_ref_cfd = 16.62
U_ref_cfd_legacy = 16.97
U_ref_expt = 16.66

#############################################################################################
#                                        Calculate profile                      #
#############################################################################################

cfd_plot_data = np.genfromtxt('/home/ziyz1701/storage/CD_airfoil/AIAA_2023_paper_plot/outlet_profile/outlet_velocity_magnitude_trip_full.csv',skip_header=1,delimiter=',')#DNS SLRT 0.12c
cfd_plot_data_2 = np.genfromtxt('/home/ziyz1701/storage/CD_airfoil/AIAA_2023_paper_plot/outlet_profile/outlet_velocity_magnitude_refined.csv',skip_header=1,delimiter=',') #15 degree 0.36c
cfd_plot_data_legacy = np.genfromtxt('/home/ziyz1701/storage/CD_airfoil/AIAA_2023_paper_plot/outlet_profile/outlet_velocity_magnitude_legacy.csv',skip_header=1,delimiter=',')
expt_plot_data = np.loadtxt('/home/ziyz1701/storage/CD_airfoil/AIAA_2023_paper_plot/outlet_profile/outlet_velocity_magnitude_expt.csv',skiprows=1,delimiter=',')

cfd_plot_2d_ff = np.genfromtxt('/home/ziyz1701/storage/CD_airfoil/AIAA_2023_paper_plot/outlet_profile/outlet_velocity_magnitude_2d_freefield.csv',skip_header=1,delimiter=',')
cfd_plot_2d_lips = np.genfromtxt('/home/ziyz1701/storage/CD_airfoil/AIAA_2023_paper_plot/outlet_profile/outlet_velocity_magnitude_2d_lips.csv',skip_header=1,delimiter=',')

print(expt_plot_data)
# plot result
plt.figure(figsize=plot_size_rectangle)
# Enable LaTeX rendering
plt.rc('text', usetex=True)
# Use plain TeX font (Computer Modern) for legend
plt.rc('legend', fontsize=fontsize_var, fancybox=True, framealpha=0.7, edgecolor='black', facecolor='white', borderaxespad=0.5, labelspacing=0.5)
plt.rc('font', family='serif', weight='normal')

plt.plot(cfd_plot_data[:,0]/0.25,cfd_plot_data[:,1]/U_inf_list[2],color = profile_color[2],label = '3D DNS-SLRT',linestyle=linestyle_list_b[2],linewidth=linewidth_list[2])
plt.plot(cfd_plot_data_legacy[:,0]/0.25,cfd_plot_data_legacy[:,1]/U_inf_list[0],color = profile_color[1],label = '3D DNS-SLR',linestyle=linestyle_list_b[1],linewidth=linewidth_list[1])
plt.plot(cfd_plot_data_2[:,0]/0.25,cfd_plot_data_2[:,1]/U_inf_list[1],color = profile_color[0],label = '3D DNS',linestyle=linestyle_list_b[0],linewidth=linewidth_list[0])
plt.plot(expt_plot_data[:,0]/0.25,expt_plot_data[:,1]/np.max(expt_plot_data[:,1]),linestyle='none',marker='s', markersize=marker_size, mec='black', mfc='none',color = 'black', label = 'Experiment',linewidth=linewidth_list[0])
errors = np.ones(len(expt_plot_data[:, 0])) * 0.01
plt.errorbar(expt_plot_data[:, 0] / 0.25, 
             expt_plot_data[:, 1] / np.max(expt_plot_data[:, 1]),
             yerr=errors,  # Specify the error bars
             linestyle='none',
             marker='s',
             markersize=marker_size,
             barsabove=True,
             mec='black',
             mfc='none',
             color='black',
             label='Experiment (ECL)',
             linewidth=linewidth_list[0])

plt.plot((cfd_plot_2d_ff[:,0]-0.02)/0.25,cfd_plot_2d_ff[:,1]/U_inf_list[4],color = profile_color[3],label = '2D free field', linestyle=linestyle_list_b[3],linewidth=linewidth_list[3])
plt.plot(cfd_plot_2d_lips[:,0]/0.25,cfd_plot_2d_lips[:,1]/U_inf_list[3],color = profile_color[4],label = '2D lips', linestyle=linestyle_list_b[4],linewidth=linewidth_list[4])

plt.xlim((-0.1, 2.1))
plt.ylim((0.8,1.05))
plt.xlabel(r'$Y/L_w$', fontsize=fontsize_var, fontweight='bold')
plt.ylabel(r'$U/U_{ref}$', fontsize=fontsize_var, fontweight='bold')
plt.grid(alpha=0.2, color='black')

plt.tick_params(axis='both', which='major', labelsize=fontsize_var)
ax = plt.gca()
ax.yaxis.offsetText.set(size=fontsize_var)
plt.yticks(fontsize=fontsize_var)
plt.xticks(fontsize=fontsize_var)

plt.legend(prop = { "size": fontsize_var })
plt.savefig(project_folder + 'outlet_profile.pdf', dpi=300)
plt.close()

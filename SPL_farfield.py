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

# Load FWH results
spectrum1 = np.loadtxt('/home/ziyz1701/storage/CD_airfoil/AIAA_2023_paper_plot/farfield_SPL/SPL_f_PSD.dat')
spectrum2 = np.loadtxt('/home/ziyz1701/storage/CD_airfoil/AIAA_2023_paper_plot/farfield_SPL/DNS_SLRT/sherfwh/psd_sherfwh_2m_90deg_current_0p05s_start.txt')
spectrum3 = np.loadtxt('/home/ziyz1701/storage/CD_airfoil/AIAA_2023_paper_plot/farfield_SPL/DNS/sherfwh/psd_sherfwh_2m_90deg_legacy_0p05s_start.txt')
spectrum4 = np.loadtxt('/home/ziyz1701/storage/CD_airfoil/AIAA_2023_paper_plot/farfield_SPL/DNS_SLR/sherfwh/psd_sherfwh_2m_90deg_shear_layer_refined_0p05s_start.txt')
spectrum5 = np.loadtxt('/home/ziyz1701/storage/CD_airfoil/AIAA_2023_paper_plot/farfield_SPL/expe_points.txt')

#Load Direct results
spectrum6 = np.loadtxt('/home/ziyz1701/storage/CD_airfoil/AIAA_2023_paper_plot/farfield_SPL/DNS_SLRT/direct/psd_pressure_r1_06.txt')
spectrum7 = np.loadtxt('/home/ziyz1701/storage/CD_airfoil/AIAA_2023_paper_plot/farfield_SPL/DNS/direct/psd_pressure_r1_06_legacy.txt')
spectrum8 = np.loadtxt('/home/ziyz1701/storage/CD_airfoil/AIAA_2023_paper_plot/farfield_SPL/DNS_SLR/direct/psd_pressure_r1_06_shear_layer_refined.txt')

#############################################################################################
#                                        Calculate profile                      #
#############################################################################################

plt.figure(figsize=plot_size_rectangle)
# Enable LaTeX rendering
plt.rc('text', usetex=True)
# Use plain TeX font (Computer Modern) for legend
plt.rc('legend', fontsize=fontsize_var, fancybox=True, framealpha=0.7, edgecolor='black', facecolor='white', borderaxespad=0.5, labelspacing=0.5)
plt.rc('font', family='serif', weight='normal')

plt.tick_params(axis='both', which='major', labelsize=fontsize_var)
ax = plt.gca()
ax.yaxis.offsetText.set(size=fontsize_var*1.3)
plt.yticks(fontsize=fontsize_var*1.3)
plt.xticks(fontsize=fontsize_var*1.3)
plt.grid(alpha=0.2, color='black')

plt.semilogx(spectrum1[:, 0], spectrum1[:, 1] + 20 * np.log10(1.21 / 2.00), color='black',label='UdeS experiment')
plt.semilogx(spectrum2[:, 0], spectrum2[:, 1], color=profile_color[2],label=name_list[2],linestyle=linestyle_list_b[2],linewidth=linewidth_list[0])
plt.semilogx(spectrum3[:, 0], spectrum3[:, 1], color=profile_color[0],label=name_list[0],linestyle=linestyle_list_b[0],linewidth=linewidth_list[0])
plt.semilogx(spectrum4[:, 0], spectrum4[:, 1], color=profile_color[1],label=name_list[1],linestyle=linestyle_list_b[1],linewidth=linewidth_list[0])
plt.semilogx(spectrum5[:, 0], spectrum5[:, 1], 'ks', markersize=marker_size, markerfacecolor='none', label='ECL experiment')

plt.legend(prop = { "size": fontsize_var })
plt.xlim([100, 10000])
plt.ylim([-30, 40])
plt.xlabel('Frequency (Hz)', fontsize=fontsize_var*1.3, fontweight='bold')
plt.ylabel('Power Spectral Density (dB/Hz)', fontsize=fontsize_var*1.3, fontweight='bold')
#plt.savefig(project_folder + 'fwh_noise.jpg', format='jpeg', dpi=300)
plt.tight_layout()
plt.savefig(project_folder + 'fwh_noise.pdf', dpi=300)
plt.close()
######################################################################################################################################################

plt.figure(figsize=plot_size_rectangle)
# Enable LaTeX rendering
plt.rc('text', usetex=True)
# Use plain TeX font (Computer Modern) for legend
plt.rc('legend', fontsize=fontsize_var, fancybox=True, framealpha=0.7, edgecolor='black', facecolor='white', borderaxespad=0.5, labelspacing=0.5)
plt.rc('font', family='serif', weight='normal')

plt.tick_params(axis='both', which='major', labelsize=fontsize_var)
ax = plt.gca()
ax.yaxis.offsetText.set(size=fontsize_var*1.3)
plt.yticks(fontsize=fontsize_var*1.3)
plt.xticks(fontsize=fontsize_var*1.3)
plt.grid(alpha=0.5, color='black')


# Compare Direct results

plt.semilogx(spectrum1[:, 0], spectrum1[:, 1] + 20 * np.log10(1.21 / 2.00), color='black',label='UdeS experiment')
plt.semilogx(spectrum6[:, 0], spectrum6[:, 1], color=profile_color[2],label=name_list[2],linestyle=linestyle_list_b[2],linewidth=linewidth_list[2])
plt.semilogx(spectrum7[:, 0], spectrum7[:, 1], color=profile_color[0],label=name_list[0],linestyle=linestyle_list_b[0],linewidth=linewidth_list[0])
plt.semilogx(spectrum8[:, 0], spectrum8[:, 1], color=profile_color[1],label=name_list[1],linestyle=linestyle_list_b[1],linewidth=linewidth_list[1])
plt.semilogx(spectrum5[:, 0], spectrum5[:, 1], 'ks', markersize=marker_size, markerfacecolor='none', label='ECL experiment')

plt.xlim([100, 10000])
plt.ylim([-30, 40])
plt.legend(prop = { "size": fontsize_var})
plt.xlabel('Frequency (Hz)', fontsize=fontsize_var*1.3, fontweight='bold')
plt.ylabel('Power Spectral Density dB/Hz', fontsize=fontsize_var*1.3, fontweight='bold')
plt.tight_layout()
plt.savefig(project_folder + 'direct_noise.pdf', dpi=300)
plt.close()

import numpy as np
#This scripts compare velocity magnitude from the simulation to velocity measurements using a single-normal probe from Doug Neal's thesis.

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

#############################################################################################
#                                        Calculate profile                      #
#############################################################################################

for probe in probe_number_list:

	probe_number = probe
	# Import data
	LegacyBL = np.genfromtxt(project_folder + 'boundary_layer_profile/' + "Bl-data-2011/probe-{}-bl-2011-extracted.csv".format(probe_number), delimiter=',')
	BL_data_refined = np.genfromtxt(project_folder + 'boundary_layer_profile/' + "BL_shear_layer_refined/probe-{}-bl-shear-layer-refined.csv".format(probe_number), delimiter=',')
	BL_data_tripped = np.genfromtxt(project_folder + 'boundary_layer_profile/' + "BL_shear_layer_refined_tripped/probe_{}_bl.csv".format(probe_number), delimiter=',')
	BL_data_2d_ff = np.genfromtxt(project_folder + 'boundary_layer_profile/' + "BL_2D_freefield/probe_{}_bl.csv".format(probe_number), delimiter=',')
	BL_data_2d_lips = np.genfromtxt(project_folder + 'boundary_layer_profile/' + "BL_2D_lips/probe_{}_bl.csv".format(probe_number), delimiter=',')
	BL_damp_corrected = np.genfromtxt(project_folder + 'boundary_layer_profile/' + "BL_DNS_SLRT_damping_corrected/probe_{}_vel.csv".format(probe_number), delimiter=',')
	ExptData = np.genfromtxt(project_folder + 'boundary_layer_profile/' + "BL-exp-data-digitized/probe-{}-bl.csv".format(probe_number), delimiter=',')

	# Non-dimensionalize boundary layer
	chord = 0.1356
	U_ref = 16.7112
	U_ref_2011 = 16.7  # recomputed (16.6114)

	BL_data_tripped[:, 0] /= chord
	BL_data_tripped[:, 1] /= U_ref_list[2]
	LegacyBL[:, 0] /= chord
	LegacyBL[:, 1] /= U_ref_list[0]
	BL_data_refined[:, 0] /= chord
	BL_data_refined[:, 1] /= U_ref_list[1]
	
	BL_data_2d_ff[:, 0] /= chord
	BL_data_2d_ff[:, 1] /= U_ref_list[4]
	BL_data_2d_lips[:, 0] /= chord
	BL_data_2d_lips[:, 1] /= U_ref_list[3]

	plt.figure(figsize=plot_size_rectangle)
	# Enable LaTeX rendering
	plt.rc('text', usetex=True)
	# Use plain TeX font (Computer Modern) for legend
	plt.rc('legend', fontsize=fontsize_var, fancybox=True, framealpha=0.7, edgecolor='black', facecolor='white', borderaxespad=0.5, labelspacing=0.5)
	plt.rc('font', family='serif', weight='normal')
	# Plot graph
	plt.plot(LegacyBL[:, 1], LegacyBL[:, 0], label=name_list[0],color = profile_color[0],linestyle=linestyle_list_b[0],linewidth=linewidth_list[0])
	plt.plot(BL_data_tripped[:, 1], BL_data_tripped[:, 0], label=name_list[2],color = profile_color[2],linestyle=linestyle_list_b[2],linewidth=linewidth_list[2])
	plt.plot(BL_data_refined[:, 1], BL_data_refined[:, 0], label=name_list[1],color = profile_color[1],linestyle=linestyle_list_b[1],linewidth=linewidth_list[1])
	plt.plot(BL_data_2d_ff[:, 1], BL_data_2d_ff[:, 0], label=name_list[3],color = profile_color[3],linestyle=linestyle_list_b[3],linewidth=linewidth_list[3])
	plt.plot(BL_data_2d_lips[:, 1], BL_data_2d_lips[:, 0], label=name_list[4],color = profile_color[4],linestyle=linestyle_list_b[4],linewidth=linewidth_list[4])
	plt.plot(ExptData[:-1, 0], ExptData[:-1, 1], '-s', linestyle='none', markersize=marker_size, markeredgecolor='black', markerfacecolor='none',label='Experiment')
	
	#plt.plot(BL_damp_corrected[:, 1], BL_damp_corrected[:, 0], label=name_list[2],color = profile_color[2],linestyle=linestyle_list_b[2],linewidth=linewidth_list[2])
	plt.xlim([0, 2.2])
	plt.ylim([0, 0.08])
	plt.legend(prop = { "size": fontsize_var })
	plt.grid(alpha=0.2, color='black')

	plt.xlabel(r'$U/U_{ref}$', fontsize=fontsize_var*1.2, fontweight='bold')
	plt.ylabel(r'$y/c$', fontsize=fontsize_var*1.2, fontweight='bold')

	plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText=True)
	ax = plt.gca()
	ax.yaxis.offsetText.set(size=fontsize_var*1.35)
	plt.yticks(fontsize=fontsize_var*1.35)
	plt.xticks(fontsize=fontsize_var*1.35)

	# Save the figure
	#plt.savefig(project_folder + 'probe_plot_{}.jpg'.format(probe_number), format='jpeg', dpi=300)
	plt.savefig(project_folder + 'probe_plot_{}.pdf'.format(probe_number), dpi=300)

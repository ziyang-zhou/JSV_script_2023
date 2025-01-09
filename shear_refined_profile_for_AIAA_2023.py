import numpy as np
#This scripts compare velocity magnitude in the shear layer from the simulation to expt

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

U_ref_list=eval(settings.at["U_ref_list", settings.columns[0]])
#############################################################################################
#                                        Calculate profile                      #
#############################################################################################

dataset = ["002mm","010mm","030mm","110mm","260mm","340mm"]
U_ref_2d_lip = U_ref_list[3]
U_ref_2d_ff = U_ref_list[4]
U_ref_cfd_legacy = U_ref_list[0]
U_ref_cfd_refined = U_ref_list[1]
U_ref_cfd_trip = U_ref_list[2]
U_ref_expt = 16.66
U_ref_area_weighted_avg = 16.4
L_w = 25 #Jet half width of the ECL wind tunnel

for distance in dataset:
  cfd_plot_data = np.genfromtxt(project_folder + 'shear_layer_profile/cfd_refined_tripped/shear_layer_{}.csv'.format(distance),skip_header=1,delimiter=',')
  cfd_plot_data_1=np.genfromtxt(project_folder +'shear_layer_profile/cfd_reference_dns/shear_layer_{}.csv'.format(distance),skip_header=1,delimiter=',')
  cfd_plot_data_2=np.genfromtxt(project_folder +'shear_layer_profile/cfd_shear_layer_refined/shear_layer_{}.csv'.format(distance),skip_header=1,delimiter=',')
  cfd_plot_data_3=np.genfromtxt(project_folder +'shear_layer_profile/shear_layer_2d_lips/shear_layer_{}.csv'.format(distance),skip_header=1,delimiter=',')
  expt_plot_data = np.loadtxt('/home/ziyz1701/storage/CD_airfoil/3D_CD_lip_simulation/no_slip_wall_0p12c_trip_shifted_full_ver2/shear_layer_profile/expt/shear_layer_{}'.format(distance))
  
  plt.figure(figsize=plot_size_rectangle)
  # Enable LaTeX rendering
  plt.rc('text', usetex=True)
  # Use plain TeX font (Computer Modern) for legend
  plt.rc('legend', fontsize=fontsize_var, fancybox=True, framealpha=0.7, edgecolor='black', facecolor='white', borderaxespad=0.5, labelspacing=0.5)
  plt.rc('font', family='serif', weight='normal')

  plt.tick_params(axis='both', which='major', labelsize=fontsize_var)
  ax = plt.gca()
  ax.yaxis.offsetText.set(size=fontsize_var*1.35)
  plt.yticks(fontsize=fontsize_var*1.35)
  plt.xticks(fontsize=fontsize_var*1.35)
  plt.grid(alpha=0.2, color='black')
  
  # plot result
  plt.plot(1000*cfd_plot_data[:,0],cfd_plot_data[:,1]/U_ref_cfd_trip, color = profile_color[2],label=name_list[2],linestyle=linestyle_list_b[2],linewidth=linewidth_list[2])
  plt.plot(1000*cfd_plot_data_1[:,0],cfd_plot_data_1[:,1]/U_ref_cfd_legacy, color = profile_color[0],label=name_list[0],linestyle=linestyle_list_b[0],linewidth=linewidth_list[0])
  plt.plot(1000*cfd_plot_data_2[:,0],cfd_plot_data_2[:,1]/U_ref_cfd_refined, color = profile_color[1],label=name_list[1],linestyle=linestyle_list_b[1],linewidth=linewidth_list[1])
  plt.plot(1000*cfd_plot_data_3[:,0],cfd_plot_data_3[:,1]/U_ref_2d_lip, color = profile_color[4],label=name_list[4],linestyle=linestyle_list_b[4],linewidth=linewidth_list[4])
  plt.plot(expt_plot_data[:,1],expt_plot_data[:,6], '-s',color = 'black', label = 'Experiment', linestyle='', markersize=marker_size, markeredgecolor='black', markerfacecolor='none')
  plt.legend(prop = { "size": fontsize_var})
  #plt.xlim((100, 10000))
  plt.ylim((0,1))
  plt.xlabel(r'$Y [mm]$', fontsize=fontsize_var, fontweight='bold')
  plt.ylabel(r'$U/U_{ref}$', fontsize=fontsize_var, fontweight='bold')
  #plt.savefig(project_folder+'shear_layer_profile_{}.jpg'.format(distance), format='jpeg', dpi=300)
  plt.savefig(project_folder+'shear_layer_profile_{}.pdf'.format(distance), dpi=300)
  plt.close()

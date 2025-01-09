# This script plots U+ against y+. 

import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import math
from scipy.ndimage import gaussian_filter
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter
from scipy.interpolate import CubicSpline
import pandas as pd
import pdb
# Import the ScalarFormatter class
from matplotlib.ticker import ScalarFormatter
sys.path.append('/home/ziyz1701/storage/CD_airfoil/tool/analysis_folder')  # Add the parent directory to the path
from functions import analysis, futils

#############################################################################################
#                                        Load Inputs                                        #
#############################################################################################
# Load the settings dataframe:
settings = pd.read_csv("/home/ziyz1701/storage/CD_airfoil/TSFP_2024/plot_script/setting.csv", index_col= 0)
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

def probe_number_converter(probe):
	# Convert RMP probe number to x coord number
	if probe == '5':
			probe_x_position = -0.12376
			probe_y_position = 0.0212125
	elif probe == '9':
			probe_x_position = -0.0632463
			probe_y_position = 0.0199974
	elif probe == '21':
			probe_x_position = -0.019227
			probe_y_position = 0.00794334
	elif probe == '3':
			probe_x_position = -0.128516
			probe_y_position = 0.0207075
	elif probe == '25':
			probe_x_position = -0.003036
			probe_y_position = 0.00191532
	elif probe == '4':
			probe_x_position = -0.126493
			probe_y_position = 0.0180604
	elif probe == '8':
			probe_x_position = -0.0814612
			probe_y_position = 0.016288
	elif probe == '10':
			probe_x_position = -0.0637522
			probe_y_position = 0.0138381
	elif probe == '29':
			probe_x_position = -0.009613
			probe_y_position = 0.00185464
	elif probe == '24':
			probe_x_position = -0.010625
			probe_y_position = 0.0047688
	elif probe == '1':
			probe_x_position = -0.13388
			probe_y_position = 0.0194865
	elif probe == '2':
			probe_x_position = -0.131552
			probe_y_position = 0.0203201
	elif probe == '7':
			probe_x_position = -0.0809552
			probe_y_position = 0.0221024
	return probe_x_position,probe_y_position

casename_list = ['DNS','DNS-SLR','DNS-SLRT']

#plot setting
h_scale_factor=0.95
fontsize_var = fontsize_var*h_scale_factor
plt.figure(figsize=plot_size_rectangle)
# Enable LaTeX rendering
plt.rc('text', usetex=True)
# Use plain TeX font (Computer Modern) for legend
#plt.rc('legend', fontsize=fontsize_var*h_scale_factor, fancybox=True, framealpha=0.7, edgecolor='black', facecolor='white', borderaxespad=0.5, labelspacing=0.5)
plt.rc('font', family='serif', weight='normal')
fig,ax = plt.subplots(1,4,figsize=(15,10),sharey=True,sharex=False)
iplot = 0

probe_number_list = ['7','9','21','24']
data_folder='../bl_stat/'
density = 1.251
viscosity = 1.44e-5

for probe_number in probe_number_list:
	for casename in casename_list:
		probe_x_position,probe_y_position = probe_number_converter(probe_number)
		coordinate_df = pd.read_csv(data_folder + '{}/L08_coordinates.csv'.format(casename))
		idxBL = np.min([np.argmin(abs(np.array(coordinate_df['xBL'])-(probe_x_position/chord + 1.0))),len(np.array(coordinate_df['xBL']))-2])
		print('idxBL',idxBL,'#############################################################################')
		bl_data = pd.read_csv(data_folder + '{}/L08_BL_{}.csv'.format(casename,str(idxBL).zfill(3)))
		surface_data = pd.read_csv(data_folder + '{}/L08_surface_parameter.csv'.format(casename))
		tau_w = surface_data['tau_wall'][idxBL]
		u_tau = np.sqrt(tau_w/density)
		y_plus = bl_data['h']*u_tau/viscosity
		u_plus = bl_data['U_t']/u_tau
		ax[iplot].plot(y_plus,u_plus, label=casename)
		ax[iplot].set_xscale('log')
	iplot += 1
	print('iplot is {}',iplot)

ax[iplot-1].legend(prop = { "size": fontsize_var})
plt.savefig(data_folder + 'bl_plot_{}'.format(probe_number))

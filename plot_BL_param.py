# This script plots all relevant BL parameter for PSD prediction. 
# These include delta_95

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

casename_list = ['DNS','DNS-SLR','DNS-SLRT']
var_list = ['wall_pressure','pressure_scale','y_plus','delta_95','delta_theta','delta_star','tau_wall','dpds','Ue','beta_c','H']
var_list_latex = ['\overline{p}','\\tau_w^2\delta^*/U_e','y+','\delta_{95} [m]','\\theta [m]','\delta^* [m]','\tau_w [Pa]','dP/ds [Pa/m]','U_e [m/s]',r'\beta_c','H']

idx_start = 1 # Start where the bl was extracted well
gaussian_sigma = 0.3

for ivar,var in enumerate(var_list):
	print('plotting {}...'.format(var))
	plt.figure(figsize=plot_size_rectangle)
	# Enable LaTeX rendering
	plt.rc('text', usetex=True)
	# Use plain TeX font (Computer Modern) for legend
	plt.rc('legend', fontsize=fontsize_var, fancybox=True, framealpha=0.7, edgecolor='black', facecolor='white', borderaxespad=0.5, labelspacing=0.5)
	plt.rc('font', family='serif', weight='normal')

	for i,case in enumerate(casename_list):
		filename = '/home/ziyz1701/storage/CD_airfoil/AIAA_2023_paper_plot/bl_stat/{}/L08_surface_parameter.csv'.format(case)
		data = pd.read_csv(filename,index_col=None,engine='python')
		if var == 'pressure_scale':
			tau_wall_smoothed = gaussian_filter(data['tau_wall'].values[idx_start:-1],sigma=3.0)
			delta_star_smoothed = gaussian_filter(data['delta_star'].values[idx_start:-1],sigma=3.0)
			pressure_scale = tau_wall_smoothed**2*delta_star_smoothed
			#data['Ue'].values[idx_start:-1]
			plt.plot(data['streamwise_location'].values[idx_start:-1]/0.1317,pressure_scale,color = profile_color[i],linewidth=linewidth_list[0],linestyle = linestyle_list_b[i],label='3D '+case)
		elif var == 'y_plus':
			y_plus = np.sqrt(data['tau_wall'].values[idx_start:-1]/1.251)*1.48e-5/1.44e-5
			plt.plot(data['streamwise_location'].values[idx_start:-1]/0.1317,y_plus,color = profile_color[i],linewidth=linewidth_list[0],linestyle = linestyle_list_b[i],label='3D '+case)
		elif var == 'H':
			delta_star_smoothed = gaussian_filter(data['delta_star'].values[idx_start:-1],sigma=3.0)
			theta_smoothed = gaussian_filter(data['delta_theta'].values[idx_start:-1],sigma=3.0)
			H = delta_star_smoothed/theta_smoothed
			plt.plot(data['streamwise_location'].values[idx_start:-1]/0.1317,H,color = profile_color[i],linewidth=linewidth_list[0],linestyle = linestyle_list_b[i],label=case)
			
		elif var in ['wall_pressure', 'dpds', 'beta_c']:
			cp_data = pd.read_csv('/home/ziyz1701/storage/CD_airfoil/AIAA_2023_paper_plot/bl_stat/{}/mean_cp_{}.csv'.format(case,case))
			cp_data = cp_data[~cp_data['Position (Up)[Length:m]'].duplicated(keep='first')]
			cp_cs = CubicSpline(cp_data['Position (Up)[Length:m]'].values-cp_data['Position (Up)[Length:m]'].iloc[0],cp_data['Value (Up)[DynamicPressure/(Density*Velocity^2):Pa*m*sec^2/kg]'].values)
			wall_pressure = cp_cs(data['streamwise_location'].values)*(0.5*1.251*16.6**2)
			if var == 'wall_pressure':
				plt.plot(data['streamwise_location'].values/0.1317,wall_pressure,color = profile_color[i],linewidth=linewidth_list[0],linestyle = linestyle_list_b[i],label='3D '+case)
			elif var in ['dpds','beta_c']:
				print('Computing pressure gradient...')
				scoor = data['streamwise_location'].values
				#wall_pressure = gaussian_filter(wall_pressure,sigma=gaussian_sigma)
				wall_pressure = savgol_filter(wall_pressure,35,2)
				dpds = np.diff(wall_pressure)/np.diff(scoor) #Note that this is taking the pressure obtained using powerviz from surface data
				dpds_interp = np.interp(scoor,scoor[:-1],dpds)
				if var == 'beta_c':
					print('Computing beta_c')
					scoor = data['streamwise_location'].values
					#wall_pressure = gaussian_filter(wall_pressure,sigma=gaussian_sigma)
					wall_pressure = savgol_filter(wall_pressure,35,2)
					dpds = np.diff(wall_pressure)/np.diff(scoor) #Note that this is taking the pressure obtained using powerviz from surface data
					dpds_interp = np.interp(scoor,scoor[:-1],dpds)
					tau_wall_smoothed = gaussian_filter(data['tau_wall'].values[idx_start:-1],sigma=gaussian_sigma)
					delta_star_smoothed = gaussian_filter(data['delta_star'].values[idx_start:-1],sigma=gaussian_sigma)
					Re_tau = data['Ue'].values[idx_start:-1]*tau_wall_smoothed/1.44e-5
					plt.plot(data['streamwise_location'].values[idx_start:-1]/0.1317,delta_star_smoothed/tau_wall_smoothed*dpds_interp[idx_start:-1],color = profile_color[i],linewidth=linewidth_list[0],linestyle = linestyle_list_b[i],label=case)
					
				else:
					plt.plot(data['streamwise_location'].values[idx_start:-1]/0.1317,dpds_interp[idx_start:-1],color = profile_color[i],linewidth=linewidth_list[0],linestyle = linestyle_list_b[i],label=case)

		else:	
			plt.plot(data['streamwise_location'].values[idx_start:-1]/0.1317,data[var].values[idx_start:-1],color = profile_color[i],linewidth=linewidth_list[0],linestyle = linestyle_list_b[i],label=case)

	plt.legend(prop = { "size": fontsize_var })
	plt.grid(alpha=0.2, color='black')
	plt.xlabel(r'$x/c$', fontsize=fontsize_var*1.2, fontweight='bold')
	plt.ylabel(r'${}$'.format(var_list_latex[ivar]), fontsize=fontsize_var*1.2, fontweight='bold')
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText=True)
	ax = plt.gca()
	ax.yaxis.offsetText.set(size=fontsize_var*1.35)
	plt.yticks(fontsize=fontsize_var*1.35)
	plt.xticks(fontsize=fontsize_var*1.35)
	#plt.xlim((0.5,1.0))
	#plt.ylim((0.0015,0.0055))
	plt.tight_layout()
	plt.savefig('../bl_stat/'+'{}.pdf'.format(var), dpi=300)

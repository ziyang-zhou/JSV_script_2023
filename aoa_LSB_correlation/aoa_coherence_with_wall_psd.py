import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
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

pressure_or_suction = 'suction' #suction or pressure
#############################################################################################
#                                        Calculate profile                      #
#############################################################################################
#This script correlates the signal in mic_data to the variable found in the probes in probe_path

save_path = 'aoa_correlation/'
read_path = 'aoa_correlation/'

#signal 1 path
filename = 'upstream_probe_aoa/temporal_plane_jet_020_010.csv'
#signal 2 path
folder_path =  'cfd_refined_tripped/{}/'.format(pressure_or_suction)

nb_chunk = 1 #was 12
tchar = 0.1356/16
tstart = 0*tchar # initial t to start from

def next_greater_power_of_2(x):  
    return 2**(x-1).bit_length()

#Load signal 1 #########################################################################################################
sig_1_data = pd.read_csv(project_folder + read_path + filename)
sig_1_x_velocity = np.array(sig_1_data['x_velocity'])
sig_1_y_velocity = np.array(sig_1_data['y_velocity'])
aoa = analysis.get_aoa(sig_1_x_velocity,sig_1_y_velocity)
keep_time = sig_1_data['time'] > tstart
aoa = aoa[:] - aoa[keep_time].mean()
sig_1 = np.transpose([sig_1_data['time'],aoa])
#Load signal 2 #########################################################################################################
# Path to the folder containing CSV files
folder_path = project_folder + read_path + 'wall_psd/' + folder_path
# Function to extract the index from the last three digits of the filename
def extract_index(file_name):
    return int(file_name.split('_')[-1].split('.')[0])
# Lists to store 'time' and 'y_velocity' data
time_list = []
y_velocity_list = []
x_velocity_list = []
indexes_array = []
# Loop through files in the folder
for file_name in os.listdir(folder_path):
    if file_name.endswith('.csv'):
        # Extract index from the last three digits of the filename
        index = extract_index(file_name)
        # Read CSV file
        file_path = os.path.join(folder_path, file_name)
        data = pd.read_csv(file_path)
        # Append data to lists
        time_list.append(data['time'])
        y_velocity_list.append(data['static_pressure'])
        indexes_array.append(index)
# Combine lists into a NumPy array with the specified structure
#Perform signal processing #########################################################################################################
#Correlation
figure_scale_factor = 1.3
plt.figure(figsize=(plot_size_rectangle[0]//figure_scale_factor,plot_size_rectangle[1]//figure_scale_factor))
# Enable LaTeX rendering
plt.rc('text', usetex=True)
# Use plain TeX font (Computer Modern) for legend
plt.rc('legend', fontsize=fontsize_var*0.7, fancybox=True, framealpha=0.7, edgecolor='black', facecolor='white', borderaxespad=0.5, labelspacing=0.5)
plt.rc('font', family='serif', weight='normal')

for n,index in enumerate(indexes_array):
	
	if (time_list[n].iloc[-1] > sig_1[-1,0]) and (time_list[n][0] < sig_1[0,0]): #truncate signal 2 to fit timeframe of signal 1
		time_list_i = (time_list[n] <= sig_1[-1,0]) & (time_list[n] >= sig_1[0,0])
		time_list[n] = time_list[n][time_list_i]
		sig_2_y_velocity = np.array(y_velocity_list[n][time_list_i]) 
	else:
		sig_1_i = (sig_1[:,0] <= time_list[n].iloc[-1]) & (sig_1[:,0] >= time_list[n][0]) #truncate signal 1 to fit timeframe of signal 2
		sig_1 = sig_1[sig_1_i]
		sig_2_y_velocity = np.array(y_velocity_list[n]) 
	
	#if len(time_list[n]) != len(sig_2_y_velocity):
	#	pdb.set_trace()
	interpolated_sig_2 = np.interp(sig_1[:,0], time_list[n], sig_2_y_velocity) #interpolate signal 2 to accomodate fs of signal 1
	print('sig_1',len(sig_1))
	print('sig_2',len(sig_2_y_velocity))
	interpolated_sig_2 = interpolated_sig_2 - np.average(interpolated_sig_2)
	dt = sig_1[1,0] - sig_1[0,0]
	time_cross_corr,cross_norm = analysis.get_pearson_corr(sig_1[:,1],interpolated_sig_2,dt) #calculate the pearson correlation of sig_1 (aoa) and each o the y_velocity
	print('index {}'.format(index))
	line_width_list = ['6','1','3','0.5']
	through_flow_time = 0.1356/16.7114
	plt.plot(time_cross_corr/through_flow_time,cross_norm,label = 'RMP {}'.format(index),color = 'black',linewidth = line_width_list[n],linestyle=linestyle_list_b[n])

plt.tick_params(axis='both', which='major', labelsize=fontsize_var)
ax = plt.gca()
ax.yaxis.offsetText.set(size=fontsize_var)
plt.yticks(fontsize=fontsize_var)
plt.xticks(fontsize=fontsize_var)

plt.tight_layout()
plt.xlabel(r'$T_{lag} \times U_{ref}/c$', fontsize=fontsize_var, fontweight='bold')
plt.ylabel(r'$R(\alpha_i\acute{,}p_{wall})$', fontsize=fontsize_var, fontweight='bold')
plt.xlim([-70,70])
plt.ylim([-1.0,1.0])
plt.grid(alpha=0.2, color='black')

# Get the handles and labels of the current axes
handles, labels = plt.gca().get_legend_handles_labels()
labels_int = []
for x in labels:
	labels_int.append(x.split(' ')[1])
sorted_labels, sorted_handles = zip(*sorted(zip(map(int, labels_int), handles)))
sorted_labels = ['Probe {}'.format(x) for x in sorted_labels]
plt.legend(sorted_handles, sorted_labels,prop = { "size": fontsize_var*0.7 })
plt.tight_layout()
#plt.savefig(project_folder + save_path + 'aoa and wall psd {} side.jpg'.format(pressure_or_suction), format='jpeg', dpi=300)
plt.savefig(project_folder + save_path + 'aoa and wall psd {} side.pdf'.format(pressure_or_suction), dpi=300)
plt.close()

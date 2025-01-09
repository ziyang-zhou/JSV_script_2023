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
#This script correlates the signal in mic_data to the variable found in the probes in probe_path

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

save_path = 'aoa_correlation/'
read_path = 'aoa_correlation/'

#signal 1 path
filename = 'x_attach.txt'
#signal 2 path
folder_path = 'upstream_probe_aoa/'

nb_chunk = 1 #was 12
tchar = 0.1356/16
tstart = 0*tchar # initial t to start from

def next_greater_power_of_2(x):  
    return 2**(x-1).bit_length()

#Load signal 1 #########################################################################################################
sig_1 = np.loadtxt(project_folder + save_path + filename,skiprows=1)
keep_time = sig_1[:,0]>tstart
sig_1[:,1] = sig_1[:,1] - sig_1[keep_time,1].mean(axis=0)
#Load signal 2 #########################################################################################################
# Path to the folder containing CSV files
folder_path = project_folder + read_path + 'upstream_probe_aoa/'
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
        y_velocity_list.append(data['y_velocity'])
        
        x_velocity_list.append(data['x_velocity'])
        indexes_array.append(index)
# Combine lists into a NumPy array with the specified structure
#Perform signal processing #########################################################################################################
figure_scale_factor = 1.3
plt.figure(figsize=(plot_size_rectangle[0]//figure_scale_factor,plot_size_rectangle[1]//figure_scale_factor))
# Enable LaTeX rendering
plt.rc('text', usetex=True)
# Use plain TeX font (Computer Modern) for legend
plt.rc('legend', fontsize=fontsize_var, fancybox=True, framealpha=0.7, edgecolor='black', facecolor='white', borderaxespad=0.5, labelspacing=0.5)
plt.rc('font', family='serif', weight='normal')

for n,index in enumerate(indexes_array):
	if n == 0:
		time_list = np.array(time_list[0])
		time_list_i = (time_list <= sig_1[-1, 0]) & (time_list >= sig_1[0, 0]) #create the index for the time_list kept to within the confines of sig_1(aoa) time frame
		time_list = time_list[time_list_i]
		interpolated_sig_1 = np.interp(time_list, sig_1[:,0], sig_1[:,1]) #interpolate sig_1 (aoa) at the time points of the y_velocity array
	y_velocity = y_velocity_list[n][time_list_i] #construct the second signal which is to be correlated against
	x_velocity = x_velocity_list[n][time_list_i]
	sig_2 = analysis.get_aoa(x_velocity,y_velocity)#reconstruct the aoa
	dt = time_list[1] - time_list[0]
	print('dt : ',dt)
	time_cross_corr,cross_norm = analysis.get_pearson_corr(-interpolated_sig_1/0.1356,sig_2,dt) #calculate the pearson correlation of sig_1 (aoa) and each o the y_velocity
	
	print('index {}'.format(index))
	through_flow_time = 0.1356/16.7114
	plt.plot(time_cross_corr/through_flow_time,cross_norm,color = 'black',linestyle=linestyle_list_b[0],linewidth=linewidth_list[0])
	
sig_2 = np.array(sig_2) #convert the AoA data into a np array

plt.tick_params(axis='both', which='major', labelsize=fontsize_var)
ax = plt.gca()
ax.yaxis.offsetText.set(size=fontsize_var)
plt.yticks(fontsize=fontsize_var)
plt.xticks(fontsize=fontsize_var)

plt.xlabel(r'$T_{lag} \times U_{ref}/c$', fontsize=fontsize_var, fontweight='bold')
plt.ylabel(r'$R(x_{LSB},\alpha_i)$', fontsize=fontsize_var, fontweight='bold')
plt.xlim([-70,70])
plt.ylim([-1.0,1.0])
plt.grid(alpha=0.2, color='black')
#plt.savefig(project_folder + save_path + 'xcorr_aoa_x_attach.jpg', format='jpeg', dpi=300)
plt.tight_layout()
plt.savefig(project_folder + save_path + 'xcorr_aoa_x_attach.pdf', dpi=300)
plt.close()

import pandas as pd
import numpy as np
from scipy import signal as sig
import scipy.io
import matplotlib.pyplot as plt
import sys
import os
import re
import pdb
sys.path.append('/home/ziyz1701/storage/CD_airfoil/tool/analysis_folder')  # Add the parent directory to the path
from functions import analysis, futils
plt.rcParams['text.usetex'] = True
#############################################################################################
#                                        Load Inputs                                        #
#############################################################################################
# Load the settings dataframe:
settings = pd.read_csv("setting.csv", index_col= 0)
# Unpack Relevant Settings:
project_folder = settings.at["project_folder", settings.columns[0]]
plot_list = eval(settings.at["plot_list", settings.columns[0]])
expt_probe_list = eval(settings.at["expt_probe_list", settings.columns[0]])
data_folder_list = eval(settings.at["data_folder_list", settings.columns[0]])
#Signal processing parameter
tstart_cfd = eval(settings.at["tstart_cfd", settings.columns[0]])
tend_cfd = eval(settings.at["tend_cfd", settings.columns[0]])
expt_fs = eval(settings.at["expt_fs", settings.columns[0]])
nb_chunks = eval(settings.at["nb_chunks", settings.columns[0]])

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
name_list = eval(settings.at["wall_name_list", settings.columns[0]])
profile_color = eval(settings.at["profile_color", settings.columns[0]])

#############################################################################################
#                                        Load result                                        #
#############################################################################################
###UdeS experimental result
# Replace 'file1.npy' and 'file2.npy' with the actual file paths
expt_file1_path = project_folder + 'wall_psd/expt/processed/' + 'wallpressurespectra_8degree_16ms_LE.mat'
expt_file2_path = project_folder + 'wall_psd/expt/processed/' + 'wallpressurespectra_8degree_16ms_TE.mat'
# Load the .npy files
f_expt = scipy.io.loadmat(expt_file1_path)['frequency']
array1 = scipy.io.loadmat(expt_file1_path)['wallpressurespectra_8degree_16ms_LE']
array2 = scipy.io.loadmat(expt_file2_path)['wallpressurespectra_8degree_16ms_TE']
# Concatenate the trimmed arrays horizontally
#pressure history at wall probes arranged in order : from file 1 : 1,2,3,5,6,9 from file 2 : 21,23,24,25,26,27,28
expt_array = np.concatenate((array1, array2), axis=1)

###ECL experimental result
exp2_file_path = project_folder + 'wall_psd/expt/processed/' + 'wall_psd_ecl.csv'
expt2_array = pd.read_csv(exp2_file_path, index_col= 0)
f_expt2 = expt2_array.index.values

###CFD result
cfd_array_list = []
#Read CFD result
for data_folder in data_folder_list:
	folder_path = project_folder + 'wall_psd/' + data_folder
	# Get a list of all files in the folder
	files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]
	# Extract and convert the numeric part of the filenames
	numeric_parts = [int(re.search(r'\d+', f).group()) for f in files]
	# Sort files based on the extracted numeric parts
	sorted_indices = np.argsort(numeric_parts)
	sorted_files = [files[i] for i in sorted_indices]
	sorted_numeric_parts = [numeric_parts[i] for i in sorted_indices]
	# Initialize an empty list to store arrays from each file
	arrays = []
	# Load each CSV file, append the array to the list, and store the numeric part
	numeric_parts_array = np.array(sorted_numeric_parts)
	for file_name in sorted_files:
		file_path = os.path.join(folder_path, file_name)
		data = np.loadtxt(file_path, delimiter=',',skiprows=1)  # Assuming CSV files are comma-separated
		arrays.append(data[:,0:2])
	# Concatenate the arrays horizontally to create the final array
	cfd_array = np.column_stack(arrays)
	
	cfd_array_list.append(cfd_array)

#############################################################################################
#                                        Plot data                                          #
#############################################################################################
cfd_fs = 1/(cfd_array_list[0][1,0] - cfd_array_list[0][0,0]) #calculate sampling frequency of cfd measurement

#Correlation
plt.figure(figsize=plot_size_rectangle)
# Enable LaTeX rendering
plt.rc('text', usetex=True)
# Use plain TeX font (Computer Modern) for legend
plt.rc('legend', fontsize=fontsize_var, fancybox=True, framealpha=0.7, edgecolor='black', facecolor='white', borderaxespad=0.5, labelspacing=0.5)
plt.rc('font', family='serif', weight='normal')

plt.tick_params(axis='both', which='major', labelsize=fontsize_var)
ax = plt.gca()
ax.yaxis.offsetText.set(size=fontsize_var)
plt.yticks(fontsize=fontsize_var)
plt.xticks(fontsize=fontsize_var)

for i,n in enumerate(plot_list):
	#Correlation
	plt.figure(figsize=plot_size_rectangle)
	# Enable LaTeX rendering
	plt.rc('text', usetex=True)
	# Use plain TeX font (Computer Modern) for legend
	plt.rc('legend', fontsize=fontsize_var*0.75, fancybox=True, framealpha=0.7, edgecolor='black', facecolor='white', borderaxespad=0.5, labelspacing=0.5)
	plt.rc('font', family='serif', weight='normal')

	plt.tick_params(axis='both', which='major', labelsize=fontsize_var)
	ax = plt.gca()
	ax.yaxis.offsetText.set(size=fontsize_var)
	plt.yticks(fontsize=fontsize_var)	
	plt.xticks(fontsize=fontsize_var)
	
	expt_i = expt_probe_list.index(n) #UdeS experimental result
	expt2_i = expt2_array.columns.get_loc('Experimental sensor {}'.format(n)) #ECL experimental result
	plt.plot(f_expt,expt_array[:,expt_i],label = 'ECL',color = 'grey',linewidth = 2.0)
	plt.plot(f_expt2,expt2_array.iloc[:,expt2_i].values,label='UdeS',color='grey',linewidth = 6.0, linestyle = '--')
	
	for j,data_name in enumerate(name_list):
		print('plotting {}'.format(data_name))
		cfd_i = sorted_numeric_parts.index(n)
		if '2d' in data_name:
			tstart_cfd = 0.5
			tend_cfd = 0.6
		cfd_fs = 1/(cfd_array_list[j][1,0]-cfd_array_list[j][0,0])
		mic_cfd = cfd_array_list[j][:,2*cfd_i+1]#read the pressure array of current probe
		keep_time=((cfd_array_list[j][:,0]>tstart_cfd) & (cfd_array_list[j][:,0]<tend_cfd)) #truncate the initialization period of cfd data
		f_cfd,Pxx_cfd = analysis.get_pwelch(cfd_fs,mic_cfd-np.average(mic_cfd),nb_chunks,keep_time)
		plt.plot(f_cfd,10*np.log10(Pxx_cfd/4.0e-10),label = data_name, color = profile_color[j],linestyle=linestyle_list_b[j],linewidth=linewidth_list[j])
		print('frequency resolution',f_cfd[1]-f_cfd[0])
		print('t0',cfd_array_list[j][0,0],'t end',cfd_array_list[j][-1,0])

	plt.xlabel('Frequency [Hz]', fontsize=fontsize_var*1.35, fontweight='bold')
	plt.ylabel('SPL [dB/Hz]', fontsize=fontsize_var*1.35, fontweight='bold')
	plt.yticks(fontsize=fontsize_var*1.35)
	plt.xticks(fontsize=fontsize_var*1.35)
	ax=plt.gca()
	ax.set_xscale('log')
	plt.xlim((100, 10000))
	plt.ylim((20,100))
	# Add grid lines in log scale for both x and y axes
	ax.grid(True, which='both', linestyle='--', linewidth=0.5, axis='x')
	plt.legend(prop = { "size": fontsize_var*0.75 })
	
	# Add grid lines
	ax.grid(True, which='both', linestyle='--', linewidth=0.5)
	
	print('saving the figure...')
	plt.savefig(project_folder + 'wall_psd/wall_psd_probe_{}.pdf'.format(n))
	plt.savefig(project_folder + 'wall_psd/wall_psd_probe_{}.jpg'.format(n))
	plt.close()
	print('probe {} completed'.format(n))

	





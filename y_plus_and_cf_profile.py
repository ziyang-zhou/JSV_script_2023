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
variable = 'y_plus' #mean_cf or y_plus
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
plot_size_square = eval(settings.at["plot_size_square", settings.columns[0]])
plot_size_rectangle = eval(settings.at["plot_size_rectangle", settings.columns[0]])
#############################################################################################
#                                        Calculate profile                      #
#############################################################################################

fontsize_var = fontsize_var*1.2

plt.figure(figsize=plot_size_rectangle)
# Enable LaTeX rendering
plt.rc('text', usetex=True)
plt.rc('axes', linewidth=2)
# Use plain TeX font (Computer Modern) for legend
plt.rc('legend', fontsize=fontsize_var, fancybox=True, framealpha=0.7, edgecolor='black', facecolor='white', borderaxespad=0.5, labelspacing=0.5)
plt.rc('font', family='serif')
#plt.rcParams['text.latex.preamble'] = r'\usepackage{sfmath} \boldmath'

for n,folder in enumerate(folder_list):
    file_location = project_folder + 'mean_cf/' + folder
    data = np.loadtxt(file_location + variable, delimiter='\t')
    mask = (data[:, 0]/0.1356 >= -0.98) & (data[:, 0]/0.1356 <= -0.01)
    plt.plot(data[mask][:,0]/0.1356, data[mask][:,1], label=name_list[n],color=profile_color[n],linestyle=linestyle_list[n],linewidth=linewidth_list[n],marker=marker_list[n],markersize=marker_size,markerfacecolor='none')

# Add labels and legend
plt.ylabel(r'$\overline{y_+}$', fontsize=fontsize_var*0.8, weight='bold')
plt.xlabel(r'$x/c$', fontsize=fontsize_var, weight='bold')
if variable == 'mean_cf':
	plt.ylim([-0.005,0.015])
if variable == 'y_plus':
	plt.ylim([0.0,1.5])
# Add legend
plt.legend(prop = { "size": fontsize_var,"weight": 'bold' })
plt.grid(alpha=0.2, color='black')

if variable == 'y_plus':
	plt.tick_params(axis='both', which='major', labelsize=fontsize_var,width=2)
	ax = plt.gca()
	ax.yaxis.offsetText.set(size=fontsize_var)
	ax.xaxis.set_tick_params(labelsize=fontsize_var)
	ax.yaxis.set_tick_params(labelsize=fontsize_var)

if variable == 'mean_cf':
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText=True)
	ax = plt.gca()
	ax.yaxis.offsetText.set(size=fontsize_var)
	plt.yticks(fontsize=fontsize_var)
	plt.xticks(fontsize=fontsize_var)

futils.mkdir(save_folder)
plt.savefig(save_folder + variable + '.pdf')
print('save in',save_folder)





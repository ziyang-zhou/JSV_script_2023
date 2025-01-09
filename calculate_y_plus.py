import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy.interpolate import interp1d,CubicSpline
import pdb
import pandas as pd


# Load the settings dataframe:
settings = pd.read_csv("setting.csv", index_col= 0)
# Load data
project_folder = '/home/ziyz1701/storage/CD_airfoil/AIAA_2023_paper_plot/shear_layer_refinement/'
probe_data = np.genfromtxt(project_folder + 'shear_layer_profile/' + 'shear_layer_-050mm.csv', delimiter=',', skip_header=1)
probe = probe_data

#plot setting
linewidth_3D_result = eval(settings.at["linewidth_3D_result", settings.columns[0]])
linewidth_2D_result = eval(settings.at["linewidth_2D_result", settings.columns[0]])
linewidth_list = eval(settings.at["linewidth_list", settings.columns[0]])
marker_list = eval(settings.at["marker_list", settings.columns[0]])
marker_size = eval(settings.at["marker_size", settings.columns[0]])
fontsize_var = eval(settings.at["fontsize_var", settings.columns[0]])
plot_size_rectangle = eval(settings.at["plot_size_rectangle", settings.columns[0]])

# Input
minimum_voxel_size = 5.92e-05
viscosity = 1.44e-5
rho = 1.251
U_ref = 15.67

#Calc variable
wall_distance = probe[:,0]
relative_velocity_magnitude = probe[:,1]

relative_velocity_magnitude[0]=0
wall_distance[0]=0

cs = CubicSpline(wall_distance,relative_velocity_magnitude)
x_0 = 0
dU_dy = cs(x_0, 1)
tau = rho*viscosity*dU_dy
u_tau = np.sqrt(tau / rho)
print('wall shear :',tau)

# Calculate local y+ based on minimum voxel height
y_plus = minimum_voxel_size * u_tau / viscosity
print('dimensionless voxel size', y_plus)
h_intp = np.linspace(0,0.1,100000)
U_intp = cs(h_intp)
#Plot the viscous sublayer

u_viscous_layer = h_intp * u_tau / viscosity * u_tau
u_log_law = 1/0.41*np.log(h_intp * u_tau / viscosity) + 5.0

plt.figure(figsize=plot_size_rectangle)
# Enable LaTeX rendering
plt.rc('text', usetex=True)
# Use plain TeX font (Computer Modern) for legend
plt.rc('legend', fontsize=fontsize_var, fancybox=True, framealpha=0.7, edgecolor='black', facecolor='white', borderaxespad=0.5, labelspacing=0.5)
plt.rc('font', family='serif', weight='normal')

plt.plot(h_intp * u_tau / viscosity,u_viscous_layer,label=r'Viscous sublayer ($u^+=y^+$)',linestyle='--',color='k',linewidth=linewidth_list[0])
plt.plot(h_intp * u_tau / viscosity,u_log_law,label=r'Log law ($u^{+} = \frac{1}{\kappa}ln(y^{+})+C^{+}$)',linestyle='-.',color='k',linewidth=linewidth_list[0])
plt.plot(h_intp * u_tau / viscosity, U_intp / u_tau,color='k',linewidth=linewidth_list[0],label='DNS-SLRT')

# Add labels and legend
plt.ylabel(r'$U^+$', fontsize=fontsize_var, weight='bold')
plt.xlabel(r'$y^+$', fontsize=fontsize_var, weight='bold')
# Add legend
plt.legend(prop = { "size": fontsize_var })
plt.tick_params(axis='both', which='major', labelsize=fontsize_var)
ax = plt.gca()
ax.yaxis.offsetText.set(size=fontsize_var)
plt.yticks(fontsize=fontsize_var)
plt.xticks(fontsize=fontsize_var)

plt.grid(alpha=0.2, color='black')
plt.xscale('log')
plt.xlim([0.1, 100])
plt.ylim([0, 20])
plt.savefig(project_folder + 'shear_layer_-050mm.pdf', dpi=300)
print('saved in',project_folder)

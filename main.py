"""
Using Fits from Matlab to Get Dissociation Rates
"""

import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy import optimize
import numpy as np
import sys
from numpy import mean as m
import os
import os.path
import sys
import statistics
import glob
from HelperFunctions import get_intact_columns_constant_pressure
from HelperFunctions import get_intact_columns_constant_speed
from HelperFunctions import get_fitted_plots_constant_temperature
from HelperFunctions import get_fitted_plots_constant_pressure
from HelperFunctions import get_MATLABFIT_dissociation_rates
from HelperFunctions import get_average_shear_normal_stress_and_average_mu_constant_temperature
from HelperFunctions import get_average_shear_normal_stress_and_average_mu_constant_pressure
from mpl_toolkits import mplot3d as Axes3D
from HelperFunctions import plot_shear_stress_vs_normal_stress, plot_variation_in_mu
from HelperFunctions import get_dissociation_rates

Temperatures = ["400K", "500K", "600K", "700K"]
Pressures = ['2GPa', '3GPa', '4GPa', '5GPa']

Big_Dataframe_1ms = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/1msComparison.csv', index_col=False)
Big_Dataframe_10ms = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/10msComparison.csv', index_col=False)
Big_Dataframe_20ms = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/20msComparison.csv', index_col=False)
Big_Dataframe_30ms = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/30msComparison.csv', index_col=False)
Big_Dataframe_40ms = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/40msComparison.csv', index_col=False)
Big_Dataframe_50ms = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/50msComparison.csv', index_col=False)

Big_Dataframe_1GPa = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/1msComparison.csv', index_col=False)
Big_Dataframe_2GPa = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/10msComparison.csv', index_col=False)
Big_Dataframe_3GPa = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/20msComparison.csv', index_col=False)
Big_Dataframe_4GPa = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/30msComparison.csv', index_col=False)
Big_Dataframe_5GPa = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/40msComparison.csv', index_col=False)

Timestep_1ms = Big_Dataframe_1ms['Timestep'].to_list()
Timestep_10ms = Big_Dataframe_10ms['Timestep'].to_list()
Timestep_20ms = Big_Dataframe_20ms['Timestep'].to_list()
Timestep_30ms = Big_Dataframe_30ms['Timestep'].to_list()
Timestep_40ms = Big_Dataframe_40ms['Timestep'].to_list()
Timestep_50ms = Big_Dataframe_50ms['Timestep'].to_list()

Timestep_1GPa = Big_Dataframe_1GPa['Timestep']
Timestep_2GPa = Big_Dataframe_1GPa['Timestep']
Timestep_3GPa = Big_Dataframe_1GPa['Timestep']
Timestep_4GPa = Big_Dataframe_1GPa['Timestep']
Timestep_5GPa = Big_Dataframe_1GPa['Timestep']

# ################# Getting Dissociation Rates for Controlled Speed Using MATLAB Fits ##########
#
#Dissociation_Rate_1ms_1GPa, LogRate_1ms_1GPa = get_MATLABFIT_dissociation_rates(Timestep_1ms, -0.4884)
Dissociation_Rate_1ms_2GPa, LogRate_1ms_2GPa = get_MATLABFIT_dissociation_rates(Timestep_1ms, -0.4788, Cutoff=0.2)
Dissociation_Rate_1ms_3GPa, LogRate_1ms_3GPa = get_MATLABFIT_dissociation_rates(Timestep_1ms, -0.8493, Cutoff=0.2)
Dissociation_Rate_1ms_4GPa, LogRate_1ms_4GPa = get_MATLABFIT_dissociation_rates(Timestep_1ms, -1.359, Cutoff=0.2)
Dissociation_Rate_1ms_5GPa, LogRate_1ms_5GPa = get_MATLABFIT_dissociation_rates(Timestep_1ms, -1.479, Cutoff=0.2)

Dissociation_Rate_10ms_1GPa, LogRate_10ms_1GPa = get_MATLABFIT_dissociation_rates(Timestep_10ms, -1.244, Cutoff=0.2)
#Dissociation_Rate_10ms_2GPa, LogRate_10ms_2GPa = get_MATLABFIT_dissociation_rates(Timestep_10ms)
Dissociation_Rate_10ms_3GPa, LogRate_10ms_3GPa = get_MATLABFIT_dissociation_rates(Timestep_10ms, -4.645, Cutoff=0.2)
Dissociation_Rate_10ms_4GPa, LogRate_10ms_4GPa = get_MATLABFIT_dissociation_rates(Timestep_10ms, -4.825, Cutoff=0.2)
Dissociation_Rate_10ms_5GPa, LogRate_10ms_5GPa = get_MATLABFIT_dissociation_rates(Timestep_10ms, -7.973, Cutoff=0.2)

Dissociation_Rate_20ms_1GPa, LogRate_20ms_1GPa = get_MATLABFIT_dissociation_rates(Timestep_20ms, -1.409, Cutoff=0.2)
Dissociation_Rate_20ms_2GPa, LogRate_20ms_2GPa = get_MATLABFIT_dissociation_rates(Timestep_20ms, -4.373, Cutoff=0.2)
Dissociation_Rate_20ms_3GPa, LogRate_20ms_3GPa = get_MATLABFIT_dissociation_rates(Timestep_20ms, -5.679, Cutoff=0.2)
Dissociation_Rate_20ms_4GPa, LogRate_20ms_4GPa = get_MATLABFIT_dissociation_rates(Timestep_20ms, -7.833, Cutoff=0.2)
Dissociation_Rate_20ms_5GPa, LogRate_20ms_5GPa = get_MATLABFIT_dissociation_rates(Timestep_20ms, -10.58, Cutoff=0.2)

#Dissociation_Rate_30ms_1GPa, LogRate_30ms_1GPa = get_MATLABFIT_dissociation_rates(Timestep_30ms, -1.41)
Dissociation_Rate_30ms_2GPa, LogRate_30ms_2GPa = get_MATLABFIT_dissociation_rates(Timestep_30ms, -5.955, Cutoff=0.2)
#Dissociation_Rate_30ms_3GPa, LogRate_30ms_3GPa = get_MATLABFIT_dissociation_rates(Timestep_30ms, -0.8493)
Dissociation_Rate_30ms_4GPa, LogRate_30ms_4GPa = get_MATLABFIT_dissociation_rates(Timestep_30ms, -10.58, Cutoff=0.2)
Dissociation_Rate_30ms_5GPa, LogRate_30ms_5GPa = get_MATLABFIT_dissociation_rates(Timestep_30ms, -11.41, Cutoff=0.2)

Dissociation_Rate_40ms_1GPa, LogRate_40ms_1GPa = get_MATLABFIT_dissociation_rates(Timestep_40ms, -0.4884, Cutoff=0.2)
Dissociation_Rate_40ms_2GPa, LogRate_40ms_2GPa = get_MATLABFIT_dissociation_rates(Timestep_40ms, -7.186, Cutoff=0.2)
Dissociation_Rate_40ms_3GPa, LogRate_40ms_3GPa = get_MATLABFIT_dissociation_rates(Timestep_40ms, -9.287, Cutoff=0.2)
Dissociation_Rate_40ms_4GPa, LogRate_40ms_4GPa = get_MATLABFIT_dissociation_rates(Timestep_40ms, -15.24, Cutoff=0.2)
Dissociation_Rate_40ms_5GPa, LogRate_40ms_5GPa = get_MATLABFIT_dissociation_rates(Timestep_40ms, -17.83, Cutoff=0.2)
#
#Dissociation_Rate_50ms_1GPa, LogRate_50ms_1GPa = get_MATLABFIT_dissociation_rates(Timestep_50ms, -0.4884)
Dissociation_Rate_50ms_2GPa, LogRate_50ms_2GPa = get_MATLABFIT_dissociation_rates(Timestep_50ms, -9.587, Cutoff=0.2)
Dissociation_Rate_50ms_3GPa, LogRate_50ms_3GPa = get_MATLABFIT_dissociation_rates(Timestep_50ms, -14.35, Cutoff=0.2)
Dissociation_Rate_50ms_4GPa, LogRate_50ms_4GPa = get_MATLABFIT_dissociation_rates(Timestep_50ms, -18.29, Cutoff=0.2)
Dissociation_Rate_50ms_5GPa, LogRate_50ms_5GPa = get_MATLABFIT_dissociation_rates(Timestep_50ms, -23.16, Cutoff=0.2)


#print(f"Dissociation Rate at 1ms, 1GPa is {Dissociation_Rate_1ms_1GPa}, log of dissociation rate is {LogRate_300K_1GPa}")
print(f"Dissociation Rate at 1ms, 2GPa is {Dissociation_Rate_1ms_2GPa}, log of dissociation rate is {LogRate_1ms_2GPa}")
print(f"Dissociation Rate at 1ms, 3GPa is {Dissociation_Rate_1ms_3GPa}, log of dissociation rate is {LogRate_1ms_3GPa}")
print(f"Dissociation Rate at 1ms, 4GPa is {Dissociation_Rate_1ms_4GPa}, log of dissociation rate is {LogRate_1ms_4GPa}")
print(f"Dissociation Rate at 1ms, 5GPa is {Dissociation_Rate_1ms_5GPa}, log of dissociation rate is {LogRate_1ms_5GPa}")

print(f"Dissociation Rate at 10ms, 1GPa is {Dissociation_Rate_10ms_1GPa}, log of dissociation rate is {LogRate_10ms_1GPa}")
#print(f"Dissociation Rate at 10ms, 2GPa is {Dissociation_Rate_10ms_2GPa}, log of dissociation rate is {LogRate_10ms_2GPa}")
print(f"Dissociation Rate at 10ms, 3GPa is {Dissociation_Rate_10ms_3GPa}, log of dissociation rate is {LogRate_10ms_3GPa}")
print(f"Dissociation Rate at 10ms, 4GPa is {Dissociation_Rate_10ms_4GPa}, log of dissociation rate is {LogRate_10ms_4GPa}")
print(f"Dissociation Rate at 10ms, 5GPa is {Dissociation_Rate_10ms_5GPa}, log of dissociation rate is {LogRate_10ms_5GPa}")
#
print(f"Dissociation Rate at 20ms, 1GPa is {Dissociation_Rate_20ms_1GPa}, log of dissociation rate is {LogRate_20ms_1GPa}")
#print(f"Dissociation Rate at 20ms, 2GPa is {Dissociation_Rate_20ms_2GPa}, log of dissociation rate is {LogRate_20ms_2GPa}")
print(f"Dissociation Rate at 20ms, 3GPa is {Dissociation_Rate_20ms_3GPa}, log of dissociation rate is {LogRate_20ms_3GPa}")
print(f"Dissociation Rate at 20ms, 4GPa is {Dissociation_Rate_20ms_4GPa}, log of dissociation rate is {LogRate_20ms_4GPa}")
print(f"Dissociation Rate at 20ms, 5GPa is {Dissociation_Rate_20ms_5GPa}, log of dissociation rate is {LogRate_20ms_5GPa}")

#print(f"Dissociation Rate at 30ms, 1GPa is {Dissociation_Rate_30ms_1GPa}, log of dissociation rate is {LogRate_30ms_1GPa}")
print(f"Dissociation Rate at 30ms, 2GPa is {Dissociation_Rate_30ms_2GPa}, log of dissociation rate is {LogRate_30ms_2GPa}")
#print(f"Dissociation Rate at 30ms, 3GPa is {Dissociation_Rate_30ms_3GPa}, log of dissociation rate is {LogRate_30ms_3GPa}")
print(f"Dissociation Rate at 30ms, 4GPa is {Dissociation_Rate_30ms_4GPa}, log of dissociation rate is {LogRate_30ms_4GPa}")
print(f"Dissociation Rate at 30ms, 5GPa is {Dissociation_Rate_30ms_5GPa}, log of dissociation rate is {LogRate_30ms_5GPa}")

print(f"Dissociation Rate at 40ms, 1GPa is {Dissociation_Rate_40ms_1GPa}, log of dissociation rate is {LogRate_40ms_1GPa}")
print(f"Dissociation Rate at 40ms, 2GPa is {Dissociation_Rate_40ms_2GPa}, log of dissociation rate is {LogRate_40ms_2GPa}")
print(f"Dissociation Rate at 40ms, 3GPa is {Dissociation_Rate_40ms_3GPa}, log of dissociation rate is {LogRate_40ms_3GPa}")
print(f"Dissociation Rate at 40ms, 4GPa is {Dissociation_Rate_40ms_4GPa}, log of dissociation rate is {LogRate_40ms_4GPa}")
print(f"Dissociation Rate at 40ms, 5GPa is {Dissociation_Rate_40ms_5GPa}, log of dissociation rate is {LogRate_40ms_5GPa}")

#print(f"Dissociation Rate at 50ms, 1GPa is {Dissociation_Rate_50ms_1GPa}, log of dissociation rate is {LogRate_50ms_1GPa}")
print(f"Dissociation Rate at 50ms, 2GPa is {Dissociation_Rate_50ms_2GPa}, log of dissociation rate is {LogRate_50ms_2GPa}")
print(f"Dissociation Rate at 50ms, 3GPa is {Dissociation_Rate_50ms_3GPa}, log of dissociation rate is {LogRate_50ms_3GPa}")
print(f"Dissociation Rate at 50ms, 4GPa is {Dissociation_Rate_50ms_4GPa}, log of dissociation rate is {LogRate_50ms_4GPa}")
print(f"Dissociation Rate at 50ms, 5GPa is {Dissociation_Rate_50ms_5GPa}, log of dissociation rate is {LogRate_50ms_5GPa}")


# ########### Plotting Log of Dissociation Rates Against 1000/T to get 2D ################
#
# Log_Dissociation_Rates_List_300K = np.array([LogRate_300K_1GPa, LogRate_300K_2GPa, LogRate_300K_3GPa,
#                                      LogRate_300K_4GPa, LogRate_300K_5GPa])
# Log_Dissociation_Rates_List_400K = np.array([LogRate_400K_1GPa, LogRate_400K_2GPa, LogRate_400K_3GPa,
#                                      LogRate_400K_4GPa, LogRate_400K_5GPa])
# Log_Dissociation_Rates_List_500K = [LogRate_500K_1GPa, LogRate_500K_2GPa, LogRate_500K_3GPa,
#                                      LogRate_500K_4GPa, LogRate_500K_5GPa]
# Log_Dissociation_Rates_List_600K = [LogRate_600K_1GPa, LogRate_600K_2GPa, LogRate_600K_3GPa,
#                                      LogRate_600K_4GPa, LogRate_600K_5GPa]
# Log_Dissociation_Rates_List_700K = [LogRate_700K_1GPa, LogRate_700K_2GPa, LogRate_700K_3GPa,
#                                      LogRate_700K_4GPa, LogRate_700K_5GPa]
#
#
#
# Temperatures = [300, 400, 500, 600, 700]
# Inverse_Temperatures = np.array([1000/x for x in Temperatures])
# #Inverse_Temperatures = sorted(Inverse_Temperatures)
#
# trend300K = np.polyfit(Inverse_Temperatures, Log_Dissociation_Rates_List_300K, 1)
# trend400K = np.polyfit(Inverse_Temperatures, Log_Dissociation_Rates_List_400K, 1)
# trend500K = np.polyfit(Inverse_Temperatures, Log_Dissociation_Rates_List_500K, 1)
# trend600K = np.polyfit(Inverse_Temperatures, Log_Dissociation_Rates_List_600K, 1)
# trend700K = np.polyfit(Inverse_Temperatures, Log_Dissociation_Rates_List_700K, 1)
#
# #def plot_lnrate_vs_1000overT():
# fig1, ax1 = plt.subplots()
# ax1.set_title('Dissociation Rates against Inverse of Temperatures')
# ax1.set_xlabel('1000/T (K-1)')
# ax1.set_ylabel('ln(Rate) (ns-1)')
# ax1.set_xlim([1.4, 3.4])
# ax1.set_ylim([0, 4.5])
# ax1.scatter(Inverse_Temperatures, Log_Dissociation_Rates_List_300K)
# ax1.scatter(Inverse_Temperatures, Log_Dissociation_Rates_List_400K)
# ax1.scatter(Inverse_Temperatures, Log_Dissociation_Rates_List_500K)
# ax1.scatter(Inverse_Temperatures, Log_Dissociation_Rates_List_600K)
# ax1.scatter(Inverse_Temperatures, Log_Dissociation_Rates_List_700K)
#
# Fit300K = np.poly1d(trend300K)
# Fit400K = np.poly1d(trend400K)
# Fit500K = np.poly1d(trend500K)
# Fit600K = np.poly1d(trend600K)
# Fit700K = np.poly1d(trend700K)
#
# ax1.plot(Inverse_Temperatures, Fit300K(Inverse_Temperatures), label='300K')
# ax1.plot(Inverse_Temperatures, Fit400K(Inverse_Temperatures), label='400K')
# ax1.plot(Inverse_Temperatures, Fit500K(Inverse_Temperatures), label='500K')
# ax1.plot(Inverse_Temperatures, Fit600K(Inverse_Temperatures), label='600K')
# ax1.plot(Inverse_Temperatures, Fit700K(Inverse_Temperatures), label='700K')
# ax1.legend()
# plt.show()

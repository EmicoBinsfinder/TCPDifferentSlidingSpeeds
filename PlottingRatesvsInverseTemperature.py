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
from HelperFunctions import get_intact_columns_constant_temperature
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

Big_Dataframe_300K = get_intact_columns_constant_temperature("300K", Pressures)
Big_Dataframe_400K = get_intact_columns_constant_temperature("400K", Pressures)
Big_Dataframe_500K = get_intact_columns_constant_temperature("500K", Pressures)
Big_Dataframe_600K = get_intact_columns_constant_temperature("600K", Pressures)
Big_Dataframe_700K = get_intact_columns_constant_temperature("700K", Pressures)

########### Getting Fitted Plots/Equations At Constant Temperatures #####################

Timestep_List_300K, fitted_function_OneGPa_300K, fitted_functionTwoGPa_300K, fitted_functionThreeGPa_300K, fitted_functionFourGPa_300K, fitted_functionFiveGPa_300K \
    = get_fitted_plots_constant_temperature(Big_Dataframe_300K, "10ms", "300K")
Timestep_List_400K, fitted_function_OneGPa_400K, fitted_functionTwoGPa_400K, fitted_functionThreeGPa_400K, fitted_functionFourGPa_400K, fitted_functionFiveGPa_400K \
    = get_fitted_plots_constant_temperature(Big_Dataframe_400K, "10ms", "400K")
Timestep_List_500K, fitted_function_OneGPa_500K, fitted_functionTwoGPa_500K, fitted_functionThreeGPa_500K, fitted_functionFourGPa_500K, fitted_functionFiveGPa_500K \
   = get_fitted_plots_constant_temperature(Big_Dataframe_500K, "10ms", "500K")
Timestep_List_600K, fitted_function_OneGPa_600K, fitted_functionTwoGPa_600K, fitted_functionThreeGPa_600K, fitted_functionFourGPa_600K, fitted_functionFiveGPa_600K \
    = get_fitted_plots_constant_temperature(Big_Dataframe_600K, "10ms", "600K")
Timestep_List_700K, fitted_function_OneGPa_700K, fitted_functionTwoGPa_700K, fitted_functionThreeGPa_700K, fitted_functionFourGPa_700K, fitted_functionFiveGPa_700K \
    = get_fitted_plots_constant_temperature(Big_Dataframe_700K, "10ms", "700K")

################# Getting Dissociation Rates for Controlled Temperature Using MATLAB Fits ##########

Dissociation_Rate_300K_1GPa, LogRate_300K_1GPa = get_MATLABFIT_dissociation_rates(Timestep_List_300K, -0.4884)
Dissociation_Rate_300K_2GPa, LogRate_300K_2GPa = get_MATLABFIT_dissociation_rates(Timestep_List_300K, -2.427)
Dissociation_Rate_300K_3GPa, LogRate_300K_3GPa = get_MATLABFIT_dissociation_rates(Timestep_List_300K, -3.233)
Dissociation_Rate_300K_4GPa, LogRate_300K_4GPa = get_MATLABFIT_dissociation_rates(Timestep_List_300K, -6.051)
Dissociation_Rate_300K_5GPa, LogRate_300K_5GPa = get_MATLABFIT_dissociation_rates(Timestep_List_300K, -8.056)

Dissociation_Rate_400K_1GPa, LogRate_400K_1GPa = get_MATLABFIT_dissociation_rates(Timestep_List_400K, -0.6863)
Dissociation_Rate_400K_2GPa, LogRate_400K_2GPa = get_MATLABFIT_dissociation_rates(Timestep_List_400K, -1.722)
Dissociation_Rate_400K_3GPa, LogRate_400K_3GPa = get_MATLABFIT_dissociation_rates(Timestep_List_400K, -3.469)
Dissociation_Rate_400K_4GPa, LogRate_400K_4GPa = get_MATLABFIT_dissociation_rates(Timestep_List_400K, -4.916)
Dissociation_Rate_400K_5GPa, LogRate_400K_5GPa = get_MATLABFIT_dissociation_rates(Timestep_List_400K, -6.849)

Dissociation_Rate_500K_1GPa, LogRate_500K_1GPa = get_MATLABFIT_dissociation_rates(Timestep_List_500K, -1.298)
Dissociation_Rate_500K_2GPa, LogRate_500K_2GPa = get_MATLABFIT_dissociation_rates(Timestep_List_500K, -3.82)
Dissociation_Rate_500K_3GPa, LogRate_500K_3GPa = get_MATLABFIT_dissociation_rates(Timestep_List_500K, -4.859)
Dissociation_Rate_500K_4GPa, LogRate_500K_4GPa = get_MATLABFIT_dissociation_rates(Timestep_List_500K, -8.563)
Dissociation_Rate_500K_5GPa, LogRate_500K_5GPa = get_MATLABFIT_dissociation_rates(Timestep_List_500K, -10.57)

Dissociation_Rate_600K_1GPa, LogRate_600K_1GPa = get_MATLABFIT_dissociation_rates(Timestep_List_600K, -4.281)
Dissociation_Rate_600K_2GPa, LogRate_600K_2GPa = get_MATLABFIT_dissociation_rates(Timestep_List_600K, -8.479)
Dissociation_Rate_600K_3GPa, LogRate_600K_3GPa = get_MATLABFIT_dissociation_rates(Timestep_List_600K, -10.09)
Dissociation_Rate_600K_4GPa, LogRate_600K_4GPa = get_MATLABFIT_dissociation_rates(Timestep_List_600K, -11.67)
Dissociation_Rate_600K_5GPa, LogRate_600K_5GPa = get_MATLABFIT_dissociation_rates(Timestep_List_600K, -13.95)

Dissociation_Rate_700K_1GPa, LogRate_700K_1GPa = get_MATLABFIT_dissociation_rates(Timestep_List_700K, -5.933)
Dissociation_Rate_700K_2GPa, LogRate_700K_2GPa = get_MATLABFIT_dissociation_rates(Timestep_List_700K, -10.27)
Dissociation_Rate_700K_3GPa, LogRate_700K_3GPa = get_MATLABFIT_dissociation_rates(Timestep_List_700K, -10.63)
Dissociation_Rate_700K_4GPa, LogRate_700K_4GPa = get_MATLABFIT_dissociation_rates(Timestep_List_700K, -16.94)
Dissociation_Rate_700K_5GPa, LogRate_700K_5GPa = get_MATLABFIT_dissociation_rates(Timestep_List_700K, -28.76)

########### Plotting Log of Dissociation Rates Against 1000/T to get 2D ################

Log_Dissociation_Rates_List_300K = np.array([LogRate_300K_1GPa, LogRate_300K_2GPa, LogRate_300K_3GPa,
                                     LogRate_300K_4GPa, LogRate_300K_5GPa])
Log_Dissociation_Rates_List_400K = np.array([LogRate_400K_1GPa, LogRate_400K_2GPa, LogRate_400K_3GPa,
                                     LogRate_400K_4GPa, LogRate_400K_5GPa])
Log_Dissociation_Rates_List_500K = [LogRate_500K_1GPa, LogRate_500K_2GPa, LogRate_500K_3GPa,
                                     LogRate_500K_4GPa, LogRate_500K_5GPa]
Log_Dissociation_Rates_List_600K = [LogRate_600K_1GPa, LogRate_600K_2GPa, LogRate_600K_3GPa,
                                     LogRate_600K_4GPa, LogRate_600K_5GPa]
Log_Dissociation_Rates_List_700K = [LogRate_700K_1GPa, LogRate_700K_2GPa, LogRate_700K_3GPa,
                                     LogRate_700K_4GPa, LogRate_700K_5GPa]



Temperatures = [300, 400, 500, 600, 700]
Inverse_Temperatures = np.array([1000/x for x in Temperatures])
#Inverse_Temperatures = sorted(Inverse_Temperatures)

trend300K = np.polyfit(Inverse_Temperatures, Log_Dissociation_Rates_List_300K, 1)
trend400K = np.polyfit(Inverse_Temperatures, Log_Dissociation_Rates_List_400K, 1)
trend500K = np.polyfit(Inverse_Temperatures, Log_Dissociation_Rates_List_500K, 1)
trend600K = np.polyfit(Inverse_Temperatures, Log_Dissociation_Rates_List_600K, 1)
trend700K = np.polyfit(Inverse_Temperatures, Log_Dissociation_Rates_List_700K, 1)

#def plot_lnrate_vs_1000overT():
fig1, ax1 = plt.subplots()
ax1.set_title('Dissociation Rates against Inverse of Temperatures')
ax1.set_xlabel('1000/T (K-1)')
ax1.set_ylabel('ln(Rate) (ns-1)')
ax1.set_xlim([1.4, 3.4])
ax1.set_ylim([0, 4.5])
ax1.scatter(Inverse_Temperatures, Log_Dissociation_Rates_List_300K)
ax1.scatter(Inverse_Temperatures, Log_Dissociation_Rates_List_400K)
ax1.scatter(Inverse_Temperatures, Log_Dissociation_Rates_List_500K)
ax1.scatter(Inverse_Temperatures, Log_Dissociation_Rates_List_600K)
ax1.scatter(Inverse_Temperatures, Log_Dissociation_Rates_List_700K)

Fit300K = np.poly1d(trend300K)
Fit400K = np.poly1d(trend400K)
Fit500K = np.poly1d(trend500K)
Fit600K = np.poly1d(trend600K)
Fit700K = np.poly1d(trend700K)

ax1.plot(Inverse_Temperatures, Fit300K(Inverse_Temperatures), label='300K')
ax1.plot(Inverse_Temperatures, Fit400K(Inverse_Temperatures), label='400K')
ax1.plot(Inverse_Temperatures, Fit500K(Inverse_Temperatures), label='500K')
ax1.plot(Inverse_Temperatures, Fit600K(Inverse_Temperatures), label='600K')
ax1.plot(Inverse_Temperatures, Fit700K(Inverse_Temperatures), label='700K')
ax1.legend()
plt.show()

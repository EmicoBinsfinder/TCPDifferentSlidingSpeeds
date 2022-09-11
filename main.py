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


#
########## Getting Formatted Dataframes to Perform Analysis On ########################
Big_Dataframe_300K = get_intact_columns_constant_temperature("300K", Pressures)
Big_Dataframe_400K = get_intact_columns_constant_temperature("400K", Pressures)
Big_Dataframe_500K = get_intact_columns_constant_temperature("500K", Pressures)
Big_Dataframe_600K = get_intact_columns_constant_temperature("600K", Pressures)
Big_Dataframe_700K = get_intact_columns_constant_temperature("700K", Pressures)

Big_Dataframe_1GPa = get_intact_columns_constant_pressure('1GPa', Temperatures)
Big_Dataframe_2GPa = get_intact_columns_constant_pressure('2GPa', Temperatures)
Big_Dataframe_3GPa = get_intact_columns_constant_pressure('3GPa', Temperatures)
Big_Dataframe_4GPa = get_intact_columns_constant_pressure('4GPa', Temperatures)
Big_Dataframe_5GPa = get_intact_columns_constant_pressure('5GPa', Temperatures)

Big_Dataframe_1GPa.to_csv('C:/Users/eeo21/Documents/PhD/TCPDecompositionExperiments/Completed/AlphaFe/OneGPaComparison.csv')
Big_Dataframe_2GPa.to_csv('C:/Users/eeo21/Documents/PhD/TCPDecompositionExperiments/Completed/AlphaFe/TwoGPaComparison.csv')
Big_Dataframe_3GPa.to_csv('C:/Users/eeo21/Documents/PhD/TCPDecompositionExperiments/Completed/AlphaFe/ThreeGPaComparison.csv')
Big_Dataframe_4GPa.to_csv('C:/Users/eeo21/Documents/PhD/TCPDecompositionExperiments/Completed/AlphaFe/FourGPaComparison.csv')
Big_Dataframe_5GPa.to_csv('C:/Users/eeo21/Documents/PhD/TCPDecompositionExperiments/Completed/AlphaFe/FiveGPaComparison.csv')

Big_Dataframe_300K.to_csv('C:/Users/eeo21/Documents/PhD/TCPDecompositionExperiments/Completed/AlphaFe/300KComparison.csv')
Big_Dataframe_400K.to_csv('C:/Users/eeo21/Documents/PhD/TCPDecompositionExperiments/Completed/AlphaFe/400KComparison.csv')
Big_Dataframe_500K.to_csv('C:/Users/eeo21/Documents/PhD/TCPDecompositionExperiments/Completed/AlphaFe/500KComparison.csv')
Big_Dataframe_600K.to_csv('C:/Users/eeo21/Documents/PhD/TCPDecompositionExperiments/Completed/AlphaFe/600KComparison.csv')
Big_Dataframe_700K.to_csv('C:/Users/eeo21/Documents/PhD/TCPDecompositionExperiments/Completed/AlphaFe/700KComparison.csv')

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

########### Getting Fitted Plots/Equations At Constant Pressures #####################

Timestep_List_1GPa, fitted_function_300K_1GPa, fitted_function_400K_1GPa, fitted_function_500K_1GPa, fitted_function_600K_1GPa, fitted_function_700K_1GPa \
    = get_fitted_plots_constant_pressure(Big_Dataframe_1GPa, "10ms", "1GPa")
Timestep_List_2GPa, fitted_function_300K_2GPa, fitted_function_400K_2GPa, fitted_function_500K_2GPa, fitted_function_600K_2GPa, fitted_function_700K_2GPa \
    = get_fitted_plots_constant_pressure(Big_Dataframe_2GPa, "10ms", "2GPa")
Timestep_List_3GPa, fitted_function_300K_3GPa, fitted_function_400K_3GPa, fitted_function_500K_3GPa, fitted_function_600K_3GPa, fitted_function_700K_3GPa \
    = get_fitted_plots_constant_pressure(Big_Dataframe_3GPa, "10ms", "3GPa")
Timestep_List_4GPa, fitted_function_300K_4GPa, fitted_function_400K_4GPa, fitted_function_500K_4GPa, fitted_function_600K_4GPa, fitted_function_700K_4GPa \
    = get_fitted_plots_constant_pressure(Big_Dataframe_4GPa, "10ms", "4GPa")
Timestep_List_5GPa, fitted_function_300K_5GPa, fitted_function_400K_5GPa, fitted_function_500K_5GPa, fitted_function_600K_5GPa, fitted_function_700K_5GPa \
    = get_fitted_plots_constant_pressure(Big_Dataframe_5GPa, "10ms", "5GPa")

######## Getting Average Shear Stress, Friction Coefficient and Normal Stress #################

Average_Mu_List_300K, Average_Shear_Stress_List_300K, NormalStressMeans_300K = get_average_shear_normal_stress_and_average_mu_constant_temperature("300K", Pressures=Pressures)
Average_Mu_List_400K, Average_Shear_Stress_List_400K, NormalStressMeans_400K = get_average_shear_normal_stress_and_average_mu_constant_temperature("400K", Pressures=Pressures)
Average_Mu_List_500K, Average_Shear_Stress_List_500K, NormalStressMeans_500K = get_average_shear_normal_stress_and_average_mu_constant_temperature("500K", Pressures=Pressures)
Average_Mu_List_600K, Average_Shear_Stress_List_600K, NormalStressMeans_600K = get_average_shear_normal_stress_and_average_mu_constant_temperature("600K", Pressures=Pressures)
Average_Mu_List_700K, Average_Shear_Stress_List_700K, NormalStressMeans_700K = get_average_shear_normal_stress_and_average_mu_constant_temperature("700K", Pressures=Pressures)

# plot_shear_stress_vs_normal_stress(Average_Shear_Stress_List_300K, Average_Shear_Stress_List_400K, Average_Shear_Stress_List_500K, Average_Shear_Stress_List_600K, Average_Shear_Stress_List_700K,
#                                   "300k", "400K", "500K", "600K", "700K")

Average_Mu_List_1GPa, Average_Shear_Stress_List_1GPa, NormalStressMeans_1GPa = get_average_shear_normal_stress_and_average_mu_constant_pressure("1GPa", Temperatures)
Average_Mu_List_2GPa, Average_Shear_Stress_List_2GPa, NormalStressMeans_2GPa = get_average_shear_normal_stress_and_average_mu_constant_pressure("2GPa", Temperatures)
Average_Mu_List_3GPa, Average_Shear_Stress_List_3GPa, NormalStressMeans_3GPa = get_average_shear_normal_stress_and_average_mu_constant_pressure("3GPa", Temperatures)
Average_Mu_List_4GPa, Average_Shear_Stress_List_4GPa, NormalStressMeans_4GPa = get_average_shear_normal_stress_and_average_mu_constant_pressure("4GPa", Temperatures)
Average_Mu_List_5GPa, Average_Shear_Stress_List_5GPa, NormalStressMeans_5GPa = get_average_shear_normal_stress_and_average_mu_constant_pressure("5GPa", Temperatures)

#############  Getting Dissociation Rates (Doesn't Matter if you do it using fits for constant pressure or temperature ################
# Dissociation_Rates_1GPa = \
#     get_dissociation_rates(Timestep_List_1GPa, fitted_function_300K_1GPa, fitted_function_400K_1GPa, fitted_function_500K_1GPa, fitted_function_600K_1GPa, fitted_function_700K_1GPa, cutoff=30)
# Dissociation_Rates_2GPa = \
#     get_dissociation_rates(Timestep_List_2GPa, fitted_function_300K_2GPa, fitted_function_400K_2GPa, fitted_function_500K_2GPa, fitted_function_600K_2GPa, fitted_function_700K_2GPa, cutoff=20)
# Dissociation_Rates_3GPa = \
#     get_dissociation_rates(Timestep_List_3GPa, fitted_function_300K_3GPa, fitted_function_400K_3GPa, fitted_function_500K_3GPa, fitted_function_600K_3GPa, fitted_function_700K_3GPa, cutoff=20)
# Dissociation_Rates_4GPa = \
#     get_dissociation_rates(Timestep_List_4GPa, fitted_function_300K_4GPa, fitted_function_400K_4GPa, fitted_function_500K_4GPa, fitted_function_600K_4GPa, fitted_function_700K_4GPa, cutoff=20)
# Dissociation_Rates_5GPa = \
#     get_dissociation_rates(Timestep_List_5GPa, fitted_function_300K_5GPa, fitted_function_400K_5GPa, fitted_function_500K_5GPa, fitted_function_600K_4GPa, fitted_function_700K_4GPa, cutoff=20)

# Dissociation_Rates_300K = \
#     get_dissociation_rates(Timestep_List_300K, fitted_function_OneGPa_300K, fitted_functionTwoGPa_300K, fitted_functionThreeGPa_300K, fitted_functionFourGPa_300K, fitted_functionFiveGPa_300K, cutoff=35)
# Dissociation_Rates_400K = \
#     get_dissociation_rates(Timestep_List_400K, fitted_function_OneGPa_400K, fitted_functionTwoGPa_400K, fitted_functionThreeGPa_400K, fitted_functionFourGPa_400K, fitted_functionFiveGPa_400K, cutoff=20)
# Dissociation_Rates_500K = \
#     get_dissociation_rates(Timestep_List_500K, fitted_function_OneGPa_500K, fitted_functionTwoGPa_500K, fitted_functionThreeGPa_500K, fitted_functionFourGPa_500K, fitted_functionFiveGPa_500K, cutoff=20)
# Dissociation_Rates_600K = \
#     get_dissociation_rates(Timestep_List_600K, fitted_function_OneGPa_600K, fitted_functionTwoGPa_600K, fitted_functionThreeGPa_600K, fitted_functionFourGPa_600K, fitted_functionFiveGPa_600K, cutoff=20)
# Dissociation_Rates_700K = \
#     get_dissociation_rates(Timestep_List_700K, fitted_function_OneGPa_700K, fitted_functionTwoGPa_700K, fitted_functionThreeGPa_700K, fitted_functionFourGPa_700K, fitted_functionFiveGPa_700K, cutoff=20)

# print(Dissociation_Rates_1GPa)
# print(Dissociation_Rates_2GPa)
# print(Dissociation_Rates_3GPa)
# print(Dissociation_Rates_4GPa)
# print(Dissociation_Rates_5GPa)

print('#######################')

Average_Shear_Stress_List_1GPa = sorted(Average_Shear_Stress_List_1GPa)
Average_Shear_Stress_List_2GPa = sorted(Average_Shear_Stress_List_2GPa)
Average_Shear_Stress_List_3GPa = sorted(Average_Shear_Stress_List_3GPa)
Average_Shear_Stress_List_4GPa = sorted(Average_Shear_Stress_List_4GPa)
Average_Shear_Stress_List_5GPa = sorted(Average_Shear_Stress_List_5GPa)

Average_Shear_Stress_List_400K = sorted(Average_Shear_Stress_List_400K)
Average_Shear_Stress_List_500K = sorted(Average_Shear_Stress_List_500K)
Average_Shear_Stress_List_600K = sorted(Average_Shear_Stress_List_600K)
Average_Shear_Stress_List_700K = sorted(Average_Shear_Stress_List_700K)

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

print(f"Dissociation Rate at 300K, 1GPa is {Dissociation_Rate_300K_1GPa}, log of dissociation rate is {LogRate_300K_1GPa}")
print(f"Dissociation Rate at 300K, 2GPa is {Dissociation_Rate_300K_2GPa}, log of dissociation rate is {LogRate_300K_2GPa}")
print(f"Dissociation Rate at 300K, 3GPa is {Dissociation_Rate_300K_3GPa}, log of dissociation rate is {LogRate_300K_3GPa}")
print(f"Dissociation Rate at 300K, 4GPa is {Dissociation_Rate_300K_4GPa}, log of dissociation rate is {LogRate_300K_4GPa}")
print(f"Dissociation Rate at 300K, 5GPa is {Dissociation_Rate_300K_5GPa}, log of dissociation rate is {LogRate_300K_5GPa}")
print(f"Dissociation Rate at 400K, 1GPa is {Dissociation_Rate_400K_1GPa}, log of dissociation rate is {LogRate_400K_1GPa}")
print(f"Dissociation Rate at 400K, 2GPa is {Dissociation_Rate_400K_2GPa}, log of dissociation rate is {LogRate_400K_2GPa}")
print(f"Dissociation Rate at 400K, 3GPa is {Dissociation_Rate_400K_3GPa}, log of dissociation rate is {LogRate_400K_3GPa}")
print(f"Dissociation Rate at 400K, 4GPa is {Dissociation_Rate_400K_4GPa}, log of dissociation rate is {LogRate_400K_4GPa}")
print(f"Dissociation Rate at 400K, 5GPa is {Dissociation_Rate_400K_5GPa}, log of dissociation rate is {LogRate_400K_5GPa}")
print(f"Dissociation Rate at 500K, 1GPa is {Dissociation_Rate_500K_1GPa}, log of dissociation rate is {LogRate_500K_1GPa}")
print(f"Dissociation Rate at 500K, 2GPa is {Dissociation_Rate_500K_2GPa}, log of dissociation rate is {LogRate_500K_2GPa}")
print(f"Dissociation Rate at 500K, 3GPa is {Dissociation_Rate_500K_3GPa}, log of dissociation rate is {LogRate_500K_3GPa}")
print(f"Dissociation Rate at 500K, 4GPa is {Dissociation_Rate_500K_4GPa}, log of dissociation rate is {LogRate_500K_4GPa}")
print(f"Dissociation Rate at 500K, 5GPa is {Dissociation_Rate_500K_5GPa}, log of dissociation rate is {LogRate_500K_5GPa}")
print(f"Dissociation Rate at 600K, 1GPa is {Dissociation_Rate_600K_1GPa}, log of dissociation rate is {LogRate_600K_1GPa}")
print(f"Dissociation Rate at 600K, 2GPa is {Dissociation_Rate_600K_2GPa}, log of dissociation rate is {LogRate_600K_2GPa}")
print(f"Dissociation Rate at 600K, 3GPa is {Dissociation_Rate_600K_3GPa}, log of dissociation rate is {LogRate_600K_3GPa}")
print(f"Dissociation Rate at 600K, 4GPa is {Dissociation_Rate_600K_4GPa}, log of dissociation rate is {LogRate_600K_4GPa}")
print(f"Dissociation Rate at 600K, 5GPa is {Dissociation_Rate_600K_5GPa}, log of dissociation rate is {LogRate_600K_5GPa}")
print(f"Dissociation Rate at 700K, 1GPa is {Dissociation_Rate_700K_1GPa}, log of dissociation rate is {LogRate_700K_1GPa}")
print(f"Dissociation Rate at 700K, 2GPa is {Dissociation_Rate_700K_2GPa}, log of dissociation rate is {LogRate_700K_2GPa}")
print(f"Dissociation Rate at 700K, 3GPa is {Dissociation_Rate_700K_3GPa}, log of dissociation rate is {LogRate_700K_3GPa}")
print(f"Dissociation Rate at 700K, 4GPa is {Dissociation_Rate_700K_4GPa}, log of dissociation rate is {LogRate_700K_4GPa}")
print(f"Dissociation Rate at 700K, 5GPa is {Dissociation_Rate_700K_5GPa}, log of dissociation rate is {LogRate_700K_5GPa}")

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

print(Log_Dissociation_Rates_List_400K)
print(Log_Dissociation_Rates_List_500K)
print(Log_Dissociation_Rates_List_600K)
print(Log_Dissociation_Rates_List_700K)

print(Average_Shear_Stress_List_400K)
print(Average_Shear_Stress_List_500K)
print(Average_Shear_Stress_List_600K)
print(Average_Shear_Stress_List_700K)

figShearvLnk, axShearvLnk = plt.subplots()
x400K = Average_Shear_Stress_List_400K
x500K = Average_Shear_Stress_List_500K
x600K = Average_Shear_Stress_List_600K
x700K = Average_Shear_Stress_List_700K

a, b = np.polyfit(x400K, Log_Dissociation_Rates_List_400K , 1)
c, d = np.polyfit(x500K, Log_Dissociation_Rates_List_500K, 1)
e, f = np.polyfit(x600K, Log_Dissociation_Rates_List_600K, 1)
g, h = np.polyfit(x700K, Log_Dissociation_Rates_List_700K, 1)

axShearvLnk.set_title('Shear Stress vs Normal Stress at Different Temperatures')
axShearvLnk.set_xlabel('Shear Stress(ns-1)')
axShearvLnk.set_ylabel('Log Dissociation Rate (GPa)')
axShearvLnk.scatter(x400K, Log_Dissociation_Rates_List_400K)
axShearvLnk.scatter(x500K, Log_Dissociation_Rates_List_500K)
axShearvLnk.scatter(x600K, Log_Dissociation_Rates_List_600K)
axShearvLnk.scatter(x700K, Log_Dissociation_Rates_List_700K)
# axShearvLnk.plot(x, c * x + d, label=Temp2)
# axShearvLnk.plot(x, a * x + b, label=Temp1)
# axShearvLnk.plot(x, e * x + f, label=Temp3)
# axShearvLnk.plot(x, g * x + h, label=Temp4)
#axShearvLnk.legend()
plt.show()
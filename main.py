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
from HelperFunctions import get_average_shear_normal_stress_and_average_mu_constant_speed
from HelperFunctions import plot_shear_stress_vs_normal_stress_different_sliding_speeds
from mpl_toolkits import mplot3d as Axes3D
from HelperFunctions import plot_shear_stress_vs_normal_stress, plot_variation_in_mu
from HelperFunctions import get_dissociation_rates
from HelperFunctions import plot_variation_in_shear_stress_constant_speed
import Dissociation_Rates

Temperatures = ["400K", "500K", "600K", "700K"]
Pressures = ['2GPa', '3GPa', '4GPa', '5GPa']
Speeds = ['10ms', '20ms', '30ms', '40ms', '50ms']

Big_Dataframe_1ms = get_intact_columns_constant_speed("1ms", Pressures)
Big_Dataframe_10ms = get_intact_columns_constant_speed("10ms", Pressures)
Big_Dataframe_20ms = get_intact_columns_constant_speed("20ms", Pressures)
Big_Dataframe_30ms = get_intact_columns_constant_speed("30ms", Pressures)
Big_Dataframe_40ms = get_intact_columns_constant_speed("40ms", Pressures)
Big_Dataframe_50ms = get_intact_columns_constant_speed("50ms", Pressures)

Big_Dataframe_1GPa = get_intact_columns_constant_pressure('1GPa', Speeds)
Big_Dataframe_2GPa = get_intact_columns_constant_pressure('2GPa', Speeds)
Big_Dataframe_3GPa = get_intact_columns_constant_pressure('3GPa', Speeds)
Big_Dataframe_4GPa = get_intact_columns_constant_pressure('4GPa', Speeds)
Big_Dataframe_5GPa = get_intact_columns_constant_pressure('5GPa', Speeds)

Big_Dataframe_1ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/1msComparison.csv')
Big_Dataframe_10ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/10msComparison.csv')
Big_Dataframe_20ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/20msComparison.csv')
Big_Dataframe_30ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/30msComparison.csv')
Big_Dataframe_40ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/40msComparison.csv')
Big_Dataframe_50ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/50msComparison.csv')

Big_Dataframe_1GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/OneGPaComparison.csv')
Big_Dataframe_2GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/TwoGPaComparison.csv')
Big_Dataframe_3GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/ThreeGPaComparison.csv')
Big_Dataframe_4GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/FourGPaComparison.csv')
Big_Dataframe_5GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/FiveGPaComparison.csv')

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

Timestep_1GPa = Big_Dataframe_1GPa['Timestep'].to_list()
Timestep_2GPa = Big_Dataframe_2GPa['Timestep'].to_list()
Timestep_3GPa = Big_Dataframe_3GPa['Timestep'].to_list()
Timestep_4GPa = Big_Dataframe_4GPa['Timestep'].to_list()
Timestep_5GPa = Big_Dataframe_5GPa['Timestep'].to_list()

# ################# Getting Dissociation Rates for Controlled Speed Using MATLAB Fits ##########
Cutoff = None
Index = 1
# Index = 0 for actual rate, 1 for upper bound, 2 for lower bound

Dissociation_Rate_1ms_1GPa, LogRate_1ms_1GPa = get_MATLABFIT_dissociation_rates(Timestep_1ms, Dissociation_Rates.Dissociation_Rate_1ms_1GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_1ms_2GPa, LogRate_1ms_2GPa = get_MATLABFIT_dissociation_rates(Timestep_1ms, Dissociation_Rates.Dissociation_Rate_1ms_2GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_1ms_3GPa, LogRate_1ms_3GPa = get_MATLABFIT_dissociation_rates(Timestep_1ms, Dissociation_Rates.Dissociation_Rate_1ms_3GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_1ms_4GPa, LogRate_1ms_4GPa = get_MATLABFIT_dissociation_rates(Timestep_1ms, Dissociation_Rates.Dissociation_Rate_1ms_4GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_1ms_5GPa, LogRate_1ms_5GPa = get_MATLABFIT_dissociation_rates(Timestep_1ms, Dissociation_Rates.Dissociation_Rate_1ms_5GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rates_1ms = [Dissociation_Rate_1ms_1GPa, Dissociation_Rate_1ms_2GPa, Dissociation_Rate_1ms_3GPa, Dissociation_Rate_1ms_4GPa, Dissociation_Rate_1ms_5GPa]
Log_Rates_1ms = [LogRate_1ms_1GPa, LogRate_1ms_2GPa, LogRate_1ms_3GPa, LogRate_1ms_4GPa, LogRate_1ms_5GPa]

Dissociation_Rate_10ms_1GPa, LogRate_10ms_1GPa = get_MATLABFIT_dissociation_rates(Timestep_10ms, Dissociation_Rates.Dissociation_Rate_10ms_1GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_10ms_2GPa, LogRate_10ms_2GPa = get_MATLABFIT_dissociation_rates(Timestep_10ms, Dissociation_Rates.Dissociation_Rate_10ms_2GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_10ms_3GPa, LogRate_10ms_3GPa = get_MATLABFIT_dissociation_rates(Timestep_10ms, Dissociation_Rates.Dissociation_Rate_10ms_3GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_10ms_4GPa, LogRate_10ms_4GPa = get_MATLABFIT_dissociation_rates(Timestep_10ms, Dissociation_Rates.Dissociation_Rate_10ms_4GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_10ms_5GPa, LogRate_10ms_5GPa = get_MATLABFIT_dissociation_rates(Timestep_10ms, Dissociation_Rates.Dissociation_Rate_10ms_5GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rates_10ms = [Dissociation_Rate_10ms_1GPa, Dissociation_Rate_10ms_2GPa, Dissociation_Rate_10ms_3GPa, Dissociation_Rate_10ms_4GPa, Dissociation_Rate_10ms_5GPa]
Log_Rates_10ms = [LogRate_10ms_1GPa, LogRate_10ms_2GPa, LogRate_10ms_3GPa, LogRate_10ms_4GPa, LogRate_10ms_5GPa]

Dissociation_Rate_20ms_1GPa, LogRate_20ms_1GPa = get_MATLABFIT_dissociation_rates(Timestep_20ms, Dissociation_Rates.Dissociation_Rate_20ms_1GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_20ms_2GPa, LogRate_20ms_2GPa = get_MATLABFIT_dissociation_rates(Timestep_20ms, Dissociation_Rates.Dissociation_Rate_20ms_2GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_20ms_3GPa, LogRate_20ms_3GPa = get_MATLABFIT_dissociation_rates(Timestep_20ms, Dissociation_Rates.Dissociation_Rate_20ms_3GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_20ms_4GPa, LogRate_20ms_4GPa = get_MATLABFIT_dissociation_rates(Timestep_20ms, Dissociation_Rates.Dissociation_Rate_20ms_4GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_20ms_5GPa, LogRate_20ms_5GPa = get_MATLABFIT_dissociation_rates(Timestep_20ms, Dissociation_Rates.Dissociation_Rate_20ms_5GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rates_20ms = [Dissociation_Rate_20ms_1GPa, Dissociation_Rate_20ms_2GPa, Dissociation_Rate_20ms_3GPa, Dissociation_Rate_20ms_4GPa, Dissociation_Rate_20ms_5GPa]
Log_Rates_20ms = [LogRate_20ms_1GPa, LogRate_20ms_2GPa, LogRate_20ms_3GPa, LogRate_20ms_4GPa, LogRate_20ms_5GPa]

Dissociation_Rate_30ms_1GPa, LogRate_30ms_1GPa = get_MATLABFIT_dissociation_rates(Timestep_30ms, Dissociation_Rates.Dissociation_Rate_30ms_1GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_30ms_2GPa, LogRate_30ms_2GPa = get_MATLABFIT_dissociation_rates(Timestep_30ms, Dissociation_Rates.Dissociation_Rate_30ms_2GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_30ms_3GPa, LogRate_30ms_3GPa = get_MATLABFIT_dissociation_rates(Timestep_30ms, Dissociation_Rates.Dissociation_Rate_30ms_3GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_30ms_4GPa, LogRate_30ms_4GPa = get_MATLABFIT_dissociation_rates(Timestep_30ms, Dissociation_Rates.Dissociation_Rate_30ms_4GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_30ms_5GPa, LogRate_30ms_5GPa = get_MATLABFIT_dissociation_rates(Timestep_30ms, Dissociation_Rates.Dissociation_Rate_30ms_5GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rates_30ms = [Dissociation_Rate_30ms_1GPa, Dissociation_Rate_30ms_2GPa, Dissociation_Rate_30ms_3GPa, Dissociation_Rate_30ms_4GPa, Dissociation_Rate_30ms_5GPa]
Log_Rates_30ms = [LogRate_30ms_1GPa, LogRate_30ms_2GPa, LogRate_30ms_3GPa, LogRate_30ms_4GPa, LogRate_30ms_5GPa]

Dissociation_Rate_40ms_1GPa, LogRate_40ms_1GPa = get_MATLABFIT_dissociation_rates(Timestep_40ms, Dissociation_Rates.Dissociation_Rate_40ms_1GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_2GPa, LogRate_40ms_2GPa = get_MATLABFIT_dissociation_rates(Timestep_40ms, Dissociation_Rates.Dissociation_Rate_40ms_2GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_3GPa, LogRate_40ms_3GPa = get_MATLABFIT_dissociation_rates(Timestep_40ms, Dissociation_Rates.Dissociation_Rate_40ms_3GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_4GPa, LogRate_40ms_4GPa = get_MATLABFIT_dissociation_rates(Timestep_40ms, Dissociation_Rates.Dissociation_Rate_40ms_4GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_5GPa, LogRate_40ms_5GPa = get_MATLABFIT_dissociation_rates(Timestep_40ms, Dissociation_Rates.Dissociation_Rate_40ms_5GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rates_40ms = [Dissociation_Rate_40ms_1GPa, Dissociation_Rate_40ms_2GPa, Dissociation_Rate_40ms_3GPa, Dissociation_Rate_40ms_4GPa, Dissociation_Rate_40ms_5GPa]
Log_Rates_40ms = [LogRate_40ms_1GPa, LogRate_40ms_2GPa, LogRate_40ms_3GPa, LogRate_40ms_4GPa, LogRate_40ms_5GPa]

Dissociation_Rate_50ms_1GPa, LogRate_50ms_1GPa = get_MATLABFIT_dissociation_rates(Timestep_50ms, Dissociation_Rates.Dissociation_Rate_50ms_1GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_50ms_2GPa, LogRate_50ms_2GPa = get_MATLABFIT_dissociation_rates(Timestep_50ms, Dissociation_Rates.Dissociation_Rate_50ms_2GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_50ms_3GPa, LogRate_50ms_3GPa = get_MATLABFIT_dissociation_rates(Timestep_50ms, Dissociation_Rates.Dissociation_Rate_50ms_3GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_50ms_4GPa, LogRate_50ms_4GPa = get_MATLABFIT_dissociation_rates(Timestep_50ms, Dissociation_Rates.Dissociation_Rate_50ms_4GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_50ms_5GPa, LogRate_50ms_5GPa = get_MATLABFIT_dissociation_rates(Timestep_50ms, Dissociation_Rates.Dissociation_Rate_50ms_5GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rates_50ms = [Dissociation_Rate_50ms_1GPa, Dissociation_Rate_50ms_2GPa, Dissociation_Rate_50ms_3GPa, Dissociation_Rate_50ms_4GPa, Dissociation_Rate_50ms_5GPa]
Log_Rates_50ms = [LogRate_50ms_1GPa, LogRate_50ms_2GPa, LogRate_50ms_3GPa, LogRate_50ms_4GPa, LogRate_50ms_5GPa]

print(f"Dissociation Rate at 1ms, 1GPa is {Dissociation_Rate_1ms_1GPa}, log of dissociation rate is {LogRate_1ms_1GPa}")
print(f"Dissociation Rate at 1ms, 2GPa is {Dissociation_Rate_1ms_2GPa}, log of dissociation rate is {LogRate_1ms_2GPa}")
print(f"Dissociation Rate at 1ms, 3GPa is {Dissociation_Rate_1ms_3GPa}, log of dissociation rate is {LogRate_1ms_3GPa}")
print(f"Dissociation Rate at 1ms, 4GPa is {Dissociation_Rate_1ms_4GPa}, log of dissociation rate is {LogRate_1ms_4GPa}")
print(f"Dissociation Rate at 1ms, 5GPa is {Dissociation_Rate_1ms_5GPa}, log of dissociation rate is {LogRate_1ms_5GPa}")

print(f"Dissociation Rate at 10ms, 1GPa is {Dissociation_Rate_10ms_1GPa}, log of dissociation rate is {LogRate_10ms_1GPa}")
print(f"Dissociation Rate at 10ms, 2GPa is {Dissociation_Rate_10ms_2GPa}, log of dissociation rate is {LogRate_10ms_2GPa}")
print(f"Dissociation Rate at 10ms, 3GPa is {Dissociation_Rate_10ms_3GPa}, log of dissociation rate is {LogRate_10ms_3GPa}")
print(f"Dissociation Rate at 10ms, 4GPa is {Dissociation_Rate_10ms_4GPa}, log of dissociation rate is {LogRate_10ms_4GPa}")
print(f"Dissociation Rate at 10ms, 5GPa is {Dissociation_Rate_10ms_5GPa}, log of dissociation rate is {LogRate_10ms_5GPa}")

print(f"Dissociation Rate at 20ms, 1GPa is {Dissociation_Rate_20ms_1GPa}, log of dissociation rate is {LogRate_20ms_1GPa}")
print(f"Dissociation Rate at 20ms, 2GPa is {Dissociation_Rate_20ms_2GPa}, log of dissociation rate is {LogRate_20ms_2GPa}")
print(f"Dissociation Rate at 20ms, 3GPa is {Dissociation_Rate_20ms_3GPa}, log of dissociation rate is {LogRate_20ms_3GPa}")
print(f"Dissociation Rate at 20ms, 4GPa is {Dissociation_Rate_20ms_4GPa}, log of dissociation rate is {LogRate_20ms_4GPa}")
print(f"Dissociation Rate at 20ms, 5GPa is {Dissociation_Rate_20ms_5GPa}, log of dissociation rate is {LogRate_20ms_5GPa}")

print(f"Dissociation Rate at 30ms, 1GPa is {Dissociation_Rate_30ms_1GPa}, log of dissociation rate is {LogRate_30ms_1GPa}")
print(f"Dissociation Rate at 30ms, 2GPa is {Dissociation_Rate_30ms_2GPa}, log of dissociation rate is {LogRate_30ms_2GPa}")
print(f"Dissociation Rate at 30ms, 3GPa is {Dissociation_Rate_30ms_3GPa}, log of dissociation rate is {LogRate_30ms_3GPa}")
print(f"Dissociation Rate at 30ms, 4GPa is {Dissociation_Rate_30ms_4GPa}, log of dissociation rate is {LogRate_30ms_4GPa}")
print(f"Dissociation Rate at 30ms, 5GPa is {Dissociation_Rate_30ms_5GPa}, log of dissociation rate is {LogRate_30ms_5GPa}")

print(f"Dissociation Rate at 40ms, 1GPa is {Dissociation_Rate_40ms_1GPa}, log of dissociation rate is {LogRate_40ms_1GPa}")
print(f"Dissociation Rate at 40ms, 2GPa is {Dissociation_Rate_40ms_2GPa}, log of dissociation rate is {LogRate_40ms_2GPa}")
print(f"Dissociation Rate at 40ms, 3GPa is {Dissociation_Rate_40ms_3GPa}, log of dissociation rate is {LogRate_40ms_3GPa}")
print(f"Dissociation Rate at 40ms, 4GPa is {Dissociation_Rate_40ms_4GPa}, log of dissociation rate is {LogRate_40ms_4GPa}")
print(f"Dissociation Rate at 40ms, 5GPa is {Dissociation_Rate_40ms_5GPa}, log of dissociation rate is {LogRate_40ms_5GPa}")

print(f"Dissociation Rate at 50ms, 1GPa is {Dissociation_Rate_50ms_1GPa}, log of dissociation rate is {LogRate_50ms_1GPa}")
print(f"Dissociation Rate at 50ms, 2GPa is {Dissociation_Rate_50ms_2GPa}, log of dissociation rate is {LogRate_50ms_2GPa}")
print(f"Dissociation Rate at 50ms, 3GPa is {Dissociation_Rate_50ms_3GPa}, log of dissociation rate is {LogRate_50ms_3GPa}")
print(f"Dissociation Rate at 50ms, 4GPa is {Dissociation_Rate_50ms_4GPa}, log of dissociation rate is {LogRate_50ms_4GPa}")
print(f"Dissociation Rate at 50ms, 5GPa is {Dissociation_Rate_50ms_5GPa}, log of dissociation rate is {LogRate_50ms_5GPa}")

Log_Rates_1GPa = [LogRate_1ms_1GPa, LogRate_10ms_1GPa, LogRate_20ms_1GPa, LogRate_30ms_1GPa, LogRate_40ms_1GPa, LogRate_50ms_1GPa]
Log_Rates_2GPa = [LogRate_1ms_2GPa, LogRate_10ms_2GPa, LogRate_20ms_2GPa, LogRate_30ms_2GPa, LogRate_40ms_2GPa, LogRate_50ms_2GPa]
Log_Rates_3GPa = [LogRate_1ms_3GPa, LogRate_10ms_3GPa, LogRate_20ms_3GPa, LogRate_30ms_3GPa, LogRate_40ms_3GPa, LogRate_50ms_3GPa]
Log_Rates_4GPa = [LogRate_1ms_4GPa, LogRate_10ms_4GPa, LogRate_20ms_4GPa, LogRate_30ms_4GPa, LogRate_40ms_4GPa, LogRate_50ms_4GPa]
Log_Rates_5GPa = [LogRate_1ms_5GPa, LogRate_10ms_5GPa, LogRate_20ms_5GPa, LogRate_30ms_5GPa, LogRate_40ms_5GPa, LogRate_50ms_5GPa]

######## Getting Average Shear Stress, Friction Coefficient and Normal Stress #################

EquilibriumFactor = 0  # How many rows (out of 99) to ignore before calculating shear stress/friction coefficient, as it won't stabilise until after a certain number of timesteps

Average_Shear_Stress_List_1ms, Average_Mu_List_1ms, NormalStressMeans_1ms = get_average_shear_normal_stress_and_average_mu_constant_speed("1ms", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor)
Average_Shear_Stress_List_10ms, Average_Mu_List_10ms, NormalStressMeans_10ms = get_average_shear_normal_stress_and_average_mu_constant_speed("10ms", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor)
Average_Shear_Stress_List_20ms, Average_Mu_List_20ms, NormalStressMeans_20ms = get_average_shear_normal_stress_and_average_mu_constant_speed("20ms", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor)
Average_Shear_Stress_List_30ms, Average_Mu_List_30ms, NormalStressMeans_30ms = get_average_shear_normal_stress_and_average_mu_constant_speed("30ms", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor)
Average_Shear_Stress_List_40ms, Average_Mu_List_40ms, NormalStressMeans_40ms = get_average_shear_normal_stress_and_average_mu_constant_speed("40ms", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor)
Average_Shear_Stress_List_50ms, Average_Mu_List_50ms, NormalStressMeans_50ms = get_average_shear_normal_stress_and_average_mu_constant_speed("50ms", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor)

Average_Mu_List_1GPa = [Average_Mu_List_1ms[0], Average_Mu_List_10ms[0], Average_Mu_List_20ms[0], Average_Mu_List_30ms[0], Average_Mu_List_40ms[0], Average_Mu_List_50ms[0]]
Average_Mu_List_2GPa = [Average_Mu_List_1ms[1], Average_Mu_List_10ms[1], Average_Mu_List_20ms[1], Average_Mu_List_30ms[1], Average_Mu_List_40ms[1], Average_Mu_List_50ms[1]]
Average_Mu_List_3GPa = [Average_Mu_List_1ms[2], Average_Mu_List_10ms[2], Average_Mu_List_20ms[2], Average_Mu_List_30ms[2], Average_Mu_List_40ms[2], Average_Mu_List_50ms[2]]
Average_Mu_List_4GPa = [Average_Mu_List_1ms[3], Average_Mu_List_10ms[3], Average_Mu_List_20ms[3], Average_Mu_List_30ms[3], Average_Mu_List_40ms[3], Average_Mu_List_50ms[3]]
Average_Mu_List_5GPa = [Average_Mu_List_1ms[4], Average_Mu_List_10ms[4], Average_Mu_List_20ms[4], Average_Mu_List_30ms[4], Average_Mu_List_40ms[4], Average_Mu_List_50ms[4]]

print(Average_Mu_List_1GPa)
print(Average_Mu_List_2GPa)
print(Average_Mu_List_3GPa)
print(Average_Mu_List_4GPa)
print(Average_Mu_List_5GPa)

plot_shear_stress_vs_normal_stress_different_sliding_speeds(Average_Shear_Stress_List_1ms, Average_Shear_Stress_List_10ms, Average_Shear_Stress_List_20ms, Average_Shear_Stress_List_30ms, Average_Shear_Stress_List_40ms, Average_Shear_Stress_List_50ms,
                                   "1ms", "10ms", "20ms", "30ms", "40ms", "50ms")

plot_variation_in_shear_stress_constant_speed('1ms', Pressures)
plot_variation_in_shear_stress_constant_speed('10ms', Pressures)
plot_variation_in_shear_stress_constant_speed('20ms', Pressures)
plot_variation_in_shear_stress_constant_speed('30ms', Pressures)
plot_variation_in_shear_stress_constant_speed('40ms', Pressures)
plot_variation_in_shear_stress_constant_speed('50ms', Pressures)

########################## Plotting ln(Rates) vs Shear Stress #####################################
x = np.array([0, 1, 2, 3, 4, 5])
params1ms = np.polyfit(Average_Shear_Stress_List_1ms, Log_Rates_1ms, 1)
params10ms = np.polyfit(Average_Shear_Stress_List_10ms, Log_Rates_10ms, 1)
params20ms = np.polyfit(Average_Shear_Stress_List_20ms, Log_Rates_20ms, 1)
params30ms = np.polyfit(Average_Shear_Stress_List_30ms, Log_Rates_30ms, 1)
params40ms = np.polyfit(Average_Shear_Stress_List_40ms, Log_Rates_40ms, 1)
params50ms = np.polyfit(Average_Shear_Stress_List_40ms, Log_Rates_50ms, 1)

RatesvsShear, RvS3 = plt.subplots()
RvS3.set_title('Log of Dissociation Rates vs Shear Stress')
RvS3.set_xlabel('Shear Stress(GPa)')
RvS3.set_ylabel('Log of Dissociation Rate (ns-1)')
#RvS3.scatter(Average_Shear_Stress_List_1ms, Log_Rates_1ms)
RvS3.scatter(Average_Shear_Stress_List_10ms, Log_Rates_10ms)
RvS3.scatter(Average_Shear_Stress_List_20ms, Log_Rates_20ms)
RvS3.scatter(Average_Shear_Stress_List_30ms, Log_Rates_30ms)
RvS3.scatter(Average_Shear_Stress_List_40ms, Log_Rates_40ms)
RvS3.scatter(Average_Shear_Stress_List_50ms, Log_Rates_50ms)
#RvS3.plot(x, params1ms[0] * x + params1ms[1], label='1ms Fitted')
RvS3.plot(x, params10ms[0] * x + params10ms[1], label='10ms Fitted')
RvS3.plot(x, params20ms[0] * x + params20ms[1], label='20ms Fitted')
RvS3.plot(x, params30ms[0] * x + params30ms[1], label='30ms Fitted')
RvS3.plot(x, params40ms[0] * x + params40ms[1], label='40ms Fitted')
RvS3.plot(x, params50ms[0] * x + params50ms[1], label='50ms Fitted')
RvS3.set_xlim(0.5, 2.25)
RvS3.set_ylim(0, 4)
RvS3.legend(loc='lower right')
plt.show()

####### Calculate Activation Volume, Using  Carlos' conversion to get in Angstrom^3 #################

actication_vol_1ms = (params1ms[0]) * (1.38065) * 400 * 1e-2
actication_vol_10ms = (params10ms[0]) * (1.38065) * 400 * 1e-2
actication_vol_20ms = (params20ms[0]) * (1.38065) * 400 * 1e-2
actication_vol_30ms = (params30ms[0]) * (1.38065) * 400 * 1e-2
actication_vol_40ms = (params40ms[0]) * (1.38065) * 400 * 1e-2
actication_vol_50ms = (params50ms[0]) * (1.38065) * 400 * 1e-2

print(actication_vol_1ms)
print(actication_vol_10ms)
print(actication_vol_20ms)
print(actication_vol_30ms)
print(actication_vol_40ms)
print(actication_vol_50ms)


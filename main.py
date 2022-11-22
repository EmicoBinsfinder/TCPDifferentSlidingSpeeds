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
from HelperFunctions import get_intact_columns_constant_pressure_and_temperature
from HelperFunctions import get_intact_columns_constant_temperature_and_speed
from HelperFunctions import get_intact_columns_constant_pressure_and_speed
from HelperFunctions import get_MATLABFIT_dissociation_rates
from HelperFunctions import get_average_shear_normal_stress_and_average_mu_constant_temperature
from HelperFunctions import get_average_shear_normal_stress_and_average_mu_constant_speed
from HelperFunctions import plot_shear_stress_vs_normal_stress_different_sliding_speeds
from mpl_toolkits import mplot3d as Axes3D
from HelperFunctions import plot_shear_stress_vs_normal_stress, plot_variation_in_mu
from HelperFunctions import get_dissociation_rates
from HelperFunctions import plot_variation_in_shear_stress_constant_speed
from HelperFunctions import get_intact_columns_constant_pressure_and_temperature
import Dissociation_Rates
from scipy.stats.distributions import t

Temperatures = ["500K", "600K", "700K"]
Pressures = ['2GPa', '3GPa', '4GPa', '5GPa']
Speeds = ['20ms', '30ms', '40ms', '50ms']

#### 10ms
Big_Dataframe_10ms_400K = get_intact_columns_constant_temperature_and_speed(Speed="10ms", Temperature="400K", Pressures=Pressures)
Big_Dataframe_10ms_500K = get_intact_columns_constant_temperature_and_speed(Speed="10ms", Temperature="500K", Pressures=Pressures)
Big_Dataframe_10ms_600K = get_intact_columns_constant_temperature_and_speed(Speed="10ms", Temperature="600K", Pressures=Pressures)
Big_Dataframe_10ms_700K = get_intact_columns_constant_temperature_and_speed(Speed="10ms", Temperature="700K", Pressures=Pressures)

Big_Dataframe_1GPa_10ms = get_intact_columns_constant_pressure_and_speed(Pressure='1GPa', Speed="10ms", Temperatures=Temperatures)
Big_Dataframe_2GPa_10ms = get_intact_columns_constant_pressure_and_speed(Pressure='2GPa', Speed="10ms", Temperatures=Temperatures)
Big_Dataframe_3GPa_10ms = get_intact_columns_constant_pressure_and_speed(Pressure='3GPa', Speed="10ms", Temperatures=Temperatures)
Big_Dataframe_4GPa_10ms = get_intact_columns_constant_pressure_and_speed(Pressure='4GPa', Speed="10ms", Temperatures=Temperatures)
Big_Dataframe_5GPa_10ms = get_intact_columns_constant_pressure_and_speed(Pressure='5GPa', Speed="10ms", Temperatures=Temperatures)

#### 20ms
Big_Dataframe_20ms_400K = get_intact_columns_constant_temperature_and_speed(Speed="20ms", Temperature="400K", Pressures=Pressures)
Big_Dataframe_20ms_500K = get_intact_columns_constant_temperature_and_speed(Speed="20ms", Temperature="500K", Pressures=Pressures)
Big_Dataframe_20ms_600K = get_intact_columns_constant_temperature_and_speed(Speed="20ms", Temperature="600K", Pressures=Pressures)
Big_Dataframe_20ms_700K = get_intact_columns_constant_temperature_and_speed(Speed="20ms", Temperature="700K", Pressures=Pressures)

Big_Dataframe_1GPa_20ms = get_intact_columns_constant_pressure_and_speed(Pressure='1GPa', Speed="20ms", Temperatures=Temperatures)
Big_Dataframe_2GPa_20ms = get_intact_columns_constant_pressure_and_speed(Pressure='2GPa', Speed="20ms", Temperatures=Temperatures)
Big_Dataframe_3GPa_20ms = get_intact_columns_constant_pressure_and_speed(Pressure='3GPa', Speed="20ms", Temperatures=Temperatures)
Big_Dataframe_4GPa_20ms = get_intact_columns_constant_pressure_and_speed(Pressure='4GPa', Speed="20ms", Temperatures=Temperatures)
Big_Dataframe_5GPa_20ms = get_intact_columns_constant_pressure_and_speed(Pressure='5GPa', Speed="20ms", Temperatures=Temperatures)

#### 30ms
Big_Dataframe_30ms_400K = get_intact_columns_constant_temperature_and_speed(Speed="30ms", Temperature="400K", Pressures=Pressures)
Big_Dataframe_30ms_500K = get_intact_columns_constant_temperature_and_speed(Speed="30ms", Temperature="500K", Pressures=Pressures)
Big_Dataframe_30ms_600K = get_intact_columns_constant_temperature_and_speed(Speed="30ms", Temperature="600K", Pressures=Pressures)
Big_Dataframe_30ms_700K = get_intact_columns_constant_temperature_and_speed(Speed="30ms", Temperature="700K", Pressures=Pressures)

Big_Dataframe_1GPa_30ms = get_intact_columns_constant_pressure_and_speed(Pressure='1GPa', Speed="30ms", Temperatures=Temperatures)
Big_Dataframe_2GPa_30ms = get_intact_columns_constant_pressure_and_speed(Pressure='2GPa', Speed="30ms", Temperatures=Temperatures)
Big_Dataframe_3GPa_30ms = get_intact_columns_constant_pressure_and_speed(Pressure='3GPa', Speed="30ms", Temperatures=Temperatures)
Big_Dataframe_4GPa_30ms = get_intact_columns_constant_pressure_and_speed(Pressure='4GPa', Speed="30ms", Temperatures=Temperatures)
Big_Dataframe_5GPa_30ms = get_intact_columns_constant_pressure_and_speed(Pressure='5GPa', Speed="30ms", Temperatures=Temperatures)
#### 40ms
Big_Dataframe_40ms_400K = get_intact_columns_constant_temperature_and_speed(Speed="40ms", Temperature="400K", Pressures=Pressures)
Big_Dataframe_40ms_500K = get_intact_columns_constant_temperature_and_speed(Speed="40ms", Temperature="500K", Pressures=Pressures)
Big_Dataframe_40ms_600K = get_intact_columns_constant_temperature_and_speed(Speed="40ms", Temperature="600K", Pressures=Pressures)
Big_Dataframe_40ms_700K = get_intact_columns_constant_temperature_and_speed(Speed="40ms", Temperature="700K", Pressures=Pressures)

Big_Dataframe_1GPa_40ms = get_intact_columns_constant_pressure_and_speed(Pressure='1GPa', Speed="40ms", Temperatures=Temperatures)
Big_Dataframe_2GPa_40ms = get_intact_columns_constant_pressure_and_speed(Pressure='2GPa', Speed="40ms", Temperatures=Temperatures)
Big_Dataframe_3GPa_40ms = get_intact_columns_constant_pressure_and_speed(Pressure='3GPa', Speed="40ms", Temperatures=Temperatures)
Big_Dataframe_4GPa_40ms = get_intact_columns_constant_pressure_and_speed(Pressure='4GPa', Speed="40ms", Temperatures=Temperatures)
Big_Dataframe_5GPa_40ms = get_intact_columns_constant_pressure_and_speed(Pressure='5GPa', Speed="40ms", Temperatures=Temperatures)

#### 50ms
Big_Dataframe_50ms_400K = get_intact_columns_constant_temperature_and_speed(Speed="50ms", Temperature="400K", Pressures=Pressures)
Big_Dataframe_50ms_500K = get_intact_columns_constant_temperature_and_speed(Speed="50ms", Temperature="500K", Pressures=Pressures)
Big_Dataframe_50ms_600K = get_intact_columns_constant_temperature_and_speed(Speed="50ms", Temperature="600K", Pressures=Pressures)
Big_Dataframe_50ms_700K = get_intact_columns_constant_temperature_and_speed(Speed="50ms", Temperature="700K", Pressures=Pressures)

Big_Dataframe_1GPa_50ms = get_intact_columns_constant_pressure_and_speed(Pressure='1GPa', Speed="50ms", Temperatures=Temperatures)
Big_Dataframe_2GPa_50ms = get_intact_columns_constant_pressure_and_speed(Pressure='2GPa', Speed="50ms", Temperatures=Temperatures)
Big_Dataframe_3GPa_50ms = get_intact_columns_constant_pressure_and_speed(Pressure='3GPa', Speed="50ms", Temperatures=Temperatures)
Big_Dataframe_4GPa_50ms = get_intact_columns_constant_pressure_and_speed(Pressure='4GPa', Speed="50ms", Temperatures=Temperatures)
Big_Dataframe_5GPa_50ms = get_intact_columns_constant_pressure_and_speed(Pressure='5GPa', Speed="50ms", Temperatures=Temperatures)

### Constant Temperature and Pressure Combinations
Big_Dataframe_400K_1GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="400K", Pressure="1GPa")
Big_Dataframe_400K_2GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="400K", Pressure="2GPa")
Big_Dataframe_400K_3GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="400K", Pressure="3GPa")
Big_Dataframe_400K_4GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="400K", Pressure="4GPa")
Big_Dataframe_400K_5GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="400K", Pressure="5GPa")
Big_Dataframe_500K_1GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="500K", Pressure="1GPa")
Big_Dataframe_500K_2GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="500K", Pressure="2GPa")
Big_Dataframe_500K_3GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="500K", Pressure="3GPa")
Big_Dataframe_500K_4GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="500K", Pressure="4GPa")
Big_Dataframe_500K_5GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="500K", Pressure="5GPa")
Big_Dataframe_600K_1GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="600K", Pressure="1GPa")
Big_Dataframe_600K_2GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="600K", Pressure="2GPa")
Big_Dataframe_600K_3GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="600K", Pressure="3GPa")
Big_Dataframe_600K_4GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="600K", Pressure="4GPa")
Big_Dataframe_600K_5GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="600K", Pressure="5GPa")
Big_Dataframe_700K_1GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="700K", Pressure="1GPa")
Big_Dataframe_700K_2GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="700K", Pressure="2GPa")
Big_Dataframe_700K_3GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="700K", Pressure="3GPa")
Big_Dataframe_700K_4GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="700K", Pressure="4GPa")
Big_Dataframe_700K_5GPa = get_intact_columns_constant_pressure_and_temperature(Speeds=Speeds, Temperature="700K", Pressure="5GPa")

#### Saving All the Dataframes ####

# Big_Dataframe_10ms_400K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/400K10msComparison.csv', index=False)
# Big_Dataframe_10ms_500K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/500K10msComparison.csv', index=False)
# Big_Dataframe_10ms_600K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/600K10msComparison.csv', index=False)
# Big_Dataframe_10ms_700K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/700K10msComparison.csv', index=False)
# Big_Dataframe_1GPa_10ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/1GPa10msComparison.csv', index=False)
# Big_Dataframe_2GPa_10ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/2GPa10msComparison.csv', index=False)
# Big_Dataframe_3GPa_10ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/3GPa10msComparison.csv', index=False)
# Big_Dataframe_4GPa_10ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/4GPa10msComparison.csv', index=False)
# Big_Dataframe_5GPa_10ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/5GPa10msComparison.csv', index=False)
#
# Big_Dataframe_20ms_400K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/400K20msComparison.csv', index=False)
# Big_Dataframe_20ms_500K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/500K20msComparison.csv', index=False)
# Big_Dataframe_20ms_600K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/600K20msComparison.csv', index=False)
# Big_Dataframe_20ms_700K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/700K20msComparison.csv', index=False)
# Big_Dataframe_1GPa_20ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/1GPa20msComparison.csv', index=False)
# Big_Dataframe_2GPa_20ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/2GPa20msComparison.csv', index=False)
# Big_Dataframe_3GPa_20ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/3GPa20msComparison.csv', index=False)
# Big_Dataframe_4GPa_20ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/4GPa20msComparison.csv', index=False)
# Big_Dataframe_5GPa_20ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/5GPa20msComparison.csv', index=False)
#
# Big_Dataframe_30ms_400K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/400K30msComparison.csv', index=False)
# Big_Dataframe_30ms_500K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/500K30msComparison.csv', index=False)
# Big_Dataframe_30ms_600K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/600K30msComparison.csv', index=False)
# Big_Dataframe_30ms_700K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/700K30msComparison.csv', index=False)
# Big_Dataframe_1GPa_30ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/1GPa30msComparison.csv', index=False)
# Big_Dataframe_2GPa_30ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/2GPa30msComparison.csv', index=False)
# Big_Dataframe_3GPa_30ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/3GPa30msComparison.csv', index=False)
# Big_Dataframe_4GPa_30ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/4GPa30msComparison.csv', index=False)
# Big_Dataframe_5GPa_30ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/5GPa30msComparison.csv', index=False)
#
# Big_Dataframe_40ms_400K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/400K40msComparison.csv', index=False)
# Big_Dataframe_40ms_500K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/500K40msComparison.csv', index=False)
# Big_Dataframe_40ms_600K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/600K40msComparison.csv', index=False)
# Big_Dataframe_40ms_700K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/700K40msComparison.csv', index=False)
# Big_Dataframe_1GPa_40ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/1GPa40msComparison.csv', index=False)
# Big_Dataframe_2GPa_40ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/2GPa40msComparison.csv', index=False)
# Big_Dataframe_3GPa_40ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/3GPa40msComparison.csv', index=False)
# Big_Dataframe_4GPa_40ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/4GPa40msComparison.csv', index=False)
# Big_Dataframe_5GPa_40ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/5GPa40msComparison.csv', index=False)
#
# Big_Dataframe_50ms_400K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/400K50msComparison.csv', index=False)
# Big_Dataframe_50ms_500K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/500K50msComparison.csv', index=False)
# Big_Dataframe_50ms_600K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/600K50msComparison.csv', index=False)
# Big_Dataframe_50ms_700K.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/700K50msComparison.csv', index=False)
# Big_Dataframe_1GPa_50ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/1GPa50msComparison.csv', index=False)
# Big_Dataframe_2GPa_50ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/2GPa50msComparison.csv', index=False)
# Big_Dataframe_3GPa_50ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/3GPa50msComparison.csv', index=False)
# Big_Dataframe_4GPa_50ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/4GPa50msComparison.csv', index=False)
# Big_Dataframe_5GPa_50ms.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/5GPa50msComparison.csv', index=False)
#
# Big_Dataframe_400K_1GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/400K1GPaComparison.csv', index = False)
# Big_Dataframe_400K_2GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/400K2GPaComparison.csv', index = False)
# Big_Dataframe_400K_3GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/400K3GPaComparison.csv', index = False)
# Big_Dataframe_400K_4GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/400K4GPaComparison.csv', index = False)
# Big_Dataframe_400K_5GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/400K5GPaComparison.csv', index = False)
# Big_Dataframe_500K_1GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/500K1GPaComparison.csv', index = False)
# Big_Dataframe_500K_2GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/500K2GPaComparison.csv', index = False)
# Big_Dataframe_500K_3GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/500K3GPaComparison.csv', index = False)
# Big_Dataframe_500K_4GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/500K4GPaComparison.csv', index = False)
# Big_Dataframe_500K_5GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/500K5GPaComparison.csv', index = False)
# Big_Dataframe_600K_1GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/600K1GPaComparison.csv', index = False)
# Big_Dataframe_600K_2GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/600K2GPaComparison.csv', index = False)
# Big_Dataframe_600K_3GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/600K3GPaComparison.csv', index = False)
# Big_Dataframe_600K_4GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/600K4GPaComparison.csv', index = False)
# Big_Dataframe_600K_5GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/600K5GPaComparison.csv', index = False)
# Big_Dataframe_700K_1GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/700K1GPaComparison.csv', index = False)
# Big_Dataframe_700K_2GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/700K2GPaComparison.csv', index = False)
# Big_Dataframe_700K_3GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/700K3GPaComparison.csv', index = False)
# Big_Dataframe_700K_4GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/700K4GPaComparison.csv', index = False)
# Big_Dataframe_700K_5GPa.to_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/700K5GPaComparison.csv', index = False)

Timestep = Big_Dataframe_10ms_400K['Timestep'].to_list()

# ################# Getting Dissociation Rates for Controlled Speed Using MATLAB Fits ##########
Cutoff = None
Index = 0 # Index = 0 for actual rate, 1 for upper bound, 2 for lower bound

# Dissociation_Rate_10ms_400K_1GPa, LogRate_10ms_400K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_400K_1GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_10ms_400K_2GPa, LogRate_10ms_400K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_400K_2GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_10ms_400K_3GPa, LogRate_10ms_400K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_400K_3GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_10ms_400K_4GPa, LogRate_10ms_400K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_400K_4GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_10ms_400K_5GPa, LogRate_10ms_400K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_400K_5GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rates_10ms_400K = [Dissociation_Rate_10ms_400K_1GPa, Dissociation_Rate_10ms_400K_2GPa, Dissociation_Rate_10ms_400K_3GPa, Dissociation_Rate_10ms_400K_4GPa, Dissociation_Rate_10ms_400K_5GPa]
# Log_Rates_400K_10ms = [LogRate_10ms_400K_1GPa, LogRate_10ms_400K_2GPa, LogRate_10ms_400K_3GPa, LogRate_10ms_400K_4GPa, LogRate_10ms_400K_5GPa]
#
# Dissociation_Rate_10ms_500K_1GPa, LogRate_10ms_500K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_500K_1GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_10ms_500K_2GPa, LogRate_10ms_500K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_500K_2GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_10ms_500K_3GPa, LogRate_10ms_500K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_500K_3GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_10ms_500K_4GPa, LogRate_10ms_500K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_500K_4GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_10ms_500K_5GPa, LogRate_10ms_500K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_500K_5GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rates_10ms_500K = [Dissociation_Rate_10ms_500K_1GPa, Dissociation_Rate_10ms_500K_2GPa, Dissociation_Rate_10ms_500K_3GPa, Dissociation_Rate_10ms_500K_4GPa, Dissociation_Rate_10ms_500K_5GPa]
# Log_Rates_500K_10ms = [LogRate_10ms_500K_1GPa, LogRate_10ms_500K_2GPa, LogRate_10ms_500K_3GPa, LogRate_10ms_500K_4GPa, LogRate_10ms_500K_5GPa]
#
# Dissociation_Rate_10ms_600K_1GPa, LogRate_10ms_600K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_600K_1GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_10ms_600K_2GPa, LogRate_10ms_600K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_600K_2GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_10ms_600K_3GPa, LogRate_10ms_600K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_600K_3GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_10ms_600K_4GPa, LogRate_10ms_600K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_600K_4GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_10ms_600K_5GPa, LogRate_10ms_600K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_600K_5GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rates_10ms_600K = [Dissociation_Rate_10ms_600K_1GPa, Dissociation_Rate_10ms_600K_2GPa, Dissociation_Rate_10ms_600K_3GPa, Dissociation_Rate_10ms_600K_4GPa, Dissociation_Rate_10ms_600K_5GPa]
# Log_Rates_600K_10ms = [LogRate_10ms_600K_1GPa, LogRate_10ms_600K_2GPa, LogRate_10ms_600K_3GPa, LogRate_10ms_600K_4GPa, LogRate_10ms_600K_5GPa]
#
# Dissociation_Rate_10ms_700K_1GPa, LogRate_10ms_700K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_700K_1GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_10ms_700K_2GPa, LogRate_10ms_700K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_700K_2GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_10ms_700K_3GPa, LogRate_10ms_700K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_700K_3GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_10ms_700K_4GPa, LogRate_10ms_700K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_700K_4GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_10ms_700K_5GPa, LogRate_10ms_700K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_10ms_700K_5GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rates_10ms_700K = [Dissociation_Rate_10ms_700K_1GPa, Dissociation_Rate_10ms_700K_2GPa, Dissociation_Rate_10ms_700K_3GPa, Dissociation_Rate_10ms_700K_4GPa, Dissociation_Rate_10ms_700K_5GPa]
# Log_Rates_700K_10ms = [LogRate_10ms_700K_1GPa, LogRate_10ms_700K_2GPa, LogRate_10ms_700K_3GPa, LogRate_10ms_700K_4GPa, LogRate_10ms_700K_5GPa]
#
# Dissociation_Rates_1GPa = [Dissociation_Rates_10ms_400K[0], Dissociation_Rates_10ms_500K[0], Dissociation_Rates_10ms_600K[0], Dissociation_Rates_10ms_700K[0]]
# Dissociation_Rates_2GPa = [Dissociation_Rates_10ms_400K[1], Dissociation_Rates_10ms_500K[1], Dissociation_Rates_10ms_600K[1], Dissociation_Rates_10ms_700K[1]]
# Dissociation_Rates_3GPa = [Dissociation_Rates_10ms_400K[2], Dissociation_Rates_10ms_500K[2], Dissociation_Rates_10ms_600K[2], Dissociation_Rates_10ms_700K[2]]
# Dissociation_Rates_4GPa = [Dissociation_Rates_10ms_400K[3], Dissociation_Rates_10ms_500K[3], Dissociation_Rates_10ms_600K[3], Dissociation_Rates_10ms_700K[3]]
# Dissociation_Rates_5GPa = [Dissociation_Rates_10ms_400K[4], Dissociation_Rates_10ms_500K[4], Dissociation_Rates_10ms_600K[4], Dissociation_Rates_10ms_700K[4]]
#
# Log_Rates_1GPa_10ms = [LogRate_10ms_400K_1GPa, LogRate_10ms_500K_1GPa, LogRate_10ms_600K_1GPa, LogRate_10ms_700K_1GPa]
# Log_Rates_2GPa_10ms = [LogRate_10ms_400K_2GPa, LogRate_10ms_500K_2GPa, LogRate_10ms_600K_2GPa, LogRate_10ms_700K_2GPa]
# Log_Rates_3GPa_10ms = [LogRate_10ms_400K_3GPa, LogRate_10ms_500K_3GPa, LogRate_10ms_600K_3GPa, LogRate_10ms_700K_3GPa]
# Log_Rates_4GPa_10ms = [LogRate_10ms_400K_4GPa, LogRate_10ms_500K_4GPa, LogRate_10ms_600K_4GPa, LogRate_10ms_700K_4GPa]
# Log_Rates_5GPa_10ms = [LogRate_10ms_400K_5GPa, LogRate_10ms_500K_5GPa, LogRate_10ms_600K_5GPa, LogRate_10ms_700K_5GPa]
#
# EquilibriumFactor = [79, 99] # How many rows (out of 99) to ignore before calculating shear stress/friction coefficient, as it won't stabilise until after a certain number of timesteps
#
# def get_average_shear_normal_stress_and_average_mu_constant_temperature(Temperature, Pressures, EquilibriumFactor, Speed):
#     Friction_Coefficient_Dataframe_Unnamed = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{Speed}/{Temperature}/1GPa/'
#                                 'fc_ave.dump'.format(Speed=Speed, Temperature=Temperature), sep=' ')
#     Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe_Unnamed.rename(columns={'v_s_bot' : 'Shear Stress 1GPa', 'v_p_bot' : 'Normal Stress 1GPa'})
#
#     for P in Pressures:
#         Dataframe = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{Speed}/{Temperature}/{P}/'
#                                 'fc_ave.dump'.format(Speed=Speed, Temperature=Temperature, P=P), sep=' ')
#         Big_DataframeP = Dataframe.rename(columns= {'Timestep': 'Timestep {}'.format(P),
#                                                         'v_s_bot': 'Shear Stress {}'.format(P),
#                                                         'v_p_bot': 'Normal Stress {}'.format(P)})
#
#         Friction_Coefficient_Dataframe = pd.concat([Friction_Coefficient_Dataframe, Big_DataframeP], axis =1)
#         Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe.dropna()
#
#
#     #print(Friction_Coefficient_Dataframe)
#     Mu_Final_Dataframe = Friction_Coefficient_Dataframe.iloc[:, [0, 1, 2, 4, 5, 7, 8, 10, 11, 13, 14]]
#     Mu_Final_Dataframe = Mu_Final_Dataframe.iloc[EquilibriumFactor[0]:EquilibriumFactor[1], :]
#     #print(Mu_Final_Dataframe)
#
#     ShearStressMeans = Mu_Final_Dataframe[['Shear Stress 1GPa', 'Shear Stress 2GPa', 'Shear Stress 3GPa', 'Shear Stress 4GPa', 'Shear Stress 5GPa']].mean()
#     Average_Shear_Stress_Dictionary = ShearStressMeans.to_dict()
#     #print(ShearStressMeans)
#     NormalStressMeans = Mu_Final_Dataframe[['Normal Stress 1GPa', 'Normal Stress 2GPa', 'Normal Stress 3GPa', 'Normal Stress 4GPa', 'Normal Stress 5GPa']].mean()
#     NormalStressMeans = NormalStressMeans.to_dict()
#     #print(NormalStressMeans)
#
#     Average_Mu_Dictionary = {}
#
#     NormalStressMeansList = list(NormalStressMeans.values())
#     Normal_Stress = NormalStressMeans.get('Normal Stress 1GPa')
#     Shear_Stress = ShearStressMeans.get('Shear Stress 1GPa')
#     Average_Mu = Shear_Stress / Normal_Stress
#     Average_Mu_Dictionary.update({'Average Mu 1GPa': Average_Mu})
#
#     for P in Pressures:
#
#         Normal_Stress = NormalStressMeans.get('Normal Stress {}'.format(P))
#         Shear_Stress = ShearStressMeans.get('Shear Stress {}'.format(P))
#         Average_Mu = Shear_Stress / Normal_Stress
#         Average_Mu_Dictionary.update({'Average Mu {}'.format(P): Average_Mu})
#
#
#     Average_Shear_Stress_List = list(Average_Shear_Stress_Dictionary.values())
#     #print(Average_Shear_Stress_List)
#     Average_Mu_List = list(Average_Mu_Dictionary.values())
#     Average_Shear_Stress_List = [x / 10000 for x in Average_Shear_Stress_List] # Conversion to GPa
#     NormalStressMeansList = [x/10000 for x in NormalStressMeansList]
#     #print(Average_Shear_Stress_List)
#
#     return Average_Shear_Stress_List, Average_Mu_List, NormalStressMeansList
#
# ######## Getting Average Shear Stress, Friction Coefficient and Normal Stress #################
# Average_Shear_Stress_List_400K_10ms, Average_Mu_List_400K_10ms, NormalStressMeans_400K_10ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("400K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="10ms")
# Average_Shear_Stress_List_500K_10ms, Average_Mu_List_500K_10ms, NormalStressMeans_500K_10ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("500K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="10ms")
# Average_Shear_Stress_List_600K_10ms, Average_Mu_List_600K_10ms, NormalStressMeans_600K_10ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("600K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="10ms")
# Average_Shear_Stress_List_700K_10ms, Average_Mu_List_700K_10ms, NormalStressMeans_700K_10ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("700K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="10ms")
#
# Average_Mu_List_1GPa_10ms = [Average_Mu_List_400K_10ms[0], Average_Mu_List_500K_10ms[0], Average_Mu_List_600K_10ms[0], Average_Mu_List_700K_10ms[0]]
# Average_Mu_List_2GPa_10ms = [Average_Mu_List_400K_10ms[1], Average_Mu_List_500K_10ms[1], Average_Mu_List_600K_10ms[1], Average_Mu_List_700K_10ms[1]]
# Average_Mu_List_3GPa_10ms = [Average_Mu_List_400K_10ms[2], Average_Mu_List_500K_10ms[2], Average_Mu_List_600K_10ms[2], Average_Mu_List_700K_10ms[2]]
# Average_Mu_List_4GPa_10ms = [Average_Mu_List_400K_10ms[3], Average_Mu_List_500K_10ms[3], Average_Mu_List_600K_10ms[3], Average_Mu_List_700K_10ms[3]]
# Average_Mu_List_5GPa_10ms = [Average_Mu_List_400K_10ms[4], Average_Mu_List_500K_10ms[4], Average_Mu_List_600K_10ms[4], Average_Mu_List_700K_10ms[4]]
#
# Average_Shear_Stress_List_1GPa_10ms = [Average_Shear_Stress_List_400K_10ms[0], Average_Shear_Stress_List_500K_10ms[0], Average_Shear_Stress_List_600K_10ms[0], Average_Shear_Stress_List_700K_10ms[0]]
# Average_Shear_Stress_List_2GPa_10ms = [Average_Shear_Stress_List_400K_10ms[1], Average_Shear_Stress_List_500K_10ms[1], Average_Shear_Stress_List_600K_10ms[1], Average_Shear_Stress_List_700K_10ms[1]]
# Average_Shear_Stress_List_3GPa_10ms = [Average_Shear_Stress_List_400K_10ms[2], Average_Shear_Stress_List_500K_10ms[2], Average_Shear_Stress_List_600K_10ms[2], Average_Shear_Stress_List_700K_10ms[2]]
# Average_Shear_Stress_List_4GPa_10ms = [Average_Shear_Stress_List_400K_10ms[3], Average_Shear_Stress_List_500K_10ms[3], Average_Shear_Stress_List_600K_10ms[3], Average_Shear_Stress_List_700K_10ms[3]]
# Average_Shear_Stress_List_5GPa_10ms = [Average_Shear_Stress_List_400K_10ms[4], Average_Shear_Stress_List_500K_10ms[4], Average_Shear_Stress_List_600K_10ms[4], Average_Shear_Stress_List_700K_10ms[4]]
#
# plot_shear_stress_vs_normal_stress(Average_Shear_Stress_List_400K_10ms, Average_Shear_Stress_List_500K_10ms, Average_Shear_Stress_List_600K_10ms, Average_Shear_Stress_List_700K_10ms,
#                                  "400K", "500K", "600K", "700K", Speed="10ms")
#
# ########################## Plotting ln(Rates) vs Shear Stress #####################################
# x = np.array([0, 1, 2, 3, 4, 5])
# params400K_10ms = np.polyfit(Average_Shear_Stress_List_400K_10ms, Log_Rates_400K_10ms, 1)
# params500K_10ms = np.polyfit(Average_Shear_Stress_List_500K_10ms, Log_Rates_500K_10ms, 1)
# params600K_10ms = np.polyfit(Average_Shear_Stress_List_600K_10ms, Log_Rates_600K_10ms, 1)
# params700K_10ms = np.polyfit(Average_Shear_Stress_List_700K_10ms, Log_Rates_700K_10ms, 1)
#
# RatesvsShear, RvS3 = plt.subplots()
# RvS3.set_title('Log of Dissociation Rates vs Shear Stress, 10ms')
# RvS3.set_xlabel('Shear Stress(GPa)')
# RvS3.set_ylabel('Log of Dissociation Rate (per nanosecond)')
# RvS3.scatter(Average_Shear_Stress_List_400K_10ms, Log_Rates_400K_10ms)
# RvS3.scatter(Average_Shear_Stress_List_500K_10ms, Log_Rates_500K_10ms)
# RvS3.scatter(Average_Shear_Stress_List_600K_10ms, Log_Rates_600K_10ms)
# RvS3.scatter(Average_Shear_Stress_List_700K_10ms, Log_Rates_700K_10ms)
# RvS3.plot(x, params400K_10ms[0] * x + params400K_10ms[1], label='400K Fitted')
# RvS3.plot(x, params500K_10ms[0] * x + params500K_10ms[1], label='500K Fitted')
# RvS3.plot(x, params600K_10ms[0] * x + params600K_10ms[1], label='600K Fitted')
# RvS3.plot(x, params700K_10ms[0] * x + params700K_10ms[1], label='700K Fitted')
# RvS3.set_xlim(0, 2)
# RvS3.set_ylim(-0.5, 4)
# RvS3.legend(loc='lower right')
# plt.show()
#
# ####### Calculate Activation Volume, Using  Carlos' conversion to get in Angstrom^3 #################
#
# activation_vol_400K_10ms = (params400K_10ms[0]) * (1.38065) * 400 * 1e-2
# activation_vol_500K_10ms = (params500K_10ms[0]) * (1.38065) * 500 * 1e-2
# activation_vol_600K_10ms = (params600K_10ms[0]) * (1.38065) * 600 * 1e-2
# activation_vol_700K_10ms = (params700K_10ms[0]) * (1.38065) * 700 * 1e-2
#
# alpha = 0.05
# def linear(x, m, n):
#     return m * x + n
#
# coef_sh_cf, coef_sh_pcov = optimize.curve_fit(linear, Average_Shear_Stress_List_700K_10ms, Log_Rates_700K_10ms)
# sigma = coef_sh_pcov[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_700K_10ms) - len(params700K_10ms))
#
# tval = t.ppf(1.0 - alpha / 2., dof)
# uncert = sigma * tval
# uncert400_10ms = uncert * (1.38065) * 400 * 1e-2
# uncert500_10ms = uncert * (1.38065) * 500 * 1e-2
# uncert600_10ms = uncert * (1.38065) * 600 * 1e-2
# uncert700_10ms = uncert * (1.38065) * 700 * 1e-2
#
# print(f'Activation Volume uncertainty at 400K, 10ms is {uncert400_10ms}')
# print(f'Activation Volume uncertainty at 500K, 10ms is {uncert500_10ms}')
# print(f'Activation Volume uncertainty at 600K, 10ms is {uncert600_10ms}')
# print(f'Activation Volume uncertainty at 700K, 10ms is {uncert700_10ms}')
#
# Activation_Volumes_10ms = [activation_vol_400K_10ms, activation_vol_500K_10ms, activation_vol_600K_10ms, activation_vol_700K_10ms]
# mean_actv_10ms = np.average(Activation_Volumes_10ms)
#
# print('Activation Volume 400K, 10ms = ' + str(activation_vol_400K_10ms))
# print('Activation Volume 500K, 10ms = ' + str(activation_vol_500K_10ms))
# print('Activation Volume 600K, 10ms = ' + str(activation_vol_600K_10ms))
# print('Activation Volume 700K, 10ms = ' + str(activation_vol_700K_10ms))
# #
# ############ Plotting lnk vs 1000/T  #########################
#
# Temperatures = [400, 500, 600, 700]
# Inverse_Temperatures = np.array([1/x for x in Temperatures])
#
# trend1GPa_10ms = np.polyfit(Inverse_Temperatures, Log_Rates_1GPa_10ms, 1)
# trend2GPa_10ms = np.polyfit(Inverse_Temperatures, Log_Rates_2GPa_10ms, 1)
# trend3GPa_10ms = np.polyfit(Inverse_Temperatures, Log_Rates_3GPa_10ms, 1)
# trend4GPa_10ms = np.polyfit(Inverse_Temperatures, Log_Rates_4GPa_10ms, 1)
# trend5GPa_10ms = np.polyfit(Inverse_Temperatures, Log_Rates_5GPa_10ms, 1)
#
# fig1, ax1 = plt.subplots()
# ax1.set_title('Log of Dissociation Rates against Inverse of Temperatures, 10ms')
# ax1.set_xlabel('1000/T (K-1)')
# ax1.set_ylabel('ln(Rate) (ns-1)')
# ax1.scatter(Inverse_Temperatures, Log_Rates_1GPa_10ms)
# ax1.scatter(Inverse_Temperatures, Log_Rates_2GPa_10ms)
# ax1.scatter(Inverse_Temperatures, Log_Rates_3GPa_10ms)
# ax1.scatter(Inverse_Temperatures, Log_Rates_4GPa_10ms)
# ax1.scatter(Inverse_Temperatures, Log_Rates_5GPa_10ms)
#
# Fit1GPa_10ms = np.poly1d(trend1GPa_10ms)
# Fit2GPa_10ms = np.poly1d(trend2GPa_10ms)
# Fit3GPa_10ms = np.poly1d(trend3GPa_10ms)
# Fit4GPa_10ms = np.poly1d(trend4GPa_10ms)
# Fit5GPa_10ms = np.poly1d(trend5GPa_10ms)
#
# ax1.plot(Inverse_Temperatures, Fit1GPa_10ms(Inverse_Temperatures), label='1GPa')
# ax1.plot(Inverse_Temperatures, Fit2GPa_10ms(Inverse_Temperatures), label='2GPa')
# ax1.plot(Inverse_Temperatures, Fit3GPa_10ms(Inverse_Temperatures), label='3GPa')
# ax1.plot(Inverse_Temperatures, Fit4GPa_10ms(Inverse_Temperatures), label='4GPa')
# ax1.plot(Inverse_Temperatures, Fit5GPa_10ms(Inverse_Temperatures), label='5GPa')
# ax1.legend()
# plt.show()
#
# ActivationEnergy_1GPa_10ms = (((1 * 1e9 * (np.average(Average_Mu_List_1GPa_10ms) * mean_actv_10ms) * 1e-30) - 1.381 * trend1GPa_10ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
# ActivationEnergy_2GPa_10ms = (((2 * 1e9 * (np.average(Average_Mu_List_2GPa_10ms) * mean_actv_10ms) * 1e-30) - 1.381 * trend2GPa_10ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
# ActivationEnergy_3GPa_10ms = (((3 * 1e9 * (np.average(Average_Mu_List_3GPa_10ms) * mean_actv_10ms) * 1e-30) - 1.381 * trend3GPa_10ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
# ActivationEnergy_4GPa_10ms = (((4 * 1e9 * (np.average(Average_Mu_List_4GPa_10ms) * mean_actv_10ms) * 1e-30) - 1.381 * trend4GPa_10ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
# ActivationEnergy_5GPa_10ms = (((5 * 1e9 * (np.average(Average_Mu_List_5GPa_10ms) * mean_actv_10ms) * 1e-30) - 1.381 * trend5GPa_10ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
#
# print(np.average(Average_Mu_List_1GPa_10ms))
# print(np.average(Average_Mu_List_2GPa_10ms))
# print(np.average(Average_Mu_List_3GPa_10ms))
# print(np.average(Average_Mu_List_4GPa_10ms))
# print(np.average(Average_Mu_List_5GPa_10ms))
#
# lnA1GPa_10ms = np.log((np.exp(trend1GPa_10ms[1]) * (10 ** 9)))
# lnA2GPa_10ms = np.log((np.exp(trend2GPa_10ms[1]) * (10 ** 9)))
# lnA3GPa_10ms = np.log((np.exp(trend3GPa_10ms[1]) * (10 ** 9)))
# lnA4GPa_10ms = np.log((np.exp(trend4GPa_10ms[1]) * (10 ** 9)))
# lnA5GPa_10ms = np.log((np.exp(trend5GPa_10ms[1]) * (10 ** 9)))
#
# print(f"ln(A) at 1GPa, 10ms is {lnA1GPa_10ms}")
# print(f"ln(A) at 2GPa, 10ms is {lnA2GPa_10ms}")
# print(f"ln(A) at 3GPa, 10ms is {lnA3GPa_10ms}")
# print(f"ln(A) at 4GPa, 10ms is {lnA4GPa_10ms}")
# print(f"ln(A) at 5GPa, 10ms is {lnA5GPa_10ms}")
#
# print('Activation Energy 1GPa, 10ms =' + str(ActivationEnergy_1GPa_10ms))
# print('Activation Energy 2GPa, 10ms =' + str(ActivationEnergy_2GPa_10ms))
# print('Activation Energy 3GPa, 10ms =' + str(ActivationEnergy_3GPa_10ms))
# print('Activation Energy 4GPa, 10ms =' + str(ActivationEnergy_4GPa_10ms))
# print('Activation Energy 5GPa, 10ms =' + str(ActivationEnergy_5GPa_10ms))
#
# def linear(x, m, n):
#     return m * x + n
#
# coef_p_cf_1GPa_10ms, coef_p_pcov_1GPa_10ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_1GPa_10ms)
# sigma_A1GPa_10ms = coef_p_pcov_1GPa_10ms[1, 1] ** 0.5
# sigma_m1GPa_10ms = coef_p_pcov_1GPa_10ms[0, 0] ** 0.5
# dof_10ms = max(0, len(Log_Rates_1GPa_10ms) - len(coef_p_cf_1GPa_10ms))
#
# alpha = 0.05
# tval_10ms = t.ppf(1.0 - alpha / 2., 4.1532)
#
# sigma_10ms = np.std(Activation_Volumes_10ms)
# error_actv_10ms = sigma_10ms * tval_10ms / np.sqrt(len(Activation_Volumes_10ms))
# uncert_A1GPa_10ms = sigma_A1GPa_10ms * tval
# uncert_m1GPa_10ms = sigma_m1GPa_10ms * tval
#
# coef_p_cf_2GPa_10ms, coef_p_pcov_2GPa_10ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_2GPa_10ms)
# sigma_A2GPa_10ms = coef_p_pcov_2GPa_10ms[1, 1] ** 0.5
# sigma_m2GPa_10ms = coef_p_pcov_2GPa_10ms[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_2GPa_10ms) - len(coef_p_cf_2GPa_10ms))
# uncert_A2GPa_10ms = sigma_A2GPa_10ms * tval_10ms
# uncert_m2GPa_10ms = sigma_m2GPa_10ms * tval_10ms
#
# coef_p_cf_3GPa_10ms, coef_p_pcov_3GPa_10ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_3GPa_10ms)
# sigma_A3GPa_10ms = coef_p_pcov_3GPa_10ms[1, 1] ** 0.5
# sigma_m3GPa_10ms = coef_p_pcov_3GPa_10ms[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_3GPa_10ms) - len(coef_p_cf_3GPa_10ms))
# uncert_A3GPa_10ms = sigma_A3GPa_10ms * tval_10ms
# uncert_m3GPa_10ms = sigma_m3GPa_10ms * tval_10ms
#
# coef_p_cf_4GPa_10ms, coef_p_pcov_4GPa_10ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_4GPa_10ms)
# sigma_A4GPa_10ms = coef_p_pcov_4GPa_10ms[1, 1] ** 0.5
# sigma_m4GPa_10ms = coef_p_pcov_4GPa_10ms[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_4GPa_10ms) - len(coef_p_cf_4GPa_10ms))
# uncert_A4GPa_10ms = sigma_A4GPa_10ms * tval_10ms
# uncert_m4GPa_10ms = sigma_m4GPa_10ms * tval_10ms
#
# coef_p_cf_5GPa_10ms, coef_p_pcov_5GPa_10ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_5GPa_10ms)
# sigma_A5GPa_10ms = coef_p_pcov_5GPa_10ms[1, 1] ** 0.5
# sigma_m5GPa_10ms = coef_p_pcov_5GPa_10ms[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_5GPa_10ms) - len(coef_p_cf_5GPa_10ms))
# uncert_A5GPa_10ms = sigma_A5GPa_10ms * tval_10ms
# uncert_m5GPa_10ms = sigma_m5GPa_10ms * tval_10ms
#
# err_E_1GPa_10ms = (np.sqrt(((1 * 1e9 * np.average(Average_Mu_List_1GPa_10ms) * error_actv_10ms * 1e-30) ** 2 + (1.381 * uncert_m1GPa_10ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
# err_E_2GPa_10ms = (np.sqrt(((2 * 1e9 * np.average(Average_Mu_List_2GPa_10ms) * error_actv_10ms * 1e-30) ** 2 + (1.381 * uncert_m2GPa_10ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
# err_E_3GPa_10ms = (np.sqrt(((3 * 1e9 * np.average(Average_Mu_List_3GPa_10ms) * error_actv_10ms * 1e-30) ** 2 + (1.381 * uncert_m3GPa_10ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
# err_E_4GPa_10ms = (np.sqrt(((4 * 1e9 * np.average(Average_Mu_List_4GPa_10ms) * error_actv_10ms * 1e-30) ** 2 + (1.381 * uncert_m4GPa_10ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
# err_E_5GPa_10ms = (np.sqrt(((5 * 1e9 * np.average(Average_Mu_List_5GPa_10ms) * error_actv_10ms * 1e-30) ** 2 + (1.381 * uncert_m5GPa_10ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
#
# print(f"Error for Activation Energy at 1GPa = {err_E_1GPa_10ms}")
# print(f"Error for Activation Energy at 2GPa = {err_E_2GPa_10ms}")
# print(f"Error for Activation Energy at 3GPa = {err_E_3GPa_10ms}")
# print(f"Error for Activation Energy at 4GPa = {err_E_4GPa_10ms}")
# print(f"Error for Activation Energy at 5GPa = {err_E_5GPa_10ms}")
#
# uncertlnA1GPa_10ms = np.log(uncert_A1GPa_10ms * (10 ** 9))
# uncertlnA2GPa_10ms = np.log(uncert_A2GPa_10ms * (10 ** 9))
# uncertlnA3GPa_10ms = np.log(uncert_A3GPa_10ms * (10 ** 9))
# uncertlnA4GPa_10ms = np.log(uncert_A4GPa_10ms * (10 ** 9))
# uncertlnA5GPa_10ms = np.log(uncert_A5GPa_10ms * (10 ** 9))
#
# print(f"Error for ln A at 1GPa = {uncert_A1GPa_10ms}")
# print(f"Error for ln A at 2GPa = {uncert_A2GPa_10ms}")
# print(f"Error for ln A at 3GPa = {uncert_A3GPa_10ms}")
# print(f"Error for ln A at 4GPa = {uncert_A4GPa_10ms}")
# print(f"Error for ln A at 5GPa = {uncert_A5GPa_10ms}")
#
# pressures = ['2GPa', '3GPa', '4GPa', '5GPa']
#
# ########################## Plotting ln(Rates) vs Normal Stress #####################################
# x = np.array([0, 1, 2, 3, 4, 5])
# params400K_10ms = np.polyfit(NormalStressMeans_400K_10ms, Log_Rates_400K_10ms, 1)
# params500K_10ms = np.polyfit(NormalStressMeans_500K_10ms, Log_Rates_500K_10ms, 1)
# params600K_10ms = np.polyfit(NormalStressMeans_600K_10ms, Log_Rates_600K_10ms, 1)
# params700K_10ms = np.polyfit(NormalStressMeans_700K_10ms, Log_Rates_700K_10ms, 1)
#
# RatesvsNormal, RvN3 = plt.subplots()
# RvN3.set_title('Log of Dissociation Rates vs Normal Stress - Alpha Fe, 10ms')
# RvN3.set_xlabel('Normal Stress(GPa)')
# RvN3.set_ylabel('Log of Dissociation Rate (ns-1)')
# RvN3.scatter(NormalStressMeans_400K_10ms, Log_Rates_400K_10ms)
# RvN3.scatter(NormalStressMeans_500K_10ms, Log_Rates_500K_10ms)
# RvN3.scatter(NormalStressMeans_600K_10ms, Log_Rates_600K_10ms)
# RvN3.scatter(NormalStressMeans_700K_10ms, Log_Rates_700K_10ms)
# RvN3.plot(x, params400K_10ms[0] * x + params400K_10ms[1], label='400K Fitted')
# RvN3.plot(x, params500K_10ms[0] * x + params500K_10ms[1], label='500K Fitted')
# RvN3.plot(x, params600K_10ms[0] * x + params600K_10ms[1], label='600K Fitted')
# RvN3.plot(x, params700K_10ms[0] * x + params700K_10ms[1], label='700K Fitted')
# #RvN3.set_xlim(1, 5)
# # RvN3.set_ylim(0, 25)
# RvN3.legend(loc='lower right')
# plt.show()
#
# def function(data, a, b, c):
#     x = data[0]
#     y = data[1]
#     return a * (x**b) * (y**c) #TODO change fitting function
#
# x_data = []
# y_data = []
# z_data = []
#
# data = [[400, Average_Shear_Stress_List_400K_10ms[0], LogRate_10ms_400K_1GPa], [400, Average_Shear_Stress_List_400K_10ms[1], LogRate_10ms_400K_2GPa], [400, Average_Shear_Stress_List_400K_10ms[2], LogRate_10ms_400K_3GPa], [400, Average_Shear_Stress_List_400K_10ms[3], LogRate_10ms_400K_4GPa], [400, Average_Shear_Stress_List_400K_10ms[4], LogRate_10ms_400K_5GPa],
#         [500, Average_Shear_Stress_List_500K_10ms[0], LogRate_10ms_500K_1GPa], [500, Average_Shear_Stress_List_500K_10ms[1], LogRate_10ms_500K_2GPa], [500, Average_Shear_Stress_List_500K_10ms[2], LogRate_10ms_500K_3GPa], [500, Average_Shear_Stress_List_500K_10ms[3], LogRate_10ms_500K_4GPa], [500, Average_Shear_Stress_List_500K_10ms[4], LogRate_10ms_500K_5GPa],
#         [600, Average_Shear_Stress_List_600K_10ms[0], LogRate_10ms_600K_1GPa], [600, Average_Shear_Stress_List_600K_10ms[1], LogRate_10ms_600K_2GPa], [600, Average_Shear_Stress_List_600K_10ms[2], LogRate_10ms_600K_3GPa], [600, Average_Shear_Stress_List_600K_10ms[3], LogRate_10ms_600K_4GPa], [600, Average_Shear_Stress_List_600K_10ms[4], LogRate_10ms_600K_5GPa],
#         [700, Average_Shear_Stress_List_700K_10ms[0], LogRate_10ms_700K_1GPa], [700, Average_Shear_Stress_List_700K_10ms[1], LogRate_10ms_700K_2GPa], [700, Average_Shear_Stress_List_700K_10ms[2], LogRate_10ms_700K_3GPa], [700, Average_Shear_Stress_List_700K_10ms[3], LogRate_10ms_700K_4GPa], [700, Average_Shear_Stress_List_700K_10ms[4], LogRate_10ms_700K_5GPa]]
#
# for item in data:
#     x_data.append(item[0])
#     y_data.append(item[1])
#     z_data.append(item[2])
# #
# #
# parameters, covariance = optimize.curve_fit(function, [x_data, y_data], z_data)
#
# # create surface function model
# # setup data points for calculating surface model
# model_x_data = np.linspace(min(x_data), max(x_data), 40)
# model_y_data = np.linspace(min(y_data), max(y_data), 40)
#
# # create coordinate arrays for vectorized evaluations
# X, Y = np.meshgrid(model_x_data, model_y_data)
# # calculate Z coordinate array
# Z = function(np.array([X, Y]), *parameters)
#
# z = []
# for row in Z:
#     row.sort()
#     z.append(row)
#
# zlogs = np.array(z)
# import matplotlib.cm
# cm = plt.get_cmap("jet")
# fig = plt.figure()
# ax4 = plt.axes(projection='3d')
# ax4.plot_surface(X, Y, Z, cmap=cm, alpha=0.5, edgecolor='black', linewidth=0.3)
# ax4.set_title('3D Plot - Variation in Log of Dissociation Rates - Alpha Fe, 10ms')
# ax4.scatter(x_data, y_data, z_data, color='black', alpha=1)
# ax4.set_xlabel('Temperature (K)')
# ax4.invert_xaxis()
# ax4.set_ylabel('Shear Stress (GPa)')
# ax4.set_zlabel('Log of Dissociation Rate (per ns)')
# plt.show()
# plt.close(fig)
#
# x_data = []
# y_data = []
# z_data = []
#
# data = [[400, Average_Shear_Stress_List_400K_10ms[0], Dissociation_Rate_10ms_400K_1GPa], [400, Average_Shear_Stress_List_400K_10ms[1], Dissociation_Rate_10ms_400K_2GPa], [400, Average_Shear_Stress_List_400K_10ms[2], Dissociation_Rate_10ms_400K_3GPa], [400, Average_Shear_Stress_List_400K_10ms[3], Dissociation_Rate_10ms_400K_4GPa], [400, Average_Shear_Stress_List_400K_10ms[4], Dissociation_Rate_10ms_400K_5GPa],
#         [500, Average_Shear_Stress_List_500K_10ms[0], Dissociation_Rate_10ms_500K_1GPa], [500, Average_Shear_Stress_List_500K_10ms[1], Dissociation_Rate_10ms_500K_2GPa], [500, Average_Shear_Stress_List_500K_10ms[2], Dissociation_Rate_10ms_500K_3GPa], [500, Average_Shear_Stress_List_500K_10ms[3], Dissociation_Rate_10ms_500K_4GPa], [500, Average_Shear_Stress_List_500K_10ms[4], Dissociation_Rate_10ms_500K_5GPa],
#         [600, Average_Shear_Stress_List_600K_10ms[0], Dissociation_Rate_10ms_600K_1GPa], [600, Average_Shear_Stress_List_600K_10ms[1], Dissociation_Rate_10ms_600K_2GPa], [600, Average_Shear_Stress_List_600K_10ms[2], Dissociation_Rate_10ms_600K_3GPa], [600, Average_Shear_Stress_List_600K_10ms[3], Dissociation_Rate_10ms_600K_4GPa], [600, Average_Shear_Stress_List_600K_10ms[4], Dissociation_Rate_10ms_600K_5GPa],
#         [700, Average_Shear_Stress_List_700K_10ms[0], Dissociation_Rate_10ms_700K_1GPa], [700, Average_Shear_Stress_List_700K_10ms[1], Dissociation_Rate_10ms_700K_2GPa], [700, Average_Shear_Stress_List_700K_10ms[2], Dissociation_Rate_10ms_700K_3GPa], [700, Average_Shear_Stress_List_700K_10ms[3], Dissociation_Rate_10ms_700K_4GPa], [700, Average_Shear_Stress_List_700K_10ms[4], Dissociation_Rate_10ms_700K_5GPa]]
#
# for item in data:
#     x_data.append(item[0])
#     y_data.append(item[1])
#     z_data.append(item[2])
#
# parameters, covariance = optimize.curve_fit(function, [x_data, y_data], z_data)
#
# # create surface function model
# # setup data points for calculating surface model
# model_x_data = np.linspace(min(x_data), max(x_data), 40)
# model_y_data = np.linspace(min(y_data), max(y_data), 40)
# # create coordinate arrays for vectorized evaluations
# X, Y = np.meshgrid(model_x_data, model_y_data)
# # calculate Z coordinate array
# Z = function(np.array([X, Y]), *parameters)
#
# z = []
# for row in Z:
#     row.sort()
#     z.append(row)
#
# z = np.array(z)
# # print(Z)
# # print('####################')
# # print(z)
#
# cm = plt.get_cmap("jet")
# fig = plt.figure()
# ax4 = plt.axes(projection='3d')
# ax4.plot_surface(X, Y, Z, cmap=cm, alpha=0.5, edgecolor= 'black', linewidth=0.3)
# ax4.set_title('3D Plot - Variation in Dissociation Rates - Alpha Fe, 10ms')
# ax4.scatter(x_data, y_data, z_data, color='black', alpha=1)
# ax4.set_xlabel('Temperature (K)')
# ax4.invert_xaxis()
# ax4.set_ylabel('Shear Stress (GPa)')
# ax4.set_zlabel('Dissociation Rate (per ns)')
# plt.show()
#
# # ########## CALCULATING THE ACTIVATION ENERGY, VOLUME AND PREFACTOR FROM 3D FIT ##########
# #
# logRates3D_10ms = zlogs
# shearstresses10ms = model_y_data
# Index = 0
# sigma = []
# params = []
# while Index < len(shearstresses10ms) - 1:
#     for logRates3Drow in logRates3D_10ms:
#         coef_sh_cf, coef_sh_pcov = optimize.curve_fit(linear, shearstresses10ms, logRates3Drow)
#         paramsrow = np.polyfit(shearstresses10ms, logRates3Drow, 1)
#         sigmarow = coef_sh_pcov[0, 0] ** 0.5
#         sigma.append(sigmarow)
#         params.append(paramsrow[0])
#         Index +=1
#
# #
# sigma = np.array(sigma)
# params = np.array(params)
# params = np.average(params)
# sigma = np.average(sigma)
#
# activation_vol_3D_10ms = (params) * (1.38065) * 500 * 1e-2
# print(f'3D Activation Volume, 10ms is {activation_vol_3D_10ms}')
#
# alpha = 0.05
# sigma = sigma ** 0.5
# dof = 2
# tval = t.ppf(1.0 - alpha / 2., dof)
# uncert = sigma * tval
# uncert_ActivationVolume_10ms = uncert * (1.38065) * 500 * 1e-2
# print(f'Activation Volume uncertainty for 3D fit at 10ms is {uncert_ActivationVolume_10ms}')
#
# logRates3D_10ms = zlogs
# temperatures = model_x_data
#
# inverse_temperatures = [1 / x for x in temperatures]
#
# Index = 0
# sigma = []
# SigmaA = []
# params = []
# interceptaverage = []
# while Index < len(inverse_temperatures) - 1:
#     for logRates3Drow in logRates3D_10ms:
#         coef_sh_cf, coef_sh_pcov = optimize.curve_fit(linear, inverse_temperatures, logRates3Drow)
#         paramsrow = np.polyfit(inverse_temperatures, logRates3Drow, 1)
#         sigmarow = coef_sh_pcov[0, 0] ** 0.5
#         sigma_A = coef_sh_pcov[1, 1] ** 0.5
#         sigma.append(sigmarow)
#         SigmaA.append(sigma_A)
#         params.append(paramsrow[0])
#         intercept = paramsrow[1]
#         interceptaverage.append((intercept))
#         Index +=1
#
# sigma = np.array(sigma)
# params = np.array(params)
# SigmaA = np.array(SigmaA)
# interceptaverage = np.array(interceptaverage)
# params_10ms = np.average(params)
# sigma = np.average(sigma)
# interceptaverage = np.average(interceptaverage)
# alpha = 0.05
# sigma = sigma ** 0.5
# tval = t.ppf(1.0 - alpha / 2., dof)
# uncert = sigma * tval
#
# Mu1GPa_10ms = np.average(Average_Mu_List_1GPa_10ms)
# Mu2GPa_10ms = np.average(Average_Mu_List_2GPa_10ms)
# Mu3GPa_10ms = np.average(Average_Mu_List_3GPa_10ms)
# Mu4GPa_10ms = np.average(Average_Mu_List_4GPa_10ms)
# Mu5GPa_10ms = np.average(Average_Mu_List_5GPa_10ms)
#
# MuAveragesDifferentPressures = np.array([Mu1GPa_10ms, Mu2GPa_10ms, Mu3GPa_10ms, Mu4GPa_10ms, Mu5GPa_10ms])
# AverageMu_10ms = np.average(MuAveragesDifferentPressures)
#
# ActivationEnergy_3D_10ms = (((3 * 1e9 * (AverageMu_10ms * activation_vol_3D_10ms) * 1e-30) - 1.381 * params_10ms * 1e-23) * 6.02214076 * (10**23)) / 1000
#
# print(f'Activation Energy for 3D fit at 10ms is {ActivationEnergy_3D_10ms}')
#
# error_3D_10ms = (np.sqrt((3 * 1e9 * np.average(AverageMu_10ms) * uncert_ActivationVolume_10ms * 1e-30) ** 2 + (1.381 * uncert * 1e-23) ** 2) * 6.02214076 * (10**23)) / 1000
# print(f"Activation_Energy Error at 10ms is {error_3D_10ms}")
#
# uncert_prefactor_3D_10ms = sigma_A * tval
# lnA_3D_10ms = np.log((np.exp(interceptaverage) * (10 ** 9)))
#
# print(f"ln(A) for 3D fit is {lnA_3D_10ms}")
# print(f"ln(A) uncertainty is {uncert_prefactor_3D_10ms}" + str(10 ** 9))
#
# ######### Plotting 3D fit vs 2D fit results along with their error margins ##########
# """
# Need to get a list with:
# - Values for Activation Volumes at different temperatures
# - Values for Activation Energies at different pressures
# - Values for ln(A) at different pressures
# - List of errors for each of the above quantities
# - Make a graph for each surface chemistry
# - Eventually will need to do the same for the different sliding speeds
#
# """
# Temperatures = [400, 500, 600, 700]
# Pressures = [1, 2, 3, 4, 5]
#
#
# Activation_Energies_10ms = [ActivationEnergy_1GPa_10ms, ActivationEnergy_2GPa_10ms, ActivationEnergy_3GPa_10ms, ActivationEnergy_4GPa_10ms, ActivationEnergy_5GPa_10ms]
# Activation_Energy_Errors_10ms = [err_E_1GPa_10ms, err_E_2GPa_10ms, err_E_3GPa_10ms, err_E_4GPa_10ms, err_E_5GPa_10ms]
# ActivationEnergy_3D_10ms = ActivationEnergy_3D_10ms
# ActivationEnergy_3D_error_10ms = error_3D_10ms
# ActivationEnergy_3D_error_UpperBound_Value_10ms = float(float(ActivationEnergy_3D_10ms) + float(ActivationEnergy_3D_error_10ms))
# ActivationEnergy_3D_error_LowerBound_Value_10ms = float(float(ActivationEnergy_3D_10ms) - float(ActivationEnergy_3D_error_10ms))
#
# Activation_Energy_Error_Plot, Ea2Dvs3d  = plt.subplots()
# Ea2Dvs3d.set_title('Comparison of Activation Energies from 2D and 3D Fits - AlphaFe, 10ms')
# Ea2Dvs3d.set_xlabel('Normal Stress (GPa)')
# Ea2Dvs3d.set_ylabel('Activation Energy')
# Ea2Dvs3d.scatter(Pressures, Activation_Energies_10ms)
# Ea2Dvs3d.errorbar(Pressures, Activation_Energies_10ms, yerr=Activation_Energy_Errors_10ms, linestyle="None", fmt='o', capsize=3)
# Ea2Dvs3d.axhline(y=ActivationEnergy_3D_10ms)
# Pressures = [0.5, 1, 2, 3, 4, 5, 5.5]
# Ea2Dvs3d.fill_between(Pressures, ActivationEnergy_3D_error_LowerBound_Value_10ms, ActivationEnergy_3D_error_UpperBound_Value_10ms, alpha=0.4)
# Ea2Dvs3d.set_xlim(0.5, 5.5)
# Ea2Dvs3d.set_ylim(0, 35)
# plt.show()
#
# ################ Activation Volume Errors #############################
#
# Activation_Volumes_10ms = [activation_vol_400K_10ms, activation_vol_500K_10ms, activation_vol_600K_10ms, activation_vol_700K_10ms]
# Activation_Volume_Errors_10ms = [uncert400_10ms, uncert500_10ms, uncert600_10ms, uncert700_10ms]
# Activation_Volume_3D_10ms = activation_vol_3D_10ms
# Activation_Volume_3D_Error_10ms = uncert_ActivationVolume_10ms
# ActivationVolume_3D_error_UpperBound_Value_10ms = float(float(Activation_Volume_3D_10ms) + float(Activation_Volume_3D_Error_10ms))
# ActivationVolume_3D_error_LowerBound_Value_10ms = float(float(Activation_Volume_3D_10ms) - float(Activation_Volume_3D_Error_10ms))
#
# Activation_Volume_Error_Plot, Av2Dvs3d  = plt.subplots()
# Av2Dvs3d.set_title('Comparison of Activation Volumes from 2D and 3D Fits - AlphaFe, 10ms')
# Av2Dvs3d.set_xlabel('Normal Stress(GPa)')
# Av2Dvs3d.set_ylabel('Activation Volume')
# Av2Dvs3d.scatter(Temperatures, Activation_Volumes_10ms)
# Av2Dvs3d.errorbar(Temperatures, Activation_Volumes_10ms, yerr=Activation_Volume_Errors_10ms, linestyle="None", fmt='o', capsize=3)
# Av2Dvs3d.axhline(y=Activation_Volume_3D_10ms)
# Temperatures = [350, 400, 500, 600, 700, 750]
# Av2Dvs3d.fill_between(Temperatures, ActivationVolume_3D_error_LowerBound_Value_10ms, ActivationVolume_3D_error_UpperBound_Value_10ms, alpha=0.4)
# Av2Dvs3d.set_xlim(350, 750)
# Av2Dvs3d.set_ylim(0, 30)
# plt.show()
#
# #################### Prefactor Errors #########################
# Pressures = [1, 2, 3, 4, 5]
# Prefactors_10ms = [lnA1GPa_10ms, lnA2GPa_10ms, lnA3GPa_10ms, lnA4GPa_10ms, lnA5GPa_10ms]
# PrefactorErrors_10ms = [uncert_A1GPa_10ms, uncert_A2GPa_10ms, uncert_A3GPa_10ms, uncert_A4GPa_10ms, uncert_A5GPa_10ms]
# lnA_3D_error_10ms = uncert_prefactor_3D_10ms
# lnA_3D_error_UpperBound_Value_10ms = float(float(lnA_3D_10ms) + float(lnA_3D_error_10ms))
# lnA_3D_error_LowerBound_Value_10ms = float(float(lnA_3D_10ms) - float(lnA_3D_error_10ms))
#
# Prefactor_Error_Plot, lnA2Dvs3d  = plt.subplots()
# lnA2Dvs3d.set_title('Comparison of Prefactors from 2D and 3D Fits - AlphaFe, 10ms')
# lnA2Dvs3d.set_xlabel('Normal Stress(GPa)')
# lnA2Dvs3d.set_ylabel('Prefactor')
# lnA2Dvs3d.scatter(Pressures, Prefactors_10ms)
# lnA2Dvs3d.errorbar(Pressures, Prefactors_10ms, yerr=PrefactorErrors_10ms, linestyle="None", fmt='o', capsize=3)
# lnA2Dvs3d.axhline(y=lnA_3D_10ms)
# Pressures = [0.5, 1, 2, 3, 4, 5, 5.5]
# lnA2Dvs3d.fill_between(Pressures, lnA_3D_error_LowerBound_Value_10ms, lnA_3D_error_UpperBound_Value_10ms, alpha=0.4)
# lnA2Dvs3d.set_xlim(0.5, 5.5)
# lnA2Dvs3d.set_ylim(20, 30)
# plt.show()
#
# ########## Checking for Kinetic Compensation Effect ##########
#
# KineticCompeEffectPlot, EavsLnA  = plt.subplots()
# EavsLnA.set_title('Activation Energy vs Prefactor - AlphaFe, 10ms')
# EavsLnA.set_xlabel('Ea')
# EavsLnA.set_ylabel('ln A')
# EavsLnA.scatter(Activation_Energies_10ms, Prefactors_10ms)
# EavsLnA.set_ylim(24, 27)
# plt.show()

#############################  20ms  #################################
# Index = 0
# Temperatures = ["500K", "600K", "700K"]
# Pressures = ['2GPa', '3GPa', '4GPa', '5GPa']
# Speeds = ['20ms', '30ms', '40ms', '50ms']
#
# Dissociation_Rate_20ms_400K_1GPa, LogRate_20ms_400K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_400K_1GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_20ms_400K_2GPa, LogRate_20ms_400K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_400K_2GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_20ms_400K_3GPa, LogRate_20ms_400K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_400K_3GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_20ms_400K_4GPa, LogRate_20ms_400K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_400K_4GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_20ms_400K_5GPa, LogRate_20ms_400K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_400K_5GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rates_20ms_400K = [Dissociation_Rate_20ms_400K_1GPa, Dissociation_Rate_20ms_400K_2GPa, Dissociation_Rate_20ms_400K_3GPa, Dissociation_Rate_20ms_400K_4GPa, Dissociation_Rate_20ms_400K_5GPa]
# Log_Rates_400K_20ms = [LogRate_20ms_400K_1GPa, LogRate_20ms_400K_2GPa, LogRate_20ms_400K_3GPa, LogRate_20ms_400K_4GPa, LogRate_20ms_400K_5GPa]
#
# Dissociation_Rate_20ms_500K_1GPa, LogRate_20ms_500K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_500K_1GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_20ms_500K_2GPa, LogRate_20ms_500K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_500K_2GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_20ms_500K_3GPa, LogRate_20ms_500K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_500K_3GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_20ms_500K_4GPa, LogRate_20ms_500K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_500K_4GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_20ms_500K_5GPa, LogRate_20ms_500K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_500K_5GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rates_20ms_500K = [Dissociation_Rate_20ms_500K_1GPa, Dissociation_Rate_20ms_500K_2GPa, Dissociation_Rate_20ms_500K_3GPa, Dissociation_Rate_20ms_500K_4GPa, Dissociation_Rate_20ms_500K_5GPa]
# Log_Rates_500K_20ms = [LogRate_20ms_500K_1GPa, LogRate_20ms_500K_2GPa, LogRate_20ms_500K_3GPa, LogRate_20ms_500K_4GPa, LogRate_20ms_500K_5GPa]
#
# Dissociation_Rate_20ms_600K_1GPa, LogRate_20ms_600K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_600K_1GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_20ms_600K_2GPa, LogRate_20ms_600K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_600K_2GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_20ms_600K_3GPa, LogRate_20ms_600K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_600K_3GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_20ms_600K_4GPa, LogRate_20ms_600K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_600K_4GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_20ms_600K_5GPa, LogRate_20ms_600K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_600K_5GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rates_20ms_600K = [Dissociation_Rate_20ms_600K_1GPa, Dissociation_Rate_20ms_600K_2GPa, Dissociation_Rate_20ms_600K_3GPa, Dissociation_Rate_20ms_600K_4GPa, Dissociation_Rate_20ms_600K_5GPa]
# Log_Rates_600K_20ms = [LogRate_20ms_600K_1GPa, LogRate_20ms_600K_2GPa, LogRate_20ms_600K_3GPa, LogRate_20ms_600K_4GPa, LogRate_20ms_600K_5GPa]
#
# Dissociation_Rate_20ms_700K_1GPa, LogRate_20ms_700K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_700K_1GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_20ms_700K_2GPa, LogRate_20ms_700K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_700K_2GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_20ms_700K_3GPa, LogRate_20ms_700K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_700K_3GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_20ms_700K_4GPa, LogRate_20ms_700K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_700K_4GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_20ms_700K_5GPa, LogRate_20ms_700K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_20ms_700K_5GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rates_20ms_700K = [Dissociation_Rate_20ms_700K_1GPa, Dissociation_Rate_20ms_700K_2GPa, Dissociation_Rate_20ms_700K_3GPa, Dissociation_Rate_20ms_700K_4GPa, Dissociation_Rate_20ms_700K_5GPa]
# Log_Rates_700K_20ms = [LogRate_20ms_700K_1GPa, LogRate_20ms_700K_2GPa, LogRate_20ms_700K_3GPa, LogRate_20ms_700K_4GPa, LogRate_20ms_700K_5GPa]
#
# Dissociation_Rates_1GPa = [Dissociation_Rates_20ms_400K[0], Dissociation_Rates_20ms_500K[0], Dissociation_Rates_20ms_600K[0], Dissociation_Rates_20ms_700K[0]]
# Dissociation_Rates_2GPa = [Dissociation_Rates_20ms_400K[1], Dissociation_Rates_20ms_500K[1], Dissociation_Rates_20ms_600K[1], Dissociation_Rates_20ms_700K[1]]
# Dissociation_Rates_3GPa = [Dissociation_Rates_20ms_400K[2], Dissociation_Rates_20ms_500K[2], Dissociation_Rates_20ms_600K[2], Dissociation_Rates_20ms_700K[2]]
# Dissociation_Rates_4GPa = [Dissociation_Rates_20ms_400K[3], Dissociation_Rates_20ms_500K[3], Dissociation_Rates_20ms_600K[3], Dissociation_Rates_20ms_700K[3]]
# Dissociation_Rates_5GPa = [Dissociation_Rates_20ms_400K[4], Dissociation_Rates_20ms_500K[4], Dissociation_Rates_20ms_600K[4], Dissociation_Rates_20ms_700K[4]]
#
# Log_Rates_1GPa_20ms = [LogRate_20ms_400K_1GPa, LogRate_20ms_500K_1GPa, LogRate_20ms_600K_1GPa, LogRate_20ms_700K_1GPa]
# Log_Rates_2GPa_20ms = [LogRate_20ms_400K_2GPa, LogRate_20ms_500K_2GPa, LogRate_20ms_600K_2GPa, LogRate_20ms_700K_2GPa]
# Log_Rates_3GPa_20ms = [LogRate_20ms_400K_3GPa, LogRate_20ms_500K_3GPa, LogRate_20ms_600K_3GPa, LogRate_20ms_700K_3GPa]
# Log_Rates_4GPa_20ms = [LogRate_20ms_400K_4GPa, LogRate_20ms_500K_4GPa, LogRate_20ms_600K_4GPa, LogRate_20ms_700K_4GPa]
# Log_Rates_5GPa_20ms = [LogRate_20ms_400K_5GPa, LogRate_20ms_500K_5GPa, LogRate_20ms_600K_5GPa, LogRate_20ms_700K_5GPa]
#
# EquilibriumFactor = [79, 99] # How many rows (out of 99) to ignore before calculating shear stress/friction coefficient, as it won't stabilise until after a certain number of timesteps
#
# def get_average_shear_normal_stress_and_average_mu_constant_temperature(Temperature, Pressures, EquilibriumFactor, Speed):
#     Friction_Coefficient_Dataframe_Unnamed = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{Speed}/{Temperature}/1GPa/'
#                                 'fc_ave.dump'.format(Speed=Speed, Temperature=Temperature), sep=' ')
#     Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe_Unnamed.rename(columns={'v_s_bot' : 'Shear Stress 1GPa', 'v_p_bot' : 'Normal Stress 1GPa'})
#
#     for P in Pressures:
#         Dataframe = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{Speed}/{Temperature}/{P}/'
#                                 'fc_ave.dump'.format(Speed=Speed, Temperature=Temperature, P=P), sep=' ')
#         Big_DataframeP = Dataframe.rename(columns= {'Timestep': 'Timestep {}'.format(P),
#                                                         'v_s_bot': 'Shear Stress {}'.format(P),
#                                                         'v_p_bot': 'Normal Stress {}'.format(P)})
#
#         Friction_Coefficient_Dataframe = pd.concat([Friction_Coefficient_Dataframe, Big_DataframeP], axis =1)
#         Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe.dropna()
#
#
#     #print(Friction_Coefficient_Dataframe)
#     Mu_Final_Dataframe = Friction_Coefficient_Dataframe.iloc[:, [0, 1, 2, 4, 5, 7, 8, 10, 11, 13, 14]]
#     Mu_Final_Dataframe = Mu_Final_Dataframe.iloc[EquilibriumFactor[0]:EquilibriumFactor[1], :]
#     #print(Mu_Final_Dataframe)
#
#     ShearStressMeans = Mu_Final_Dataframe[['Shear Stress 1GPa', 'Shear Stress 2GPa', 'Shear Stress 3GPa', 'Shear Stress 4GPa', 'Shear Stress 5GPa']].mean()
#     Average_Shear_Stress_Dictionary = ShearStressMeans.to_dict()
#     #print(ShearStressMeans)
#     NormalStressMeans = Mu_Final_Dataframe[['Normal Stress 1GPa', 'Normal Stress 2GPa', 'Normal Stress 3GPa', 'Normal Stress 4GPa', 'Normal Stress 5GPa']].mean()
#     NormalStressMeans = NormalStressMeans.to_dict()
#     #print(NormalStressMeans)
#
#     Average_Mu_Dictionary = {}
#
#     NormalStressMeansList = list(NormalStressMeans.values())
#     Normal_Stress = NormalStressMeans.get('Normal Stress 1GPa')
#     Shear_Stress = ShearStressMeans.get('Shear Stress 1GPa')
#     Average_Mu = Shear_Stress / Normal_Stress
#     Average_Mu_Dictionary.update({'Average Mu 1GPa': Average_Mu})
#
#     for P in Pressures:
#
#         Normal_Stress = NormalStressMeans.get('Normal Stress {}'.format(P))
#         Shear_Stress = ShearStressMeans.get('Shear Stress {}'.format(P))
#         Average_Mu = Shear_Stress / Normal_Stress
#         Average_Mu_Dictionary.update({'Average Mu {}'.format(P): Average_Mu})
#
#
#     Average_Shear_Stress_List = list(Average_Shear_Stress_Dictionary.values())
#     #print(Average_Shear_Stress_List)
#     Average_Mu_List = list(Average_Mu_Dictionary.values())
#     Average_Shear_Stress_List = [x / 10000 for x in Average_Shear_Stress_List] # Conversion to GPa
#     NormalStressMeansList = [x/10000 for x in NormalStressMeansList]
#     #print(Average_Shear_Stress_List)
#
#     return Average_Shear_Stress_List, Average_Mu_List, NormalStressMeansList
#
# ######## Getting Average Shear Stress, Friction Coefficient and Normal Stress #################
# Average_Shear_Stress_List_400K_20ms, Average_Mu_List_400K_20ms, NormalStressMeans_400K_20ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("400K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="20ms")
# Average_Shear_Stress_List_500K_20ms, Average_Mu_List_500K_20ms, NormalStressMeans_500K_20ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("500K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="20ms")
# Average_Shear_Stress_List_600K_20ms, Average_Mu_List_600K_20ms, NormalStressMeans_600K_20ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("600K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="20ms")
# Average_Shear_Stress_List_700K_20ms, Average_Mu_List_700K_20ms, NormalStressMeans_700K_20ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("700K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="20ms")
#
# Average_Mu_List_1GPa_20ms = [Average_Mu_List_400K_20ms[0], Average_Mu_List_500K_20ms[0], Average_Mu_List_600K_20ms[0], Average_Mu_List_700K_20ms[0]]
# Average_Mu_List_2GPa_20ms = [Average_Mu_List_400K_20ms[1], Average_Mu_List_500K_20ms[1], Average_Mu_List_600K_20ms[1], Average_Mu_List_700K_20ms[1]]
# Average_Mu_List_3GPa_20ms = [Average_Mu_List_400K_20ms[2], Average_Mu_List_500K_20ms[2], Average_Mu_List_600K_20ms[2], Average_Mu_List_700K_20ms[2]]
# Average_Mu_List_4GPa_20ms = [Average_Mu_List_400K_20ms[3], Average_Mu_List_500K_20ms[3], Average_Mu_List_600K_20ms[3], Average_Mu_List_700K_20ms[3]]
# Average_Mu_List_5GPa_20ms = [Average_Mu_List_400K_20ms[4], Average_Mu_List_500K_20ms[4], Average_Mu_List_600K_20ms[4], Average_Mu_List_700K_20ms[4]]
#
# Average_Shear_Stress_List_1GPa_20ms = [Average_Shear_Stress_List_400K_20ms[0], Average_Shear_Stress_List_500K_20ms[0], Average_Shear_Stress_List_600K_20ms[0], Average_Shear_Stress_List_700K_20ms[0]]
# Average_Shear_Stress_List_2GPa_20ms = [Average_Shear_Stress_List_400K_20ms[1], Average_Shear_Stress_List_500K_20ms[1], Average_Shear_Stress_List_600K_20ms[1], Average_Shear_Stress_List_700K_20ms[1]]
# Average_Shear_Stress_List_3GPa_20ms = [Average_Shear_Stress_List_400K_20ms[2], Average_Shear_Stress_List_500K_20ms[2], Average_Shear_Stress_List_600K_20ms[2], Average_Shear_Stress_List_700K_20ms[2]]
# Average_Shear_Stress_List_4GPa_20ms = [Average_Shear_Stress_List_400K_20ms[3], Average_Shear_Stress_List_500K_20ms[3], Average_Shear_Stress_List_600K_20ms[3], Average_Shear_Stress_List_700K_20ms[3]]
# Average_Shear_Stress_List_5GPa_20ms = [Average_Shear_Stress_List_400K_20ms[4], Average_Shear_Stress_List_500K_20ms[4], Average_Shear_Stress_List_600K_20ms[4], Average_Shear_Stress_List_700K_20ms[4]]
#
# plot_shear_stress_vs_normal_stress(Average_Shear_Stress_List_400K_20ms, Average_Shear_Stress_List_500K_20ms, Average_Shear_Stress_List_600K_20ms, Average_Shear_Stress_List_700K_20ms,
#                                  "400K", "500K", "600K", "700K", Speed="20ms")
#
# ########################## Plotting ln(Rates) vs Shear Stress #####################################
# x = np.array([0, 1, 2, 3, 4, 5])
# params400K_20ms = np.polyfit(Average_Shear_Stress_List_400K_20ms, Log_Rates_400K_20ms, 1)
# params500K_20ms = np.polyfit(Average_Shear_Stress_List_500K_20ms, Log_Rates_500K_20ms, 1)
# params600K_20ms = np.polyfit(Average_Shear_Stress_List_600K_20ms, Log_Rates_600K_20ms, 1)
# params700K_20ms = np.polyfit(Average_Shear_Stress_List_700K_20ms, Log_Rates_700K_20ms, 1)
#
# RatesvsShear, RvS3 = plt.subplots()
# RvS3.set_title('Log of Dissociation Rates vs Shear Stress, 20ms')
# RvS3.set_xlabel('Shear Stress(GPa)')
# RvS3.set_ylabel('Log of Dissociation Rate (per nanosecond)')
# RvS3.scatter(Average_Shear_Stress_List_400K_20ms, Log_Rates_400K_20ms)
# RvS3.scatter(Average_Shear_Stress_List_500K_20ms, Log_Rates_500K_20ms)
# RvS3.scatter(Average_Shear_Stress_List_600K_20ms, Log_Rates_600K_20ms)
# RvS3.scatter(Average_Shear_Stress_List_700K_20ms, Log_Rates_700K_20ms)
# RvS3.plot(x, params400K_20ms[0] * x + params400K_20ms[1], label='400K Fitted')
# RvS3.plot(x, params500K_20ms[0] * x + params500K_20ms[1], label='500K Fitted')
# RvS3.plot(x, params600K_20ms[0] * x + params600K_20ms[1], label='600K Fitted')
# RvS3.plot(x, params700K_20ms[0] * x + params700K_20ms[1], label='700K Fitted')
# RvS3.set_xlim(0, 2)
# RvS3.set_ylim(-0.5, 4)
# RvS3.legend(loc='lower right')
# plt.show()
#
# ####### Calculate Activation Volume, Using  Carlos' conversion to get in Angstrom^3 #################
#
# activation_vol_400K_20ms = (params400K_20ms[0]) * (1.38065) * 400 * 1e-2
# activation_vol_500K_20ms = (params500K_20ms[0]) * (1.38065) * 500 * 1e-2
# activation_vol_600K_20ms = (params600K_20ms[0]) * (1.38065) * 600 * 1e-2
# activation_vol_700K_20ms = (params700K_20ms[0]) * (1.38065) * 700 * 1e-2
#
# alpha = 0.05
# def linear(x, m, n):
#     return m * x + n
#
# coef_sh_cf, coef_sh_pcov = optimize.curve_fit(linear, Average_Shear_Stress_List_700K_20ms, Log_Rates_700K_20ms)
# sigma = coef_sh_pcov[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_700K_20ms) - len(params700K_20ms))
#
# tval = t.ppf(1.0 - alpha / 2., dof)
# uncert = sigma * tval
# uncert400_20ms = uncert * (1.38065) * 400 * 1e-2
# uncert500_20ms = uncert * (1.38065) * 500 * 1e-2
# uncert600_20ms = uncert * (1.38065) * 600 * 1e-2
# uncert700_20ms = uncert * (1.38065) * 700 * 1e-2
#
# print(f'Activation Volume uncertainty at 400K, 20ms is {uncert400_20ms}')
# print(f'Activation Volume uncertainty at 500K, 20ms is {uncert500_20ms}')
# print(f'Activation Volume uncertainty at 600K, 20ms is {uncert600_20ms}')
# print(f'Activation Volume uncertainty at 700K, 20ms is {uncert700_20ms}')
#
# Activation_Volumes_20ms = [activation_vol_400K_20ms, activation_vol_500K_20ms, activation_vol_600K_20ms, activation_vol_700K_20ms]
# mean_actv_20ms = np.average(Activation_Volumes_20ms)
#
# print('Activation Volume 400K, 20ms = ' + str(activation_vol_400K_20ms))
# print('Activation Volume 500K, 20ms = ' + str(activation_vol_500K_20ms))
# print('Activation Volume 600K, 20ms = ' + str(activation_vol_600K_20ms))
# print('Activation Volume 700K, 20ms = ' + str(activation_vol_700K_20ms))
# #
# ############ Plotting lnk vs 1000/T  #########################
#
# Temperatures = [400, 500, 600, 700]
# Inverse_Temperatures = np.array([1/x for x in Temperatures])
#
# trend1GPa_20ms = np.polyfit(Inverse_Temperatures, Log_Rates_1GPa_20ms, 1)
# trend2GPa_20ms = np.polyfit(Inverse_Temperatures, Log_Rates_2GPa_20ms, 1)
# trend3GPa_20ms = np.polyfit(Inverse_Temperatures, Log_Rates_3GPa_20ms, 1)
# trend4GPa_20ms = np.polyfit(Inverse_Temperatures, Log_Rates_4GPa_20ms, 1)
# trend5GPa_20ms = np.polyfit(Inverse_Temperatures, Log_Rates_5GPa_20ms, 1)
#
# fig1, ax1 = plt.subplots()
# ax1.set_title('Log of Dissociation Rates against Inverse of Temperatures, 20ms')
# ax1.set_xlabel('1000/T (K-1)')
# ax1.set_ylabel('ln(Rate) (ns-1)')
# ax1.scatter(Inverse_Temperatures, Log_Rates_1GPa_20ms)
# ax1.scatter(Inverse_Temperatures, Log_Rates_2GPa_20ms)
# ax1.scatter(Inverse_Temperatures, Log_Rates_3GPa_20ms)
# ax1.scatter(Inverse_Temperatures, Log_Rates_4GPa_20ms)
# ax1.scatter(Inverse_Temperatures, Log_Rates_5GPa_20ms)
#
# Fit1GPa_20ms = np.poly1d(trend1GPa_20ms)
# Fit2GPa_20ms = np.poly1d(trend2GPa_20ms)
# Fit3GPa_20ms = np.poly1d(trend3GPa_20ms)
# Fit4GPa_20ms = np.poly1d(trend4GPa_20ms)
# Fit5GPa_20ms = np.poly1d(trend5GPa_20ms)
#
# ax1.plot(Inverse_Temperatures, Fit1GPa_20ms(Inverse_Temperatures), label='1GPa')
# ax1.plot(Inverse_Temperatures, Fit2GPa_20ms(Inverse_Temperatures), label='2GPa')
# ax1.plot(Inverse_Temperatures, Fit3GPa_20ms(Inverse_Temperatures), label='3GPa')
# ax1.plot(Inverse_Temperatures, Fit4GPa_20ms(Inverse_Temperatures), label='4GPa')
# ax1.plot(Inverse_Temperatures, Fit5GPa_20ms(Inverse_Temperatures), label='5GPa')
# ax1.legend()
# plt.show()
#
# ActivationEnergy_1GPa_20ms = (((1 * 1e9 * (np.average(Average_Mu_List_1GPa_20ms) * mean_actv_20ms) * 1e-30) - 1.381 * trend1GPa_20ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
# ActivationEnergy_2GPa_20ms = (((2 * 1e9 * (np.average(Average_Mu_List_2GPa_20ms) * mean_actv_20ms) * 1e-30) - 1.381 * trend2GPa_20ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
# ActivationEnergy_3GPa_20ms = (((3 * 1e9 * (np.average(Average_Mu_List_3GPa_20ms) * mean_actv_20ms) * 1e-30) - 1.381 * trend3GPa_20ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
# ActivationEnergy_4GPa_20ms = (((4 * 1e9 * (np.average(Average_Mu_List_4GPa_20ms) * mean_actv_20ms) * 1e-30) - 1.381 * trend4GPa_20ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
# ActivationEnergy_5GPa_20ms = (((5 * 1e9 * (np.average(Average_Mu_List_5GPa_20ms) * mean_actv_20ms) * 1e-30) - 1.381 * trend5GPa_20ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
#
# print(np.average(Average_Mu_List_1GPa_20ms))
# print(np.average(Average_Mu_List_2GPa_20ms))
# print(np.average(Average_Mu_List_3GPa_20ms))
# print(np.average(Average_Mu_List_4GPa_20ms))
# print(np.average(Average_Mu_List_5GPa_20ms))
#
# lnA1GPa_20ms = np.log((np.exp(trend1GPa_20ms[1]) * (10 ** 9)))
# lnA2GPa_20ms = np.log((np.exp(trend2GPa_20ms[1]) * (10 ** 9)))
# lnA3GPa_20ms = np.log((np.exp(trend3GPa_20ms[1]) * (10 ** 9)))
# lnA4GPa_20ms = np.log((np.exp(trend4GPa_20ms[1]) * (10 ** 9)))
# lnA5GPa_20ms = np.log((np.exp(trend5GPa_20ms[1]) * (10 ** 9)))
#
# print(f"ln(A) at 1GPa, 20ms is {lnA1GPa_20ms}")
# print(f"ln(A) at 2GPa, 20ms is {lnA2GPa_20ms}")
# print(f"ln(A) at 3GPa, 20ms is {lnA3GPa_20ms}")
# print(f"ln(A) at 4GPa, 20ms is {lnA4GPa_20ms}")
# print(f"ln(A) at 5GPa, 20ms is {lnA5GPa_20ms}")
#
# print('Activation Energy 1GPa, 20ms =' + str(ActivationEnergy_1GPa_20ms))
# print('Activation Energy 2GPa, 20ms =' + str(ActivationEnergy_2GPa_20ms))
# print('Activation Energy 3GPa, 20ms =' + str(ActivationEnergy_3GPa_20ms))
# print('Activation Energy 4GPa, 20ms =' + str(ActivationEnergy_4GPa_20ms))
# print('Activation Energy 5GPa, 20ms =' + str(ActivationEnergy_5GPa_20ms))
#
# def linear(x, m, n):
#     return m * x + n
#
# coef_p_cf_1GPa_20ms, coef_p_pcov_1GPa_20ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_1GPa_20ms)
# sigma_A1GPa_20ms = coef_p_pcov_1GPa_20ms[1, 1] ** 0.5
# sigma_m1GPa_20ms = coef_p_pcov_1GPa_20ms[0, 0] ** 0.5
# dof_20ms = max(0, len(Log_Rates_1GPa_20ms) - len(coef_p_cf_1GPa_20ms))
#
# alpha = 0.05
# tval_20ms = t.ppf(1.0 - alpha / 2., 4.1532)
#
# sigma_20ms = np.std(Activation_Volumes_20ms)
# error_actv_20ms = sigma_20ms * tval_20ms / np.sqrt(len(Activation_Volumes_20ms))
# uncert_A1GPa_20ms = sigma_A1GPa_20ms * tval
# uncert_m1GPa_20ms = sigma_m1GPa_20ms * tval
#
# coef_p_cf_2GPa_20ms, coef_p_pcov_2GPa_20ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_2GPa_20ms)
# sigma_A2GPa_20ms = coef_p_pcov_2GPa_20ms[1, 1] ** 0.5
# sigma_m2GPa_20ms = coef_p_pcov_2GPa_20ms[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_2GPa_20ms) - len(coef_p_cf_2GPa_20ms))
# uncert_A2GPa_20ms = sigma_A2GPa_20ms * tval_20ms
# uncert_m2GPa_20ms = sigma_m2GPa_20ms * tval_20ms
#
# coef_p_cf_3GPa_20ms, coef_p_pcov_3GPa_20ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_3GPa_20ms)
# sigma_A3GPa_20ms = coef_p_pcov_3GPa_20ms[1, 1] ** 0.5
# sigma_m3GPa_20ms = coef_p_pcov_3GPa_20ms[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_3GPa_20ms) - len(coef_p_cf_3GPa_20ms))
# uncert_A3GPa_20ms = sigma_A3GPa_20ms * tval_20ms
# uncert_m3GPa_20ms = sigma_m3GPa_20ms * tval_20ms
#
# coef_p_cf_4GPa_20ms, coef_p_pcov_4GPa_20ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_4GPa_20ms)
# sigma_A4GPa_20ms = coef_p_pcov_4GPa_20ms[1, 1] ** 0.5
# sigma_m4GPa_20ms = coef_p_pcov_4GPa_20ms[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_4GPa_20ms) - len(coef_p_cf_4GPa_20ms))
# uncert_A4GPa_20ms = sigma_A4GPa_20ms * tval_20ms
# uncert_m4GPa_20ms = sigma_m4GPa_20ms * tval_20ms
#
# coef_p_cf_5GPa_20ms, coef_p_pcov_5GPa_20ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_5GPa_20ms)
# sigma_A5GPa_20ms = coef_p_pcov_5GPa_20ms[1, 1] ** 0.5
# sigma_m5GPa_20ms = coef_p_pcov_5GPa_20ms[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_5GPa_20ms) - len(coef_p_cf_5GPa_20ms))
# uncert_A5GPa_20ms = sigma_A5GPa_20ms * tval_20ms
# uncert_m5GPa_20ms = sigma_m5GPa_20ms * tval_20ms
#
# err_E_1GPa_20ms = (np.sqrt(((1 * 1e9 * np.average(Average_Mu_List_1GPa_20ms) * error_actv_20ms * 1e-30) ** 2 + (1.381 * uncert_m1GPa_20ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
# err_E_2GPa_20ms = (np.sqrt(((2 * 1e9 * np.average(Average_Mu_List_2GPa_20ms) * error_actv_20ms * 1e-30) ** 2 + (1.381 * uncert_m2GPa_20ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
# err_E_3GPa_20ms = (np.sqrt(((3 * 1e9 * np.average(Average_Mu_List_3GPa_20ms) * error_actv_20ms * 1e-30) ** 2 + (1.381 * uncert_m3GPa_20ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
# err_E_4GPa_20ms = (np.sqrt(((4 * 1e9 * np.average(Average_Mu_List_4GPa_20ms) * error_actv_20ms * 1e-30) ** 2 + (1.381 * uncert_m4GPa_20ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
# err_E_5GPa_20ms = (np.sqrt(((5 * 1e9 * np.average(Average_Mu_List_5GPa_20ms) * error_actv_20ms * 1e-30) ** 2 + (1.381 * uncert_m5GPa_20ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
#
# print(f"Error for Activation Energy at 1GPa = {err_E_1GPa_20ms}")
# print(f"Error for Activation Energy at 2GPa = {err_E_2GPa_20ms}")
# print(f"Error for Activation Energy at 3GPa = {err_E_3GPa_20ms}")
# print(f"Error for Activation Energy at 4GPa = {err_E_4GPa_20ms}")
# print(f"Error for Activation Energy at 5GPa = {err_E_5GPa_20ms}")
#
# uncertlnA1GPa_20ms = np.log(uncert_A1GPa_20ms * (10 ** 9))
# uncertlnA2GPa_20ms = np.log(uncert_A2GPa_20ms * (10 ** 9))
# uncertlnA3GPa_20ms = np.log(uncert_A3GPa_20ms * (10 ** 9))
# uncertlnA4GPa_20ms = np.log(uncert_A4GPa_20ms * (10 ** 9))
# uncertlnA5GPa_20ms = np.log(uncert_A5GPa_20ms * (10 ** 9))
#
# print(f"Error for ln A at 1GPa = {uncert_A1GPa_20ms}")
# print(f"Error for ln A at 2GPa = {uncert_A2GPa_20ms}")
# print(f"Error for ln A at 3GPa = {uncert_A3GPa_20ms}")
# print(f"Error for ln A at 4GPa = {uncert_A4GPa_20ms}")
# print(f"Error for ln A at 5GPa = {uncert_A5GPa_20ms}")
#
# pressures = ['2GPa', '3GPa', '4GPa', '5GPa']
#
# ########################## Plotting ln(Rates) vs Normal Stress #####################################
# x = np.array([0, 1, 2, 3, 4, 5])
# params400K_20ms = np.polyfit(NormalStressMeans_400K_20ms, Log_Rates_400K_20ms, 1)
# params500K_20ms = np.polyfit(NormalStressMeans_500K_20ms, Log_Rates_500K_20ms, 1)
# params600K_20ms = np.polyfit(NormalStressMeans_600K_20ms, Log_Rates_600K_20ms, 1)
# params700K_20ms = np.polyfit(NormalStressMeans_700K_20ms, Log_Rates_700K_20ms, 1)
#
# RatesvsNormal, RvN3 = plt.subplots()
# RvN3.set_title('Log of Dissociation Rates vs Normal Stress - Alpha Fe, 20ms')
# RvN3.set_xlabel('Normal Stress(GPa)')
# RvN3.set_ylabel('Log of Dissociation Rate (ns-1)')
# RvN3.scatter(NormalStressMeans_400K_20ms, Log_Rates_400K_20ms)
# RvN3.scatter(NormalStressMeans_500K_20ms, Log_Rates_500K_20ms)
# RvN3.scatter(NormalStressMeans_600K_20ms, Log_Rates_600K_20ms)
# RvN3.scatter(NormalStressMeans_700K_20ms, Log_Rates_700K_20ms)
# RvN3.plot(x, params400K_20ms[0] * x + params400K_20ms[1], label='400K Fitted')
# RvN3.plot(x, params500K_20ms[0] * x + params500K_20ms[1], label='500K Fitted')
# RvN3.plot(x, params600K_20ms[0] * x + params600K_20ms[1], label='600K Fitted')
# RvN3.plot(x, params700K_20ms[0] * x + params700K_20ms[1], label='700K Fitted')
# #RvN3.set_xlim(1, 5)
# # RvN3.set_ylim(0, 25)
# RvN3.legend(loc='lower right')
# plt.show()
#
# def function(data, a, b, c):
#     x = data[0]
#     y = data[1]
#     return a * (x**b) * (y**c) #TODO change fitting function
#
# x_data = []
# y_data = []
# z_data = []
#
# data = [[400, Average_Shear_Stress_List_400K_20ms[0], LogRate_20ms_400K_1GPa], [400, Average_Shear_Stress_List_400K_20ms[1], LogRate_20ms_400K_2GPa], [400, Average_Shear_Stress_List_400K_20ms[2], LogRate_20ms_400K_3GPa], [400, Average_Shear_Stress_List_400K_20ms[3], LogRate_20ms_400K_4GPa], [400, Average_Shear_Stress_List_400K_20ms[4], LogRate_20ms_400K_5GPa],
#         [500, Average_Shear_Stress_List_500K_20ms[0], LogRate_20ms_500K_1GPa], [500, Average_Shear_Stress_List_500K_20ms[1], LogRate_20ms_500K_2GPa], [500, Average_Shear_Stress_List_500K_20ms[2], LogRate_20ms_500K_3GPa], [500, Average_Shear_Stress_List_500K_20ms[3], LogRate_20ms_500K_4GPa], [500, Average_Shear_Stress_List_500K_20ms[4], LogRate_20ms_500K_5GPa],
#         [600, Average_Shear_Stress_List_600K_20ms[0], LogRate_20ms_600K_1GPa], [600, Average_Shear_Stress_List_600K_20ms[1], LogRate_20ms_600K_2GPa], [600, Average_Shear_Stress_List_600K_20ms[2], LogRate_20ms_600K_3GPa], [600, Average_Shear_Stress_List_600K_20ms[3], LogRate_20ms_600K_4GPa], [600, Average_Shear_Stress_List_600K_20ms[4], LogRate_20ms_600K_5GPa],
#         [700, Average_Shear_Stress_List_700K_20ms[0], LogRate_20ms_700K_1GPa], [700, Average_Shear_Stress_List_700K_20ms[1], LogRate_20ms_700K_2GPa], [700, Average_Shear_Stress_List_700K_20ms[2], LogRate_20ms_700K_3GPa], [700, Average_Shear_Stress_List_700K_20ms[3], LogRate_20ms_700K_4GPa], [700, Average_Shear_Stress_List_700K_20ms[4], LogRate_20ms_700K_5GPa]]
#
# for item in data:
#     x_data.append(item[0])
#     y_data.append(item[1])
#     z_data.append(item[2])
# #
# #
# parameters, covariance = optimize.curve_fit(function, [x_data, y_data], z_data)
#
# # create surface function model
# # setup data points for calculating surface model
# model_x_data = np.linspace(min(x_data), max(x_data), 40)
# model_y_data = np.linspace(min(y_data), max(y_data), 40)
#
# # create coordinate arrays for vectorized evaluations
# X, Y = np.meshgrid(model_x_data, model_y_data)
# # calculate Z coordinate array
# Z = function(np.array([X, Y]), *parameters)
#
# z = []
# for row in Z:
#     row.sort()
#     z.append(row)
#
# zlogs = np.array(z)
# import matplotlib.cm
# cm = plt.get_cmap("jet")
# fig = plt.figure()
# ax4 = plt.axes(projection='3d')
# ax4.plot_surface(X, Y, Z, cmap=cm, alpha=0.5, edgecolor='black', linewidth=0.3)
# ax4.set_title('3D Plot - Variation in Log of Dissociation Rates - Alpha Fe, 20ms')
# ax4.scatter(x_data, y_data, z_data, color='black', alpha=1)
# ax4.set_xlabel('Temperature (K)')
# ax4.invert_xaxis()
# ax4.set_ylabel('Shear Stress (GPa)')
# ax4.set_zlabel('Log of Dissociation Rate (per ns)')
# plt.show()
# plt.close(fig)
#
# x_data = []
# y_data = []
# z_data = []
#
# data = [[400, Average_Shear_Stress_List_400K_20ms[0], Dissociation_Rate_20ms_400K_1GPa], [400, Average_Shear_Stress_List_400K_20ms[1], Dissociation_Rate_20ms_400K_2GPa], [400, Average_Shear_Stress_List_400K_20ms[2], Dissociation_Rate_20ms_400K_3GPa], [400, Average_Shear_Stress_List_400K_20ms[3], Dissociation_Rate_20ms_400K_4GPa], [400, Average_Shear_Stress_List_400K_20ms[4], Dissociation_Rate_20ms_400K_5GPa],
#         [500, Average_Shear_Stress_List_500K_20ms[0], Dissociation_Rate_20ms_500K_1GPa], [500, Average_Shear_Stress_List_500K_20ms[1], Dissociation_Rate_20ms_500K_2GPa], [500, Average_Shear_Stress_List_500K_20ms[2], Dissociation_Rate_20ms_500K_3GPa], [500, Average_Shear_Stress_List_500K_20ms[3], Dissociation_Rate_20ms_500K_4GPa], [500, Average_Shear_Stress_List_500K_20ms[4], Dissociation_Rate_20ms_500K_5GPa],
#         [600, Average_Shear_Stress_List_600K_20ms[0], Dissociation_Rate_20ms_600K_1GPa], [600, Average_Shear_Stress_List_600K_20ms[1], Dissociation_Rate_20ms_600K_2GPa], [600, Average_Shear_Stress_List_600K_20ms[2], Dissociation_Rate_20ms_600K_3GPa], [600, Average_Shear_Stress_List_600K_20ms[3], Dissociation_Rate_20ms_600K_4GPa], [600, Average_Shear_Stress_List_600K_20ms[4], Dissociation_Rate_20ms_600K_5GPa],
#         [700, Average_Shear_Stress_List_700K_20ms[0], Dissociation_Rate_20ms_700K_1GPa], [700, Average_Shear_Stress_List_700K_20ms[1], Dissociation_Rate_20ms_700K_2GPa], [700, Average_Shear_Stress_List_700K_20ms[2], Dissociation_Rate_20ms_700K_3GPa], [700, Average_Shear_Stress_List_700K_20ms[3], Dissociation_Rate_20ms_700K_4GPa], [700, Average_Shear_Stress_List_700K_20ms[4], Dissociation_Rate_20ms_700K_5GPa]]
#
# for item in data:
#     x_data.append(item[0])
#     y_data.append(item[1])
#     z_data.append(item[2])
#
# parameters, covariance = optimize.curve_fit(function, [x_data, y_data], z_data)
#
# # create surface function model
# # setup data points for calculating surface model
# model_x_data = np.linspace(min(x_data), max(x_data), 40)
# model_y_data = np.linspace(min(y_data), max(y_data), 40)
# # create coordinate arrays for vectorized evaluations
# X, Y = np.meshgrid(model_x_data, model_y_data)
# # calculate Z coordinate array
# Z = function(np.array([X, Y]), *parameters)
#
# z = []
# for row in Z:
#     row.sort()
#     z.append(row)
#
# z = np.array(z)
# # print(Z)
# # print('####################')
# # print(z)
#
# cm = plt.get_cmap("jet")
# fig = plt.figure()
# ax4 = plt.axes(projection='3d')
# ax4.plot_surface(X, Y, Z, cmap=cm, alpha=0.5, edgecolor= 'black', linewidth=0.3)
# ax4.set_title('3D Plot - Variation in Dissociation Rates - Alpha Fe, 20ms')
# ax4.scatter(x_data, y_data, z_data, color='black', alpha=1)
# ax4.set_xlabel('Temperature (K)')
# ax4.invert_xaxis()
# ax4.set_ylabel('Shear Stress (GPa)')
# ax4.set_zlabel('Dissociation Rate (per ns)')
# plt.show()
#
# # ########## CALCULATING THE ACTIVATION ENERGY, VOLUME AND PREFACTOR FROM 3D FIT ##########
# #
# logRates3D_20ms = zlogs
# shearstresses20ms = model_y_data
# Index = 0
# sigma = []
# params = []
# while Index < len(shearstresses20ms) - 1:
#     for logRates3Drow in logRates3D_20ms:
#         coef_sh_cf, coef_sh_pcov = optimize.curve_fit(linear, shearstresses20ms, logRates3Drow)
#         paramsrow = np.polyfit(shearstresses20ms, logRates3Drow, 1)
#         sigmarow = coef_sh_pcov[0, 0] ** 0.5
#         sigma.append(sigmarow)
#         params.append(paramsrow[0])
#         Index +=1
#
# #
# sigma = np.array(sigma)
# params = np.array(params)
# params = np.average(params)
# sigma = np.average(sigma)
#
# activation_vol_3D_20ms = (params) * (1.38065) * 500 * 1e-2
# print(f'3D Activation Volume, 20ms is {activation_vol_3D_20ms}')
#
# alpha = 0.05
# sigma = sigma ** 0.5
# dof = 2
# tval = t.ppf(1.0 - alpha / 2., dof)
# uncert = sigma * tval
# uncert_ActivationVolume_20ms = uncert * (1.38065) * 500 * 1e-2
# print(f'Activation Volume uncertainty for 3D fit at 20ms is {uncert_ActivationVolume_20ms}')
#
# logRates3D_20ms = zlogs
# temperatures = model_x_data
#
# inverse_temperatures = [1 / x for x in temperatures]
#
# Index = 0
# sigma = []
# SigmaA = []
# params = []
# interceptaverage = []
# while Index < len(inverse_temperatures) - 1:
#     for logRates3Drow in logRates3D_20ms:
#         coef_sh_cf, coef_sh_pcov = optimize.curve_fit(linear, inverse_temperatures, logRates3Drow)
#         paramsrow = np.polyfit(inverse_temperatures, logRates3Drow, 1)
#         sigmarow = coef_sh_pcov[0, 0] ** 0.5
#         sigma_A = coef_sh_pcov[1, 1] ** 0.5
#         sigma.append(sigmarow)
#         SigmaA.append(sigma_A)
#         params.append(paramsrow[0])
#         intercept = paramsrow[1]
#         interceptaverage.append((intercept))
#         Index +=1
#
# sigma = np.array(sigma)
# params = np.array(params)
# SigmaA = np.array(SigmaA)
# interceptaverage = np.array(interceptaverage)
# params_20ms = np.average(params)
# sigma = np.average(sigma)
# interceptaverage = np.average(interceptaverage)
# alpha = 0.05
# sigma = sigma ** 0.5
# tval = t.ppf(1.0 - alpha / 2., dof)
# uncert = sigma * tval
#
# Mu1GPa_20ms = np.average(Average_Mu_List_1GPa_20ms)
# Mu2GPa_20ms = np.average(Average_Mu_List_2GPa_20ms)
# Mu3GPa_20ms = np.average(Average_Mu_List_3GPa_20ms)
# Mu4GPa_20ms = np.average(Average_Mu_List_4GPa_20ms)
# Mu5GPa_20ms = np.average(Average_Mu_List_5GPa_20ms)
#
# MuAveragesDifferentPressures = np.array([Mu1GPa_20ms, Mu2GPa_20ms, Mu3GPa_20ms, Mu4GPa_20ms, Mu5GPa_20ms])
# AverageMu_20ms = np.average(MuAveragesDifferentPressures)
#
# ActivationEnergy_3D_20ms = (((3 * 1e9 * (AverageMu_20ms * activation_vol_3D_20ms) * 1e-30) - 1.381 * params_20ms * 1e-23) * 6.02214076 * (10**23)) / 1000
#
# print(f'Activation Energy for 3D fit at 20ms is {ActivationEnergy_3D_20ms}')
#
# error_3D_20ms = (np.sqrt((3 * 1e9 * np.average(AverageMu_20ms) * uncert_ActivationVolume_20ms * 1e-30) ** 2 + (1.381 * uncert * 1e-23) ** 2) * 6.02214076 * (10**23)) / 1000
# print(f"Activation_Energy Error at 20ms is {error_3D_20ms}")
#
# uncert_prefactor_3D_20ms = sigma_A * tval
# lnA_3D_20ms = np.log((np.exp(interceptaverage) * (10 ** 9)))
#
# print(f"ln(A) for 3D fit is {lnA_3D_20ms}")
# print(f"ln(A) uncertainty is {uncert_prefactor_3D_20ms}" + str(10 ** 9))
#
# ######### Plotting 3D fit vs 2D fit results along with their error margins ##########
# """
# Need to get a list with:
# - Values for Activation Volumes at different temperatures
# - Values for Activation Energies at different pressures
# - Values for ln(A) at different pressures
# - List of errors for each of the above quantities
# - Make a graph for each surface chemistry
# - Eventually will need to do the same for the different sliding speeds
#
# """
# Temperatures = [400, 500, 600, 700]
# Pressures = [1, 2, 3, 4, 5]
#
#
# Activation_Energies_20ms = [ActivationEnergy_1GPa_20ms, ActivationEnergy_2GPa_20ms, ActivationEnergy_3GPa_20ms, ActivationEnergy_4GPa_20ms, ActivationEnergy_5GPa_20ms]
# Activation_Energy_Errors_20ms = [err_E_1GPa_20ms, err_E_2GPa_20ms, err_E_3GPa_20ms, err_E_4GPa_20ms, err_E_5GPa_20ms]
# ActivationEnergy_3D_20ms = ActivationEnergy_3D_20ms
# ActivationEnergy_3D_error_20ms = error_3D_20ms
# ActivationEnergy_3D_error_UpperBound_Value_20ms = float(float(ActivationEnergy_3D_20ms) + float(ActivationEnergy_3D_error_20ms))
# ActivationEnergy_3D_error_LowerBound_Value_20ms = float(float(ActivationEnergy_3D_20ms) - float(ActivationEnergy_3D_error_20ms))
#
# Activation_Energy_Error_Plot, Ea2Dvs3d  = plt.subplots()
# Ea2Dvs3d.set_title('Comparison of Activation Energies from 2D and 3D Fits - AlphaFe, 20ms')
# Ea2Dvs3d.set_xlabel('Normal Stress (GPa)')
# Ea2Dvs3d.set_ylabel('Activation Energy')
# Ea2Dvs3d.scatter(Pressures, Activation_Energies_20ms)
# Ea2Dvs3d.errorbar(Pressures, Activation_Energies_20ms, yerr=Activation_Energy_Errors_20ms, linestyle="None", fmt='o', capsize=3)
# Ea2Dvs3d.axhline(y=ActivationEnergy_3D_20ms)
# Pressures = [0.5, 1, 2, 3, 4, 5, 5.5]
# Ea2Dvs3d.fill_between(Pressures, ActivationEnergy_3D_error_LowerBound_Value_20ms, ActivationEnergy_3D_error_UpperBound_Value_20ms, alpha=0.4)
# Ea2Dvs3d.set_xlim(0.5, 5.5)
# Ea2Dvs3d.set_ylim(0, 35)
# plt.show()
#
# ################ Activation Volume Errors #############################
#
# Activation_Volumes_20ms = [activation_vol_400K_20ms, activation_vol_500K_20ms, activation_vol_600K_20ms, activation_vol_700K_20ms]
# Activation_Volume_Errors_20ms = [uncert400_20ms, uncert500_20ms, uncert600_20ms, uncert700_20ms]
# Activation_Volume_3D_20ms = activation_vol_3D_20ms
# Activation_Volume_3D_Error_20ms = uncert_ActivationVolume_20ms
# ActivationVolume_3D_error_UpperBound_Value_20ms = float(float(Activation_Volume_3D_20ms) + float(Activation_Volume_3D_Error_20ms))
# ActivationVolume_3D_error_LowerBound_Value_20ms = float(float(Activation_Volume_3D_20ms) - float(Activation_Volume_3D_Error_20ms))
#
# Activation_Volume_Error_Plot, Av2Dvs3d  = plt.subplots()
# Av2Dvs3d.set_title('Comparison of Activation Volumes from 2D and 3D Fits - AlphaFe, 20ms')
# Av2Dvs3d.set_xlabel('Normal Stress(GPa)')
# Av2Dvs3d.set_ylabel('Activation Volume')
# Av2Dvs3d.scatter(Temperatures, Activation_Volumes_20ms)
# Av2Dvs3d.errorbar(Temperatures, Activation_Volumes_20ms, yerr=Activation_Volume_Errors_20ms, linestyle="None", fmt='o', capsize=3)
# Av2Dvs3d.axhline(y=Activation_Volume_3D_20ms)
# Temperatures = [350, 400, 500, 600, 700, 750]
# Av2Dvs3d.fill_between(Temperatures, ActivationVolume_3D_error_LowerBound_Value_20ms, ActivationVolume_3D_error_UpperBound_Value_20ms, alpha=0.4)
# Av2Dvs3d.set_xlim(350, 750)
# Av2Dvs3d.set_ylim(0, 30)
# plt.show()
#
# #################### Prefactor Errors #########################
# Pressures = [1, 2, 3, 4, 5]
# Prefactors_20ms = [lnA1GPa_20ms, lnA2GPa_20ms, lnA3GPa_20ms, lnA4GPa_20ms, lnA5GPa_20ms]
# PrefactorErrors_20ms = [uncert_A1GPa_20ms, uncert_A2GPa_20ms, uncert_A3GPa_20ms, uncert_A4GPa_20ms, uncert_A5GPa_20ms]
# lnA_3D_error_20ms = uncert_prefactor_3D_20ms
# lnA_3D_error_UpperBound_Value_20ms = float(float(lnA_3D_20ms) + float(lnA_3D_error_20ms))
# lnA_3D_error_LowerBound_Value_20ms = float(float(lnA_3D_20ms) - float(lnA_3D_error_20ms))
#
# Prefactor_Error_Plot, lnA2Dvs3d  = plt.subplots()
# lnA2Dvs3d.set_title('Comparison of Prefactors from 2D and 3D Fits - AlphaFe, 20ms')
# lnA2Dvs3d.set_xlabel('Normal Stress(GPa)')
# lnA2Dvs3d.set_ylabel('Prefactor')
# lnA2Dvs3d.scatter(Pressures, Prefactors_20ms)
# lnA2Dvs3d.errorbar(Pressures, Prefactors_20ms, yerr=PrefactorErrors_20ms, linestyle="None", fmt='o', capsize=3)
# lnA2Dvs3d.axhline(y=lnA_3D_20ms)
# Pressures = [0.5, 1, 2, 3, 4, 5, 5.5]
# lnA2Dvs3d.fill_between(Pressures, lnA_3D_error_LowerBound_Value_20ms, lnA_3D_error_UpperBound_Value_20ms, alpha=0.4)
# lnA2Dvs3d.set_xlim(0.5, 5.5)
# lnA2Dvs3d.set_ylim(20, 30)
# plt.show()
#
# ########## Checking for Kinetic Compensation Effect ##########
#
# KineticCompeEffectPlot, EavsLnA  = plt.subplots()
# EavsLnA.set_title('Activation Energy vs Prefactor - AlphaFe, 20ms')
# EavsLnA.set_xlabel('Ea')
# EavsLnA.set_ylabel('ln A')
# EavsLnA.scatter(Activation_Energies_20ms, Prefactors_20ms)
# EavsLnA.set_ylim(24, 27)
# plt.show()
#
# ################### 30ms ####################################
#
# Index = 0
# Temperatures = ["500K", "600K", "700K"]
# Pressures = ['2GPa', '3GPa', '4GPa', '5GPa']
# Speeds = ['20ms', '30ms', '40ms', '50ms']
#
# Dissociation_Rate_30ms_400K_1GPa, LogRate_30ms_400K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_400K_1GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_30ms_400K_2GPa, LogRate_30ms_400K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_400K_2GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_30ms_400K_3GPa, LogRate_30ms_400K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_400K_3GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_30ms_400K_4GPa, LogRate_30ms_400K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_400K_4GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_30ms_400K_5GPa, LogRate_30ms_400K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_400K_5GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rates_30ms_400K = [Dissociation_Rate_30ms_400K_1GPa, Dissociation_Rate_30ms_400K_2GPa, Dissociation_Rate_30ms_400K_3GPa, Dissociation_Rate_30ms_400K_4GPa, Dissociation_Rate_30ms_400K_5GPa]
# Log_Rates_400K_30ms = [LogRate_30ms_400K_1GPa, LogRate_30ms_400K_2GPa, LogRate_30ms_400K_3GPa, LogRate_30ms_400K_4GPa, LogRate_30ms_400K_5GPa]
#
# Dissociation_Rate_30ms_500K_1GPa, LogRate_30ms_500K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_500K_1GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_30ms_500K_2GPa, LogRate_30ms_500K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_500K_2GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_30ms_500K_3GPa, LogRate_30ms_500K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_500K_3GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_30ms_500K_4GPa, LogRate_30ms_500K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_500K_4GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_30ms_500K_5GPa, LogRate_30ms_500K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_500K_5GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rates_30ms_500K = [Dissociation_Rate_30ms_500K_1GPa, Dissociation_Rate_30ms_500K_2GPa, Dissociation_Rate_30ms_500K_3GPa, Dissociation_Rate_30ms_500K_4GPa, Dissociation_Rate_30ms_500K_5GPa]
# Log_Rates_500K_30ms = [LogRate_30ms_500K_1GPa, LogRate_30ms_500K_2GPa, LogRate_30ms_500K_3GPa, LogRate_30ms_500K_4GPa, LogRate_30ms_500K_5GPa]
#
# Dissociation_Rate_30ms_600K_1GPa, LogRate_30ms_600K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_600K_1GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_30ms_600K_2GPa, LogRate_30ms_600K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_600K_2GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_30ms_600K_3GPa, LogRate_30ms_600K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_600K_3GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_30ms_600K_4GPa, LogRate_30ms_600K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_600K_4GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_30ms_600K_5GPa, LogRate_30ms_600K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_600K_5GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rates_30ms_600K = [Dissociation_Rate_30ms_600K_1GPa, Dissociation_Rate_30ms_600K_2GPa, Dissociation_Rate_30ms_600K_3GPa, Dissociation_Rate_30ms_600K_4GPa, Dissociation_Rate_30ms_600K_5GPa]
# Log_Rates_600K_30ms = [LogRate_30ms_600K_1GPa, LogRate_30ms_600K_2GPa, LogRate_30ms_600K_3GPa, LogRate_30ms_600K_4GPa, LogRate_30ms_600K_5GPa]
#
# Dissociation_Rate_30ms_700K_1GPa, LogRate_30ms_700K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_700K_1GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_30ms_700K_2GPa, LogRate_30ms_700K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_700K_2GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_30ms_700K_3GPa, LogRate_30ms_700K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_700K_3GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_30ms_700K_4GPa, LogRate_30ms_700K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_700K_4GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_30ms_700K_5GPa, LogRate_30ms_700K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_30ms_700K_5GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rates_30ms_700K = [Dissociation_Rate_30ms_700K_1GPa, Dissociation_Rate_30ms_700K_2GPa, Dissociation_Rate_30ms_700K_3GPa, Dissociation_Rate_30ms_700K_4GPa, Dissociation_Rate_30ms_700K_5GPa]
# Log_Rates_700K_30ms = [LogRate_30ms_700K_1GPa, LogRate_30ms_700K_2GPa, LogRate_30ms_700K_3GPa, LogRate_30ms_700K_4GPa, LogRate_30ms_700K_5GPa]
#
# Dissociation_Rates_1GPa = [Dissociation_Rates_30ms_400K[0], Dissociation_Rates_30ms_500K[0], Dissociation_Rates_30ms_600K[0], Dissociation_Rates_30ms_700K[0]]
# Dissociation_Rates_2GPa = [Dissociation_Rates_30ms_400K[1], Dissociation_Rates_30ms_500K[1], Dissociation_Rates_30ms_600K[1], Dissociation_Rates_30ms_700K[1]]
# Dissociation_Rates_3GPa = [Dissociation_Rates_30ms_400K[2], Dissociation_Rates_30ms_500K[2], Dissociation_Rates_30ms_600K[2], Dissociation_Rates_30ms_700K[2]]
# Dissociation_Rates_4GPa = [Dissociation_Rates_30ms_400K[3], Dissociation_Rates_30ms_500K[3], Dissociation_Rates_30ms_600K[3], Dissociation_Rates_30ms_700K[3]]
# Dissociation_Rates_5GPa = [Dissociation_Rates_30ms_400K[4], Dissociation_Rates_30ms_500K[4], Dissociation_Rates_30ms_600K[4], Dissociation_Rates_30ms_700K[4]]
#
# Log_Rates_1GPa_30ms = [LogRate_30ms_400K_1GPa, LogRate_30ms_500K_1GPa, LogRate_30ms_600K_1GPa, LogRate_30ms_700K_1GPa]
# Log_Rates_2GPa_30ms = [LogRate_30ms_400K_2GPa, LogRate_30ms_500K_2GPa, LogRate_30ms_600K_2GPa, LogRate_30ms_700K_2GPa]
# Log_Rates_3GPa_30ms = [LogRate_30ms_400K_3GPa, LogRate_30ms_500K_3GPa, LogRate_30ms_600K_3GPa, LogRate_30ms_700K_3GPa]
# Log_Rates_4GPa_30ms = [LogRate_30ms_400K_4GPa, LogRate_30ms_500K_4GPa, LogRate_30ms_600K_4GPa, LogRate_30ms_700K_4GPa]
# Log_Rates_5GPa_30ms = [LogRate_30ms_400K_5GPa, LogRate_30ms_500K_5GPa, LogRate_30ms_600K_5GPa, LogRate_30ms_700K_5GPa]
#
# EquilibriumFactor = [79, 99] # How many rows (out of 99) to ignore before calculating shear stress/friction coefficient, as it won't stabilise until after a certain number of timesteps
#
# def get_average_shear_normal_stress_and_average_mu_constant_temperature(Temperature, Pressures, EquilibriumFactor, Speed):
#     Friction_Coefficient_Dataframe_Unnamed = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{Speed}/{Temperature}/1GPa/'
#                                 'fc_ave.dump'.format(Speed=Speed, Temperature=Temperature), sep=' ')
#     Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe_Unnamed.rename(columns={'v_s_bot' : 'Shear Stress 1GPa', 'v_p_bot' : 'Normal Stress 1GPa'})
#
#     for P in Pressures:
#         Dataframe = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{Speed}/{Temperature}/{P}/'
#                                 'fc_ave.dump'.format(Speed=Speed, Temperature=Temperature, P=P), sep=' ')
#         Big_DataframeP = Dataframe.rename(columns= {'Timestep': 'Timestep {}'.format(P),
#                                                         'v_s_bot': 'Shear Stress {}'.format(P),
#                                                         'v_p_bot': 'Normal Stress {}'.format(P)})
#
#         Friction_Coefficient_Dataframe = pd.concat([Friction_Coefficient_Dataframe, Big_DataframeP], axis =1)
#         Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe.dropna()
#
#
#     #print(Friction_Coefficient_Dataframe)
#     Mu_Final_Dataframe = Friction_Coefficient_Dataframe.iloc[:, [0, 1, 2, 4, 5, 7, 8, 10, 11, 13, 14]]
#     Mu_Final_Dataframe = Mu_Final_Dataframe.iloc[EquilibriumFactor[0]:EquilibriumFactor[1], :]
#     #print(Mu_Final_Dataframe)
#
#     ShearStressMeans = Mu_Final_Dataframe[['Shear Stress 1GPa', 'Shear Stress 2GPa', 'Shear Stress 3GPa', 'Shear Stress 4GPa', 'Shear Stress 5GPa']].mean()
#     Average_Shear_Stress_Dictionary = ShearStressMeans.to_dict()
#     #print(ShearStressMeans)
#     NormalStressMeans = Mu_Final_Dataframe[['Normal Stress 1GPa', 'Normal Stress 2GPa', 'Normal Stress 3GPa', 'Normal Stress 4GPa', 'Normal Stress 5GPa']].mean()
#     NormalStressMeans = NormalStressMeans.to_dict()
#     #print(NormalStressMeans)
#
#     Average_Mu_Dictionary = {}
#
#     NormalStressMeansList = list(NormalStressMeans.values())
#     Normal_Stress = NormalStressMeans.get('Normal Stress 1GPa')
#     Shear_Stress = ShearStressMeans.get('Shear Stress 1GPa')
#     Average_Mu = Shear_Stress / Normal_Stress
#     Average_Mu_Dictionary.update({'Average Mu 1GPa': Average_Mu})
#
#     for P in Pressures:
#
#         Normal_Stress = NormalStressMeans.get('Normal Stress {}'.format(P))
#         Shear_Stress = ShearStressMeans.get('Shear Stress {}'.format(P))
#         Average_Mu = Shear_Stress / Normal_Stress
#         Average_Mu_Dictionary.update({'Average Mu {}'.format(P): Average_Mu})
#
#
#     Average_Shear_Stress_List = list(Average_Shear_Stress_Dictionary.values())
#     #print(Average_Shear_Stress_List)
#     Average_Mu_List = list(Average_Mu_Dictionary.values())
#     Average_Shear_Stress_List = [x / 10000 for x in Average_Shear_Stress_List] # Conversion to GPa
#     NormalStressMeansList = [x/10000 for x in NormalStressMeansList]
#     #print(Average_Shear_Stress_List)
#
#     return Average_Shear_Stress_List, Average_Mu_List, NormalStressMeansList
#
# ######## Getting Average Shear Stress, Friction Coefficient and Normal Stress #################
# Average_Shear_Stress_List_400K_30ms, Average_Mu_List_400K_30ms, NormalStressMeans_400K_30ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("400K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="30ms")
# Average_Shear_Stress_List_500K_30ms, Average_Mu_List_500K_30ms, NormalStressMeans_500K_30ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("500K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="30ms")
# Average_Shear_Stress_List_600K_30ms, Average_Mu_List_600K_30ms, NormalStressMeans_600K_30ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("600K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="30ms")
# Average_Shear_Stress_List_700K_30ms, Average_Mu_List_700K_30ms, NormalStressMeans_700K_30ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("700K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="30ms")
#
# Average_Mu_List_1GPa_30ms = [Average_Mu_List_400K_30ms[0], Average_Mu_List_500K_30ms[0], Average_Mu_List_600K_30ms[0], Average_Mu_List_700K_30ms[0]]
# Average_Mu_List_2GPa_30ms = [Average_Mu_List_400K_30ms[1], Average_Mu_List_500K_30ms[1], Average_Mu_List_600K_30ms[1], Average_Mu_List_700K_30ms[1]]
# Average_Mu_List_3GPa_30ms = [Average_Mu_List_400K_30ms[2], Average_Mu_List_500K_30ms[2], Average_Mu_List_600K_30ms[2], Average_Mu_List_700K_30ms[2]]
# Average_Mu_List_4GPa_30ms = [Average_Mu_List_400K_30ms[3], Average_Mu_List_500K_30ms[3], Average_Mu_List_600K_30ms[3], Average_Mu_List_700K_30ms[3]]
# Average_Mu_List_5GPa_30ms = [Average_Mu_List_400K_30ms[4], Average_Mu_List_500K_30ms[4], Average_Mu_List_600K_30ms[4], Average_Mu_List_700K_30ms[4]]
#
# Average_Shear_Stress_List_1GPa_30ms = [Average_Shear_Stress_List_400K_30ms[0], Average_Shear_Stress_List_500K_30ms[0], Average_Shear_Stress_List_600K_30ms[0], Average_Shear_Stress_List_700K_30ms[0]]
# Average_Shear_Stress_List_2GPa_30ms = [Average_Shear_Stress_List_400K_30ms[1], Average_Shear_Stress_List_500K_30ms[1], Average_Shear_Stress_List_600K_30ms[1], Average_Shear_Stress_List_700K_30ms[1]]
# Average_Shear_Stress_List_3GPa_30ms = [Average_Shear_Stress_List_400K_30ms[2], Average_Shear_Stress_List_500K_30ms[2], Average_Shear_Stress_List_600K_30ms[2], Average_Shear_Stress_List_700K_30ms[2]]
# Average_Shear_Stress_List_4GPa_30ms = [Average_Shear_Stress_List_400K_30ms[3], Average_Shear_Stress_List_500K_30ms[3], Average_Shear_Stress_List_600K_30ms[3], Average_Shear_Stress_List_700K_30ms[3]]
# Average_Shear_Stress_List_5GPa_30ms = [Average_Shear_Stress_List_400K_30ms[4], Average_Shear_Stress_List_500K_30ms[4], Average_Shear_Stress_List_600K_30ms[4], Average_Shear_Stress_List_700K_30ms[4]]
#
# plot_shear_stress_vs_normal_stress(Average_Shear_Stress_List_400K_30ms, Average_Shear_Stress_List_500K_30ms, Average_Shear_Stress_List_600K_30ms, Average_Shear_Stress_List_700K_30ms,
#                                  "400K", "500K", "600K", "700K", Speed="30ms")
#
# ########################## Plotting ln(Rates) vs Shear Stress #####################################
# x = np.array([0, 1, 2, 3, 4, 5])
# params400K_30ms = np.polyfit(Average_Shear_Stress_List_400K_30ms, Log_Rates_400K_30ms, 1)
# params500K_30ms = np.polyfit(Average_Shear_Stress_List_500K_30ms, Log_Rates_500K_30ms, 1)
# params600K_30ms = np.polyfit(Average_Shear_Stress_List_600K_30ms, Log_Rates_600K_30ms, 1)
# params700K_30ms = np.polyfit(Average_Shear_Stress_List_700K_30ms, Log_Rates_700K_30ms, 1)
#
# RatesvsShear, RvS3 = plt.subplots()
# RvS3.set_title('Log of Dissociation Rates vs Shear Stress, 30ms')
# RvS3.set_xlabel('Shear Stress(GPa)')
# RvS3.set_ylabel('Log of Dissociation Rate (per nanosecond)')
# RvS3.scatter(Average_Shear_Stress_List_400K_30ms, Log_Rates_400K_30ms)
# RvS3.scatter(Average_Shear_Stress_List_500K_30ms, Log_Rates_500K_30ms)
# RvS3.scatter(Average_Shear_Stress_List_600K_30ms, Log_Rates_600K_30ms)
# RvS3.scatter(Average_Shear_Stress_List_700K_30ms, Log_Rates_700K_30ms)
# RvS3.plot(x, params400K_30ms[0] * x + params400K_30ms[1], label='400K Fitted')
# RvS3.plot(x, params500K_30ms[0] * x + params500K_30ms[1], label='500K Fitted')
# RvS3.plot(x, params600K_30ms[0] * x + params600K_30ms[1], label='600K Fitted')
# RvS3.plot(x, params700K_30ms[0] * x + params700K_30ms[1], label='700K Fitted')
# RvS3.set_xlim(0, 2)
# RvS3.set_ylim(-0.5, 4)
# RvS3.legend(loc='lower right')
# plt.show()
#
# ####### Calculate Activation Volume, Using  Carlos' conversion to get in Angstrom^3 #################
#
# activation_vol_400K_30ms = (params400K_30ms[0]) * (1.38065) * 400 * 1e-2
# activation_vol_500K_30ms = (params500K_30ms[0]) * (1.38065) * 500 * 1e-2
# activation_vol_600K_30ms = (params600K_30ms[0]) * (1.38065) * 600 * 1e-2
# activation_vol_700K_30ms = (params700K_30ms[0]) * (1.38065) * 700 * 1e-2
#
# alpha = 0.05
# def linear(x, m, n):
#     return m * x + n
#
# coef_sh_cf, coef_sh_pcov = optimize.curve_fit(linear, Average_Shear_Stress_List_700K_30ms, Log_Rates_700K_30ms)
# sigma = coef_sh_pcov[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_700K_30ms) - len(params700K_30ms))
#
# tval = t.ppf(1.0 - alpha / 2., dof)
# uncert = sigma * tval
# uncert400_30ms = uncert * (1.38065) * 400 * 1e-2
# uncert500_30ms = uncert * (1.38065) * 500 * 1e-2
# uncert600_30ms = uncert * (1.38065) * 600 * 1e-2
# uncert700_30ms = uncert * (1.38065) * 700 * 1e-2
#
# print(f'Activation Volume uncertainty at 400K, 30ms is {uncert400_30ms}')
# print(f'Activation Volume uncertainty at 500K, 30ms is {uncert500_30ms}')
# print(f'Activation Volume uncertainty at 600K, 30ms is {uncert600_30ms}')
# print(f'Activation Volume uncertainty at 700K, 30ms is {uncert700_30ms}')
#
# Activation_Volumes_30ms = [activation_vol_400K_30ms, activation_vol_500K_30ms, activation_vol_600K_30ms, activation_vol_700K_30ms]
# mean_actv_30ms = np.average(Activation_Volumes_30ms)
#
# print('Activation Volume 400K, 30ms = ' + str(activation_vol_400K_30ms))
# print('Activation Volume 500K, 30ms = ' + str(activation_vol_500K_30ms))
# print('Activation Volume 600K, 30ms = ' + str(activation_vol_600K_30ms))
# print('Activation Volume 700K, 30ms = ' + str(activation_vol_700K_30ms))
# #
# ############ Plotting lnk vs 1000/T  #########################
#
# Temperatures = [400, 500, 600, 700]
# Inverse_Temperatures = np.array([1/x for x in Temperatures])
#
# trend1GPa_30ms = np.polyfit(Inverse_Temperatures, Log_Rates_1GPa_30ms, 1)
# trend2GPa_30ms = np.polyfit(Inverse_Temperatures, Log_Rates_2GPa_30ms, 1)
# trend3GPa_30ms = np.polyfit(Inverse_Temperatures, Log_Rates_3GPa_30ms, 1)
# trend4GPa_30ms = np.polyfit(Inverse_Temperatures, Log_Rates_4GPa_30ms, 1)
# trend5GPa_30ms = np.polyfit(Inverse_Temperatures, Log_Rates_5GPa_30ms, 1)
#
# fig1, ax1 = plt.subplots()
# ax1.set_title('Log of Dissociation Rates against Inverse of Temperatures, 30ms')
# ax1.set_xlabel('1000/T (K-1)')
# ax1.set_ylabel('ln(Rate) (ns-1)')
# ax1.scatter(Inverse_Temperatures, Log_Rates_1GPa_30ms)
# ax1.scatter(Inverse_Temperatures, Log_Rates_2GPa_30ms)
# ax1.scatter(Inverse_Temperatures, Log_Rates_3GPa_30ms)
# ax1.scatter(Inverse_Temperatures, Log_Rates_4GPa_30ms)
# ax1.scatter(Inverse_Temperatures, Log_Rates_5GPa_30ms)
#
# Fit1GPa_30ms = np.poly1d(trend1GPa_30ms)
# Fit2GPa_30ms = np.poly1d(trend2GPa_30ms)
# Fit3GPa_30ms = np.poly1d(trend3GPa_30ms)
# Fit4GPa_30ms = np.poly1d(trend4GPa_30ms)
# Fit5GPa_30ms = np.poly1d(trend5GPa_30ms)
#
# ax1.plot(Inverse_Temperatures, Fit1GPa_30ms(Inverse_Temperatures), label='1GPa')
# ax1.plot(Inverse_Temperatures, Fit2GPa_30ms(Inverse_Temperatures), label='2GPa')
# ax1.plot(Inverse_Temperatures, Fit3GPa_30ms(Inverse_Temperatures), label='3GPa')
# ax1.plot(Inverse_Temperatures, Fit4GPa_30ms(Inverse_Temperatures), label='4GPa')
# ax1.plot(Inverse_Temperatures, Fit5GPa_30ms(Inverse_Temperatures), label='5GPa')
# ax1.legend()
# plt.show()
#
# ActivationEnergy_1GPa_30ms = (((1 * 1e9 * (np.average(Average_Mu_List_1GPa_30ms) * mean_actv_30ms) * 1e-30) - 1.381 * trend1GPa_30ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
# ActivationEnergy_2GPa_30ms = (((2 * 1e9 * (np.average(Average_Mu_List_2GPa_30ms) * mean_actv_30ms) * 1e-30) - 1.381 * trend2GPa_30ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
# ActivationEnergy_3GPa_30ms = (((3 * 1e9 * (np.average(Average_Mu_List_3GPa_30ms) * mean_actv_30ms) * 1e-30) - 1.381 * trend3GPa_30ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
# ActivationEnergy_4GPa_30ms = (((4 * 1e9 * (np.average(Average_Mu_List_4GPa_30ms) * mean_actv_30ms) * 1e-30) - 1.381 * trend4GPa_30ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
# ActivationEnergy_5GPa_30ms = (((5 * 1e9 * (np.average(Average_Mu_List_5GPa_30ms) * mean_actv_30ms) * 1e-30) - 1.381 * trend5GPa_30ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
#
# print(np.average(Average_Mu_List_1GPa_30ms))
# print(np.average(Average_Mu_List_2GPa_30ms))
# print(np.average(Average_Mu_List_3GPa_30ms))
# print(np.average(Average_Mu_List_4GPa_30ms))
# print(np.average(Average_Mu_List_5GPa_30ms))
#
# lnA1GPa_30ms = np.log((np.exp(trend1GPa_30ms[1]) * (10 ** 9)))
# lnA2GPa_30ms = np.log((np.exp(trend2GPa_30ms[1]) * (10 ** 9)))
# lnA3GPa_30ms = np.log((np.exp(trend3GPa_30ms[1]) * (10 ** 9)))
# lnA4GPa_30ms = np.log((np.exp(trend4GPa_30ms[1]) * (10 ** 9)))
# lnA5GPa_30ms = np.log((np.exp(trend5GPa_30ms[1]) * (10 ** 9)))
#
# print(f"ln(A) at 1GPa, 30ms is {lnA1GPa_30ms}")
# print(f"ln(A) at 2GPa, 30ms is {lnA2GPa_30ms}")
# print(f"ln(A) at 3GPa, 30ms is {lnA3GPa_30ms}")
# print(f"ln(A) at 4GPa, 30ms is {lnA4GPa_30ms}")
# print(f"ln(A) at 5GPa, 30ms is {lnA5GPa_30ms}")
#
# print('Activation Energy 1GPa, 30ms =' + str(ActivationEnergy_1GPa_30ms))
# print('Activation Energy 2GPa, 30ms =' + str(ActivationEnergy_2GPa_30ms))
# print('Activation Energy 3GPa, 30ms =' + str(ActivationEnergy_3GPa_30ms))
# print('Activation Energy 4GPa, 30ms =' + str(ActivationEnergy_4GPa_30ms))
# print('Activation Energy 5GPa, 30ms =' + str(ActivationEnergy_5GPa_30ms))
#
# def linear(x, m, n):
#     return m * x + n
#
# coef_p_cf_1GPa_30ms, coef_p_pcov_1GPa_30ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_1GPa_30ms)
# sigma_A1GPa_30ms = coef_p_pcov_1GPa_30ms[1, 1] ** 0.5
# sigma_m1GPa_30ms = coef_p_pcov_1GPa_30ms[0, 0] ** 0.5
# dof_30ms = max(0, len(Log_Rates_1GPa_30ms) - len(coef_p_cf_1GPa_30ms))
#
# alpha = 0.05
# tval_30ms = t.ppf(1.0 - alpha / 2., 4.1532)
#
# sigma_30ms = np.std(Activation_Volumes_30ms)
# error_actv_30ms = sigma_30ms * tval_30ms / np.sqrt(len(Activation_Volumes_30ms))
# uncert_A1GPa_30ms = sigma_A1GPa_30ms * tval
# uncert_m1GPa_30ms = sigma_m1GPa_30ms * tval
#
# coef_p_cf_2GPa_30ms, coef_p_pcov_2GPa_30ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_2GPa_30ms)
# sigma_A2GPa_30ms = coef_p_pcov_2GPa_30ms[1, 1] ** 0.5
# sigma_m2GPa_30ms = coef_p_pcov_2GPa_30ms[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_2GPa_30ms) - len(coef_p_cf_2GPa_30ms))
# uncert_A2GPa_30ms = sigma_A2GPa_30ms * tval_30ms
# uncert_m2GPa_30ms = sigma_m2GPa_30ms * tval_30ms
#
# coef_p_cf_3GPa_30ms, coef_p_pcov_3GPa_30ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_3GPa_30ms)
# sigma_A3GPa_30ms = coef_p_pcov_3GPa_30ms[1, 1] ** 0.5
# sigma_m3GPa_30ms = coef_p_pcov_3GPa_30ms[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_3GPa_30ms) - len(coef_p_cf_3GPa_30ms))
# uncert_A3GPa_30ms = sigma_A3GPa_30ms * tval_30ms
# uncert_m3GPa_30ms = sigma_m3GPa_30ms * tval_30ms
#
# coef_p_cf_4GPa_30ms, coef_p_pcov_4GPa_30ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_4GPa_30ms)
# sigma_A4GPa_30ms = coef_p_pcov_4GPa_30ms[1, 1] ** 0.5
# sigma_m4GPa_30ms = coef_p_pcov_4GPa_30ms[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_4GPa_30ms) - len(coef_p_cf_4GPa_30ms))
# uncert_A4GPa_30ms = sigma_A4GPa_30ms * tval_30ms
# uncert_m4GPa_30ms = sigma_m4GPa_30ms * tval_30ms
#
# coef_p_cf_5GPa_30ms, coef_p_pcov_5GPa_30ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_5GPa_30ms)
# sigma_A5GPa_30ms = coef_p_pcov_5GPa_30ms[1, 1] ** 0.5
# sigma_m5GPa_30ms = coef_p_pcov_5GPa_30ms[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_5GPa_30ms) - len(coef_p_cf_5GPa_30ms))
# uncert_A5GPa_30ms = sigma_A5GPa_30ms * tval_30ms
# uncert_m5GPa_30ms = sigma_m5GPa_30ms * tval_30ms
#
# err_E_1GPa_30ms = (np.sqrt(((1 * 1e9 * np.average(Average_Mu_List_1GPa_30ms) * error_actv_30ms * 1e-30) ** 2 + (1.381 * uncert_m1GPa_30ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
# err_E_2GPa_30ms = (np.sqrt(((2 * 1e9 * np.average(Average_Mu_List_2GPa_30ms) * error_actv_30ms * 1e-30) ** 2 + (1.381 * uncert_m2GPa_30ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
# err_E_3GPa_30ms = (np.sqrt(((3 * 1e9 * np.average(Average_Mu_List_3GPa_30ms) * error_actv_30ms * 1e-30) ** 2 + (1.381 * uncert_m3GPa_30ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
# err_E_4GPa_30ms = (np.sqrt(((4 * 1e9 * np.average(Average_Mu_List_4GPa_30ms) * error_actv_30ms * 1e-30) ** 2 + (1.381 * uncert_m4GPa_30ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
# err_E_5GPa_30ms = (np.sqrt(((5 * 1e9 * np.average(Average_Mu_List_5GPa_30ms) * error_actv_30ms * 1e-30) ** 2 + (1.381 * uncert_m5GPa_30ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
#
# print(f"Error for Activation Energy at 1GPa = {err_E_1GPa_30ms}")
# print(f"Error for Activation Energy at 2GPa = {err_E_2GPa_30ms}")
# print(f"Error for Activation Energy at 3GPa = {err_E_3GPa_30ms}")
# print(f"Error for Activation Energy at 4GPa = {err_E_4GPa_30ms}")
# print(f"Error for Activation Energy at 5GPa = {err_E_5GPa_30ms}")
#
# uncertlnA1GPa_30ms = np.log(uncert_A1GPa_30ms * (10 ** 9))
# uncertlnA2GPa_30ms = np.log(uncert_A2GPa_30ms * (10 ** 9))
# uncertlnA3GPa_30ms = np.log(uncert_A3GPa_30ms * (10 ** 9))
# uncertlnA4GPa_30ms = np.log(uncert_A4GPa_30ms * (10 ** 9))
# uncertlnA5GPa_30ms = np.log(uncert_A5GPa_30ms * (10 ** 9))
#
# print(f"Error for ln A at 1GPa = {uncert_A1GPa_30ms}")
# print(f"Error for ln A at 2GPa = {uncert_A2GPa_30ms}")
# print(f"Error for ln A at 3GPa = {uncert_A3GPa_30ms}")
# print(f"Error for ln A at 4GPa = {uncert_A4GPa_30ms}")
# print(f"Error for ln A at 5GPa = {uncert_A5GPa_30ms}")
#
# pressures = ['2GPa', '3GPa', '4GPa', '5GPa']
#
# ########################## Plotting ln(Rates) vs Normal Stress #####################################
# x = np.array([0, 1, 2, 3, 4, 5])
# params400K_30ms = np.polyfit(NormalStressMeans_400K_30ms, Log_Rates_400K_30ms, 1)
# params500K_30ms = np.polyfit(NormalStressMeans_500K_30ms, Log_Rates_500K_30ms, 1)
# params600K_30ms = np.polyfit(NormalStressMeans_600K_30ms, Log_Rates_600K_30ms, 1)
# params700K_30ms = np.polyfit(NormalStressMeans_700K_30ms, Log_Rates_700K_30ms, 1)
#
# RatesvsNormal, RvN3 = plt.subplots()
# RvN3.set_title('Log of Dissociation Rates vs Normal Stress - Alpha Fe, 30ms')
# RvN3.set_xlabel('Normal Stress(GPa)')
# RvN3.set_ylabel('Log of Dissociation Rate (ns-1)')
# RvN3.scatter(NormalStressMeans_400K_30ms, Log_Rates_400K_30ms)
# RvN3.scatter(NormalStressMeans_500K_30ms, Log_Rates_500K_30ms)
# RvN3.scatter(NormalStressMeans_600K_30ms, Log_Rates_600K_30ms)
# RvN3.scatter(NormalStressMeans_700K_30ms, Log_Rates_700K_30ms)
# RvN3.plot(x, params400K_30ms[0] * x + params400K_30ms[1], label='400K Fitted')
# RvN3.plot(x, params500K_30ms[0] * x + params500K_30ms[1], label='500K Fitted')
# RvN3.plot(x, params600K_30ms[0] * x + params600K_30ms[1], label='600K Fitted')
# RvN3.plot(x, params700K_30ms[0] * x + params700K_30ms[1], label='700K Fitted')
# #RvN3.set_xlim(1, 5)
# # RvN3.set_ylim(0, 25)
# RvN3.legend(loc='lower right')
# plt.show()
#
# def function(data, a, b, c):
#     x = data[0]
#     y = data[1]
#     return a * (x**b) * (y**c) #TODO change fitting function
#
# x_data = []
# y_data = []
# z_data = []
#
# data = [[400, Average_Shear_Stress_List_400K_30ms[0], LogRate_30ms_400K_1GPa], [400, Average_Shear_Stress_List_400K_30ms[1], LogRate_30ms_400K_2GPa], [400, Average_Shear_Stress_List_400K_30ms[2], LogRate_30ms_400K_3GPa], [400, Average_Shear_Stress_List_400K_30ms[3], LogRate_30ms_400K_4GPa], [400, Average_Shear_Stress_List_400K_30ms[4], LogRate_30ms_400K_5GPa],
#         [500, Average_Shear_Stress_List_500K_30ms[0], LogRate_30ms_500K_1GPa], [500, Average_Shear_Stress_List_500K_30ms[1], LogRate_30ms_500K_2GPa], [500, Average_Shear_Stress_List_500K_30ms[2], LogRate_30ms_500K_3GPa], [500, Average_Shear_Stress_List_500K_30ms[3], LogRate_30ms_500K_4GPa], [500, Average_Shear_Stress_List_500K_30ms[4], LogRate_30ms_500K_5GPa],
#         [600, Average_Shear_Stress_List_600K_30ms[0], LogRate_30ms_600K_1GPa], [600, Average_Shear_Stress_List_600K_30ms[1], LogRate_30ms_600K_2GPa], [600, Average_Shear_Stress_List_600K_30ms[2], LogRate_30ms_600K_3GPa], [600, Average_Shear_Stress_List_600K_30ms[3], LogRate_30ms_600K_4GPa], [600, Average_Shear_Stress_List_600K_30ms[4], LogRate_30ms_600K_5GPa],
#         [700, Average_Shear_Stress_List_700K_30ms[0], LogRate_30ms_700K_1GPa], [700, Average_Shear_Stress_List_700K_30ms[1], LogRate_30ms_700K_2GPa], [700, Average_Shear_Stress_List_700K_30ms[2], LogRate_30ms_700K_3GPa], [700, Average_Shear_Stress_List_700K_30ms[3], LogRate_30ms_700K_4GPa], [700, Average_Shear_Stress_List_700K_30ms[4], LogRate_30ms_700K_5GPa]]
#
# for item in data:
#     x_data.append(item[0])
#     y_data.append(item[1])
#     z_data.append(item[2])
# #
# #
# parameters, covariance = optimize.curve_fit(function, [x_data, y_data], z_data)
#
# # create surface function model
# # setup data points for calculating surface model
# model_x_data = np.linspace(min(x_data), max(x_data), 40)
# model_y_data = np.linspace(min(y_data), max(y_data), 40)
#
# # create coordinate arrays for vectorized evaluations
# X, Y = np.meshgrid(model_x_data, model_y_data)
# # calculate Z coordinate array
# Z = function(np.array([X, Y]), *parameters)
#
# z = []
# for row in Z:
#     row.sort()
#     z.append(row)
#
# zlogs = np.array(z)
# import matplotlib.cm
# cm = plt.get_cmap("jet")
# fig = plt.figure()
# ax4 = plt.axes(projection='3d')
# ax4.plot_surface(X, Y, Z, cmap=cm, alpha=0.5, edgecolor='black', linewidth=0.3)
# ax4.set_title('3D Plot - Variation in Log of Dissociation Rates - Alpha Fe, 30ms')
# ax4.scatter(x_data, y_data, z_data, color='black', alpha=1)
# ax4.set_xlabel('Temperature (K)')
# ax4.invert_xaxis()
# ax4.set_ylabel('Shear Stress (GPa)')
# ax4.set_zlabel('Log of Dissociation Rate (per ns)')
# plt.show()
# plt.close(fig)
#
# x_data = []
# y_data = []
# z_data = []
#
# data = [[400, Average_Shear_Stress_List_400K_30ms[0], Dissociation_Rate_30ms_400K_1GPa], [400, Average_Shear_Stress_List_400K_30ms[1], Dissociation_Rate_30ms_400K_2GPa], [400, Average_Shear_Stress_List_400K_30ms[2], Dissociation_Rate_30ms_400K_3GPa], [400, Average_Shear_Stress_List_400K_30ms[3], Dissociation_Rate_30ms_400K_4GPa], [400, Average_Shear_Stress_List_400K_30ms[4], Dissociation_Rate_30ms_400K_5GPa],
#         [500, Average_Shear_Stress_List_500K_30ms[0], Dissociation_Rate_30ms_500K_1GPa], [500, Average_Shear_Stress_List_500K_30ms[1], Dissociation_Rate_30ms_500K_2GPa], [500, Average_Shear_Stress_List_500K_30ms[2], Dissociation_Rate_30ms_500K_3GPa], [500, Average_Shear_Stress_List_500K_30ms[3], Dissociation_Rate_30ms_500K_4GPa], [500, Average_Shear_Stress_List_500K_30ms[4], Dissociation_Rate_30ms_500K_5GPa],
#         [600, Average_Shear_Stress_List_600K_30ms[0], Dissociation_Rate_30ms_600K_1GPa], [600, Average_Shear_Stress_List_600K_30ms[1], Dissociation_Rate_30ms_600K_2GPa], [600, Average_Shear_Stress_List_600K_30ms[2], Dissociation_Rate_30ms_600K_3GPa], [600, Average_Shear_Stress_List_600K_30ms[3], Dissociation_Rate_30ms_600K_4GPa], [600, Average_Shear_Stress_List_600K_30ms[4], Dissociation_Rate_30ms_600K_5GPa],
#         [700, Average_Shear_Stress_List_700K_30ms[0], Dissociation_Rate_30ms_700K_1GPa], [700, Average_Shear_Stress_List_700K_30ms[1], Dissociation_Rate_30ms_700K_2GPa], [700, Average_Shear_Stress_List_700K_30ms[2], Dissociation_Rate_30ms_700K_3GPa], [700, Average_Shear_Stress_List_700K_30ms[3], Dissociation_Rate_30ms_700K_4GPa], [700, Average_Shear_Stress_List_700K_30ms[4], Dissociation_Rate_30ms_700K_5GPa]]
#
# for item in data:
#     x_data.append(item[0])
#     y_data.append(item[1])
#     z_data.append(item[2])
#
# parameters, covariance = optimize.curve_fit(function, [x_data, y_data], z_data)
#
# # create surface function model
# # setup data points for calculating surface model
# model_x_data = np.linspace(min(x_data), max(x_data), 40)
# model_y_data = np.linspace(min(y_data), max(y_data), 40)
# # create coordinate arrays for vectorized evaluations
# X, Y = np.meshgrid(model_x_data, model_y_data)
# # calculate Z coordinate array
# Z = function(np.array([X, Y]), *parameters)
#
# z = []
# for row in Z:
#     row.sort()
#     z.append(row)
#
# z = np.array(z)
# # print(Z)
# # print('####################')
# # print(z)
#
# cm = plt.get_cmap("jet")
# fig = plt.figure()
# ax4 = plt.axes(projection='3d')
# ax4.plot_surface(X, Y, Z, cmap=cm, alpha=0.5, edgecolor= 'black', linewidth=0.3)
# ax4.set_title('3D Plot - Variation in Dissociation Rates - Alpha Fe, 30ms')
# ax4.scatter(x_data, y_data, z_data, color='black', alpha=1)
# ax4.set_xlabel('Temperature (K)')
# ax4.invert_xaxis()
# ax4.set_ylabel('Shear Stress (GPa)')
# ax4.set_zlabel('Dissociation Rate (per ns)')
# plt.show()
#
# # ########## CALCULATING THE ACTIVATION ENERGY, VOLUME AND PREFACTOR FROM 3D FIT ##########
# #
# logRates3D_30ms = zlogs
# shearstresses30ms = model_y_data
# Index = 0
# sigma = []
# params = []
# while Index < len(shearstresses30ms) - 1:
#     for logRates3Drow in logRates3D_30ms:
#         coef_sh_cf, coef_sh_pcov = optimize.curve_fit(linear, shearstresses30ms, logRates3Drow)
#         paramsrow = np.polyfit(shearstresses30ms, logRates3Drow, 1)
#         sigmarow = coef_sh_pcov[0, 0] ** 0.5
#         sigma.append(sigmarow)
#         params.append(paramsrow[0])
#         Index +=1
#
# #
# sigma = np.array(sigma)
# params = np.array(params)
# params = np.average(params)
# sigma = np.average(sigma)
#
# activation_vol_3D_30ms = (params) * (1.38065) * 500 * 1e-2
# print(f'3D Activation Volume, 30ms is {activation_vol_3D_30ms}')
#
# alpha = 0.05
# sigma = sigma ** 0.5
# dof = 2
# tval = t.ppf(1.0 - alpha / 2., dof)
# uncert = sigma * tval
# uncert_ActivationVolume_30ms = uncert * (1.38065) * 500 * 1e-2
# print(f'Activation Volume uncertainty for 3D fit at 30ms is {uncert_ActivationVolume_30ms}')
#
# logRates3D_30ms = zlogs
# temperatures = model_x_data
#
# inverse_temperatures = [1 / x for x in temperatures]
#
# Index = 0
# sigma = []
# SigmaA = []
# params = []
# interceptaverage = []
# while Index < len(inverse_temperatures) - 1:
#     for logRates3Drow in logRates3D_30ms:
#         coef_sh_cf, coef_sh_pcov = optimize.curve_fit(linear, inverse_temperatures, logRates3Drow)
#         paramsrow = np.polyfit(inverse_temperatures, logRates3Drow, 1)
#         sigmarow = coef_sh_pcov[0, 0] ** 0.5
#         sigma_A = coef_sh_pcov[1, 1] ** 0.5
#         sigma.append(sigmarow)
#         SigmaA.append(sigma_A)
#         params.append(paramsrow[0])
#         intercept = paramsrow[1]
#         interceptaverage.append((intercept))
#         Index +=1
#
# sigma = np.array(sigma)
# params = np.array(params)
# SigmaA = np.array(SigmaA)
# interceptaverage = np.array(interceptaverage)
# params_30ms = np.average(params)
# sigma = np.average(sigma)
# interceptaverage = np.average(interceptaverage)
# alpha = 0.05
# sigma = sigma ** 0.5
# tval = t.ppf(1.0 - alpha / 2., dof)
# uncert = sigma * tval
#
# Mu1GPa_30ms = np.average(Average_Mu_List_1GPa_30ms)
# Mu2GPa_30ms = np.average(Average_Mu_List_2GPa_30ms)
# Mu3GPa_30ms = np.average(Average_Mu_List_3GPa_30ms)
# Mu4GPa_30ms = np.average(Average_Mu_List_4GPa_30ms)
# Mu5GPa_30ms = np.average(Average_Mu_List_5GPa_30ms)
#
# MuAveragesDifferentPressures = np.array([Mu1GPa_30ms, Mu2GPa_30ms, Mu3GPa_30ms, Mu4GPa_30ms, Mu5GPa_30ms])
# AverageMu_30ms = np.average(MuAveragesDifferentPressures)
#
# ActivationEnergy_3D_30ms = (((3 * 1e9 * (AverageMu_30ms * activation_vol_3D_30ms) * 1e-30) - 1.381 * params_30ms * 1e-23) * 6.02214076 * (10**23)) / 1000
#
# print(f'Activation Energy for 3D fit at 30ms is {ActivationEnergy_3D_30ms}')
#
# error_3D_30ms = (np.sqrt((3 * 1e9 * np.average(AverageMu_30ms) * uncert_ActivationVolume_30ms * 1e-30) ** 2 + (1.381 * uncert * 1e-23) ** 2) * 6.02214076 * (10**23)) / 1000
# print(f"Activation_Energy Error at 30ms is {error_3D_30ms}")
#
# uncert_prefactor_3D_30ms = sigma_A * tval
# lnA_3D_30ms = np.log((np.exp(interceptaverage) * (10 ** 9)))
#
# print(f"ln(A) for 3D fit is {lnA_3D_30ms}")
# print(f"ln(A) uncertainty is {uncert_prefactor_3D_30ms}" + str(10 ** 9))
#
# ######### Plotting 3D fit vs 2D fit results along with their error margins ##########
# """
# Need to get a list with:
# - Values for Activation Volumes at different temperatures
# - Values for Activation Energies at different pressures
# - Values for ln(A) at different pressures
# - List of errors for each of the above quantities
# - Make a graph for each surface chemistry
# - Eventually will need to do the same for the different sliding speeds
#
# """
# Temperatures = [400, 500, 600, 700]
# Pressures = [1, 2, 3, 4, 5]
#
#
# Activation_Energies_30ms = [ActivationEnergy_1GPa_30ms, ActivationEnergy_2GPa_30ms, ActivationEnergy_3GPa_30ms, ActivationEnergy_4GPa_30ms, ActivationEnergy_5GPa_30ms]
# Activation_Energy_Errors_30ms = [err_E_1GPa_30ms, err_E_2GPa_30ms, err_E_3GPa_30ms, err_E_4GPa_30ms, err_E_5GPa_30ms]
# ActivationEnergy_3D_30ms = ActivationEnergy_3D_30ms
# ActivationEnergy_3D_error_30ms = error_3D_30ms
# ActivationEnergy_3D_error_UpperBound_Value_30ms = float(float(ActivationEnergy_3D_30ms) + float(ActivationEnergy_3D_error_30ms))
# ActivationEnergy_3D_error_LowerBound_Value_30ms = float(float(ActivationEnergy_3D_30ms) - float(ActivationEnergy_3D_error_30ms))
#
# Activation_Energy_Error_Plot, Ea2Dvs3d  = plt.subplots()
# Ea2Dvs3d.set_title('Comparison of Activation Energies from 2D and 3D Fits - AlphaFe, 30ms')
# Ea2Dvs3d.set_xlabel('Normal Stress (GPa)')
# Ea2Dvs3d.set_ylabel('Activation Energy')
# Ea2Dvs3d.scatter(Pressures, Activation_Energies_30ms)
# Ea2Dvs3d.errorbar(Pressures, Activation_Energies_30ms, yerr=Activation_Energy_Errors_30ms, linestyle="None", fmt='o', capsize=3)
# Ea2Dvs3d.axhline(y=ActivationEnergy_3D_30ms)
# Pressures = [0.5, 1, 2, 3, 4, 5, 5.5]
# Ea2Dvs3d.fill_between(Pressures, ActivationEnergy_3D_error_LowerBound_Value_30ms, ActivationEnergy_3D_error_UpperBound_Value_30ms, alpha=0.4)
# Ea2Dvs3d.set_xlim(0.5, 5.5)
# Ea2Dvs3d.set_ylim(0, 35)
# plt.show()
#
# ################ Activation Volume Errors #############################
#
# Activation_Volumes_30ms = [activation_vol_400K_30ms, activation_vol_500K_30ms, activation_vol_600K_30ms, activation_vol_700K_30ms]
# Activation_Volume_Errors_30ms = [uncert400_30ms, uncert500_30ms, uncert600_30ms, uncert700_30ms]
# Activation_Volume_3D_30ms = activation_vol_3D_30ms
# Activation_Volume_3D_Error_30ms = uncert_ActivationVolume_30ms
# ActivationVolume_3D_error_UpperBound_Value_30ms = float(float(Activation_Volume_3D_30ms) + float(Activation_Volume_3D_Error_30ms))
# ActivationVolume_3D_error_LowerBound_Value_30ms = float(float(Activation_Volume_3D_30ms) - float(Activation_Volume_3D_Error_30ms))
#
# Activation_Volume_Error_Plot, Av2Dvs3d  = plt.subplots()
# Av2Dvs3d.set_title('Comparison of Activation Volumes from 2D and 3D Fits - AlphaFe, 30ms')
# Av2Dvs3d.set_xlabel('Normal Stress(GPa)')
# Av2Dvs3d.set_ylabel('Activation Volume')
# Av2Dvs3d.scatter(Temperatures, Activation_Volumes_30ms)
# Av2Dvs3d.errorbar(Temperatures, Activation_Volumes_30ms, yerr=Activation_Volume_Errors_30ms, linestyle="None", fmt='o', capsize=3)
# Av2Dvs3d.axhline(y=Activation_Volume_3D_30ms)
# Temperatures = [350, 400, 500, 600, 700, 750]
# Av2Dvs3d.fill_between(Temperatures, ActivationVolume_3D_error_LowerBound_Value_30ms, ActivationVolume_3D_error_UpperBound_Value_30ms, alpha=0.4)
# Av2Dvs3d.set_xlim(350, 750)
# Av2Dvs3d.set_ylim(0, 30)
# plt.show()
#
# #################### Prefactor Errors #########################
# Pressures = [1, 2, 3, 4, 5]
# Prefactors_30ms = [lnA1GPa_30ms, lnA2GPa_30ms, lnA3GPa_30ms, lnA4GPa_30ms, lnA5GPa_30ms]
# PrefactorErrors_30ms = [uncert_A1GPa_30ms, uncert_A2GPa_30ms, uncert_A3GPa_30ms, uncert_A4GPa_30ms, uncert_A5GPa_30ms]
# lnA_3D_error_30ms = uncert_prefactor_3D_30ms
# lnA_3D_error_UpperBound_Value_30ms = float(float(lnA_3D_30ms) + float(lnA_3D_error_30ms))
# lnA_3D_error_LowerBound_Value_30ms = float(float(lnA_3D_30ms) - float(lnA_3D_error_30ms))
#
# Prefactor_Error_Plot, lnA2Dvs3d  = plt.subplots()
# lnA2Dvs3d.set_title('Comparison of Prefactors from 2D and 3D Fits - AlphaFe, 30ms')
# lnA2Dvs3d.set_xlabel('Normal Stress(GPa)')
# lnA2Dvs3d.set_ylabel('Prefactor')
# lnA2Dvs3d.scatter(Pressures, Prefactors_30ms)
# lnA2Dvs3d.errorbar(Pressures, Prefactors_30ms, yerr=PrefactorErrors_30ms, linestyle="None", fmt='o', capsize=3)
# lnA2Dvs3d.axhline(y=lnA_3D_30ms)
# Pressures = [0.5, 1, 2, 3, 4, 5, 5.5]
# lnA2Dvs3d.fill_between(Pressures, lnA_3D_error_LowerBound_Value_30ms, lnA_3D_error_UpperBound_Value_30ms, alpha=0.4)
# lnA2Dvs3d.set_xlim(0.5, 5.5)
# lnA2Dvs3d.set_ylim(20, 30)
# plt.show()
#
# ########## Checking for Kinetic Compensation Effect ##########
#
# KineticCompeEffectPlot, EavsLnA  = plt.subplots()
# EavsLnA.set_title('Activation Energy vs Prefactor - AlphaFe, 30ms')
# EavsLnA.set_xlabel('Ea')
# EavsLnA.set_ylabel('ln A')
# EavsLnA.scatter(Activation_Energies_30ms, Prefactors_20ms)
# EavsLnA.set_ylim(24, 27)
# plt.show()
#
# ############################################# 40ms #####################################
#
#
# Index = 0
# Temperatures = ["500K", "600K", "700K"]
# Pressures = ['2GPa', '3GPa', '4GPa', '5GPa']
# Speeds = ['20ms', '30ms', '40ms', '50ms']
#
# Dissociation_Rate_40ms_400K_1GPa, LogRate_40ms_400K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_400K_1GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_40ms_400K_2GPa, LogRate_40ms_400K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_400K_2GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_40ms_400K_3GPa, LogRate_40ms_400K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_400K_3GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_40ms_400K_4GPa, LogRate_40ms_400K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_400K_4GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_40ms_400K_5GPa, LogRate_40ms_400K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_400K_5GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rates_40ms_400K = [Dissociation_Rate_40ms_400K_1GPa, Dissociation_Rate_40ms_400K_2GPa, Dissociation_Rate_40ms_400K_3GPa, Dissociation_Rate_40ms_400K_4GPa, Dissociation_Rate_40ms_400K_5GPa]
# Log_Rates_400K_40ms = [LogRate_40ms_400K_1GPa, LogRate_40ms_400K_2GPa, LogRate_40ms_400K_3GPa, LogRate_40ms_400K_4GPa, LogRate_40ms_400K_5GPa]
#
# Dissociation_Rate_40ms_500K_1GPa, LogRate_40ms_500K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_500K_1GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_40ms_500K_2GPa, LogRate_40ms_500K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_500K_2GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_40ms_500K_3GPa, LogRate_40ms_500K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_500K_3GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_40ms_500K_4GPa, LogRate_40ms_500K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_500K_4GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_40ms_500K_5GPa, LogRate_40ms_500K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_500K_5GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rates_40ms_500K = [Dissociation_Rate_40ms_500K_1GPa, Dissociation_Rate_40ms_500K_2GPa, Dissociation_Rate_40ms_500K_3GPa, Dissociation_Rate_40ms_500K_4GPa, Dissociation_Rate_40ms_500K_5GPa]
# Log_Rates_500K_40ms = [LogRate_40ms_500K_1GPa, LogRate_40ms_500K_2GPa, LogRate_40ms_500K_3GPa, LogRate_40ms_500K_4GPa, LogRate_40ms_500K_5GPa]
#
# Dissociation_Rate_40ms_600K_1GPa, LogRate_40ms_600K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_600K_1GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_40ms_600K_2GPa, LogRate_40ms_600K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_600K_2GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_40ms_600K_3GPa, LogRate_40ms_600K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_600K_3GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_40ms_600K_4GPa, LogRate_40ms_600K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_600K_4GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_40ms_600K_5GPa, LogRate_40ms_600K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_600K_5GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rates_40ms_600K = [Dissociation_Rate_40ms_600K_1GPa, Dissociation_Rate_40ms_600K_2GPa, Dissociation_Rate_40ms_600K_3GPa, Dissociation_Rate_40ms_600K_4GPa, Dissociation_Rate_40ms_600K_5GPa]
# Log_Rates_600K_40ms = [LogRate_40ms_600K_1GPa, LogRate_40ms_600K_2GPa, LogRate_40ms_600K_3GPa, LogRate_40ms_600K_4GPa, LogRate_40ms_600K_5GPa]
#
# Dissociation_Rate_40ms_700K_1GPa, LogRate_40ms_700K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_700K_1GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_40ms_700K_2GPa, LogRate_40ms_700K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_700K_2GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_40ms_700K_3GPa, LogRate_40ms_700K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_700K_3GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_40ms_700K_4GPa, LogRate_40ms_700K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_700K_4GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rate_40ms_700K_5GPa, LogRate_40ms_700K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_700K_5GPa_Coefficients[Index], Cutoff=Cutoff)
# Dissociation_Rates_40ms_700K = [Dissociation_Rate_40ms_700K_1GPa, Dissociation_Rate_40ms_700K_2GPa, Dissociation_Rate_40ms_700K_3GPa, Dissociation_Rate_40ms_700K_4GPa, Dissociation_Rate_40ms_700K_5GPa]
# Log_Rates_700K_40ms = [LogRate_40ms_700K_1GPa, LogRate_40ms_700K_2GPa, LogRate_40ms_700K_3GPa, LogRate_40ms_700K_4GPa, LogRate_40ms_700K_5GPa]
#
# Dissociation_Rates_1GPa = [Dissociation_Rates_40ms_400K[0], Dissociation_Rates_40ms_500K[0], Dissociation_Rates_40ms_600K[0], Dissociation_Rates_40ms_700K[0]]
# Dissociation_Rates_2GPa = [Dissociation_Rates_40ms_400K[1], Dissociation_Rates_40ms_500K[1], Dissociation_Rates_40ms_600K[1], Dissociation_Rates_40ms_700K[1]]
# Dissociation_Rates_3GPa = [Dissociation_Rates_40ms_400K[2], Dissociation_Rates_40ms_500K[2], Dissociation_Rates_40ms_600K[2], Dissociation_Rates_40ms_700K[2]]
# Dissociation_Rates_4GPa = [Dissociation_Rates_40ms_400K[3], Dissociation_Rates_40ms_500K[3], Dissociation_Rates_40ms_600K[3], Dissociation_Rates_40ms_700K[3]]
# Dissociation_Rates_5GPa = [Dissociation_Rates_40ms_400K[4], Dissociation_Rates_40ms_500K[4], Dissociation_Rates_40ms_600K[4], Dissociation_Rates_40ms_700K[4]]
#
# Log_Rates_1GPa_40ms = [LogRate_40ms_400K_1GPa, LogRate_40ms_500K_1GPa, LogRate_40ms_600K_1GPa, LogRate_40ms_700K_1GPa]
# Log_Rates_2GPa_40ms = [LogRate_40ms_400K_2GPa, LogRate_40ms_500K_2GPa, LogRate_40ms_600K_2GPa, LogRate_40ms_700K_2GPa]
# Log_Rates_3GPa_40ms = [LogRate_40ms_400K_3GPa, LogRate_40ms_500K_3GPa, LogRate_40ms_600K_3GPa, LogRate_40ms_700K_3GPa]
# Log_Rates_4GPa_40ms = [LogRate_40ms_400K_4GPa, LogRate_40ms_500K_4GPa, LogRate_40ms_600K_4GPa, LogRate_40ms_700K_4GPa]
# Log_Rates_5GPa_40ms = [LogRate_40ms_400K_5GPa, LogRate_40ms_500K_5GPa, LogRate_40ms_600K_5GPa, LogRate_40ms_700K_5GPa]
#
# EquilibriumFactor = [79, 99] # How many rows (out of 99) to ignore before calculating shear stress/friction coefficient, as it won't stabilise until after a certain number of timesteps
#
# def get_average_shear_normal_stress_and_average_mu_constant_temperature(Temperature, Pressures, EquilibriumFactor, Speed):
#     Friction_Coefficient_Dataframe_Unnamed = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{Speed}/{Temperature}/1GPa/'
#                                 'fc_ave.dump'.format(Speed=Speed, Temperature=Temperature), sep=' ')
#     Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe_Unnamed.rename(columns={'v_s_bot' : 'Shear Stress 1GPa', 'v_p_bot' : 'Normal Stress 1GPa'})
#
#     for P in Pressures:
#         Dataframe = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{Speed}/{Temperature}/{P}/'
#                                 'fc_ave.dump'.format(Speed=Speed, Temperature=Temperature, P=P), sep=' ')
#         Big_DataframeP = Dataframe.rename(columns= {'Timestep': 'Timestep {}'.format(P),
#                                                         'v_s_bot': 'Shear Stress {}'.format(P),
#                                                         'v_p_bot': 'Normal Stress {}'.format(P)})
#
#         Friction_Coefficient_Dataframe = pd.concat([Friction_Coefficient_Dataframe, Big_DataframeP], axis =1)
#         Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe.dropna()
#
#
#     #print(Friction_Coefficient_Dataframe)
#     Mu_Final_Dataframe = Friction_Coefficient_Dataframe.iloc[:, [0, 1, 2, 4, 5, 7, 8, 10, 11, 13, 14]]
#     Mu_Final_Dataframe = Mu_Final_Dataframe.iloc[EquilibriumFactor[0]:EquilibriumFactor[1], :]
#     #print(Mu_Final_Dataframe)
#
#     ShearStressMeans = Mu_Final_Dataframe[['Shear Stress 1GPa', 'Shear Stress 2GPa', 'Shear Stress 3GPa', 'Shear Stress 4GPa', 'Shear Stress 5GPa']].mean()
#     Average_Shear_Stress_Dictionary = ShearStressMeans.to_dict()
#     #print(ShearStressMeans)
#     NormalStressMeans = Mu_Final_Dataframe[['Normal Stress 1GPa', 'Normal Stress 2GPa', 'Normal Stress 3GPa', 'Normal Stress 4GPa', 'Normal Stress 5GPa']].mean()
#     NormalStressMeans = NormalStressMeans.to_dict()
#     #print(NormalStressMeans)
#
#     Average_Mu_Dictionary = {}
#
#     NormalStressMeansList = list(NormalStressMeans.values())
#     Normal_Stress = NormalStressMeans.get('Normal Stress 1GPa')
#     Shear_Stress = ShearStressMeans.get('Shear Stress 1GPa')
#     Average_Mu = Shear_Stress / Normal_Stress
#     Average_Mu_Dictionary.update({'Average Mu 1GPa': Average_Mu})
#
#     for P in Pressures:
#
#         Normal_Stress = NormalStressMeans.get('Normal Stress {}'.format(P))
#         Shear_Stress = ShearStressMeans.get('Shear Stress {}'.format(P))
#         Average_Mu = Shear_Stress / Normal_Stress
#         Average_Mu_Dictionary.update({'Average Mu {}'.format(P): Average_Mu})
#
#
#     Average_Shear_Stress_List = list(Average_Shear_Stress_Dictionary.values())
#     #print(Average_Shear_Stress_List)
#     Average_Mu_List = list(Average_Mu_Dictionary.values())
#     Average_Shear_Stress_List = [x / 10000 for x in Average_Shear_Stress_List] # Conversion to GPa
#     NormalStressMeansList = [x/10000 for x in NormalStressMeansList]
#     #print(Average_Shear_Stress_List)
#
#     return Average_Shear_Stress_List, Average_Mu_List, NormalStressMeansList
#
# ######## Getting Average Shear Stress, Friction Coefficient and Normal Stress #################
# Average_Shear_Stress_List_400K_40ms, Average_Mu_List_400K_40ms, NormalStressMeans_400K_40ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("400K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="40ms")
# Average_Shear_Stress_List_500K_40ms, Average_Mu_List_500K_40ms, NormalStressMeans_500K_40ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("500K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="40ms")
# Average_Shear_Stress_List_600K_40ms, Average_Mu_List_600K_40ms, NormalStressMeans_600K_40ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("600K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="40ms")
# Average_Shear_Stress_List_700K_40ms, Average_Mu_List_700K_40ms, NormalStressMeans_700K_40ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("700K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="40ms")
#
# Average_Mu_List_1GPa_40ms = [Average_Mu_List_400K_40ms[0], Average_Mu_List_500K_40ms[0], Average_Mu_List_600K_40ms[0], Average_Mu_List_700K_40ms[0]]
# Average_Mu_List_2GPa_40ms = [Average_Mu_List_400K_40ms[1], Average_Mu_List_500K_40ms[1], Average_Mu_List_600K_40ms[1], Average_Mu_List_700K_40ms[1]]
# Average_Mu_List_3GPa_40ms = [Average_Mu_List_400K_40ms[2], Average_Mu_List_500K_40ms[2], Average_Mu_List_600K_40ms[2], Average_Mu_List_700K_40ms[2]]
# Average_Mu_List_4GPa_40ms = [Average_Mu_List_400K_40ms[3], Average_Mu_List_500K_40ms[3], Average_Mu_List_600K_40ms[3], Average_Mu_List_700K_40ms[3]]
# Average_Mu_List_5GPa_40ms = [Average_Mu_List_400K_40ms[4], Average_Mu_List_500K_40ms[4], Average_Mu_List_600K_40ms[4], Average_Mu_List_700K_40ms[4]]
#
# Average_Shear_Stress_List_1GPa_40ms = [Average_Shear_Stress_List_400K_40ms[0], Average_Shear_Stress_List_500K_40ms[0], Average_Shear_Stress_List_600K_40ms[0], Average_Shear_Stress_List_700K_40ms[0]]
# Average_Shear_Stress_List_2GPa_40ms = [Average_Shear_Stress_List_400K_40ms[1], Average_Shear_Stress_List_500K_40ms[1], Average_Shear_Stress_List_600K_40ms[1], Average_Shear_Stress_List_700K_40ms[1]]
# Average_Shear_Stress_List_3GPa_40ms = [Average_Shear_Stress_List_400K_40ms[2], Average_Shear_Stress_List_500K_40ms[2], Average_Shear_Stress_List_600K_40ms[2], Average_Shear_Stress_List_700K_40ms[2]]
# Average_Shear_Stress_List_4GPa_40ms = [Average_Shear_Stress_List_400K_40ms[3], Average_Shear_Stress_List_500K_40ms[3], Average_Shear_Stress_List_600K_40ms[3], Average_Shear_Stress_List_700K_40ms[3]]
# Average_Shear_Stress_List_5GPa_40ms = [Average_Shear_Stress_List_400K_40ms[4], Average_Shear_Stress_List_500K_40ms[4], Average_Shear_Stress_List_600K_40ms[4], Average_Shear_Stress_List_700K_40ms[4]]
#
# plot_shear_stress_vs_normal_stress(Average_Shear_Stress_List_400K_40ms, Average_Shear_Stress_List_500K_40ms, Average_Shear_Stress_List_600K_40ms, Average_Shear_Stress_List_700K_40ms,
#                                  "400K", "500K", "600K", "700K", Speed="40ms")
#
# ########################## Plotting ln(Rates) vs Shear Stress #####################################
# x = np.array([0, 1, 2, 3, 4, 5])
# params400K_40ms = np.polyfit(Average_Shear_Stress_List_400K_40ms, Log_Rates_400K_40ms, 1)
# params500K_40ms = np.polyfit(Average_Shear_Stress_List_500K_40ms, Log_Rates_500K_40ms, 1)
# params600K_40ms = np.polyfit(Average_Shear_Stress_List_600K_40ms, Log_Rates_600K_40ms, 1)
# params700K_40ms = np.polyfit(Average_Shear_Stress_List_700K_40ms, Log_Rates_700K_40ms, 1)
#
# RatesvsShear, RvS3 = plt.subplots()
# RvS3.set_title('Log of Dissociation Rates vs Shear Stress, 40ms')
# RvS3.set_xlabel('Shear Stress(GPa)')
# RvS3.set_ylabel('Log of Dissociation Rate (per nanosecond)')
# RvS3.scatter(Average_Shear_Stress_List_400K_40ms, Log_Rates_400K_40ms)
# RvS3.scatter(Average_Shear_Stress_List_500K_40ms, Log_Rates_500K_40ms)
# RvS3.scatter(Average_Shear_Stress_List_600K_40ms, Log_Rates_600K_40ms)
# RvS3.scatter(Average_Shear_Stress_List_700K_40ms, Log_Rates_700K_40ms)
# RvS3.plot(x, params400K_40ms[0] * x + params400K_40ms[1], label='400K Fitted')
# RvS3.plot(x, params500K_40ms[0] * x + params500K_40ms[1], label='500K Fitted')
# RvS3.plot(x, params600K_40ms[0] * x + params600K_40ms[1], label='600K Fitted')
# RvS3.plot(x, params700K_40ms[0] * x + params700K_40ms[1], label='700K Fitted')
# RvS3.set_xlim(0, 2)
# RvS3.set_ylim(-0.5, 4)
# RvS3.legend(loc='lower right')
# plt.show()
#
# ####### Calculate Activation Volume, Using  Carlos' conversion to get in Angstrom^3 #################
#
# activation_vol_400K_40ms = (params400K_40ms[0]) * (1.38065) * 400 * 1e-2
# activation_vol_500K_40ms = (params500K_40ms[0]) * (1.38065) * 500 * 1e-2
# activation_vol_600K_40ms = (params600K_40ms[0]) * (1.38065) * 600 * 1e-2
# activation_vol_700K_40ms = (params700K_40ms[0]) * (1.38065) * 700 * 1e-2
#
# alpha = 0.05
# def linear(x, m, n):
#     return m * x + n
#
# coef_sh_cf, coef_sh_pcov = optimize.curve_fit(linear, Average_Shear_Stress_List_700K_40ms, Log_Rates_700K_40ms)
# sigma = coef_sh_pcov[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_700K_40ms) - len(params700K_40ms))
#
# tval = t.ppf(1.0 - alpha / 2., dof)
# uncert = sigma * tval
# uncert400_40ms = uncert * (1.38065) * 400 * 1e-2
# uncert500_40ms = uncert * (1.38065) * 500 * 1e-2
# uncert600_40ms = uncert * (1.38065) * 600 * 1e-2
# uncert700_40ms = uncert * (1.38065) * 700 * 1e-2
#
# print(f'Activation Volume uncertainty at 400K, 40ms is {uncert400_40ms}')
# print(f'Activation Volume uncertainty at 500K, 40ms is {uncert500_40ms}')
# print(f'Activation Volume uncertainty at 600K, 40ms is {uncert600_40ms}')
# print(f'Activation Volume uncertainty at 700K, 40ms is {uncert700_40ms}')
#
# Activation_Volumes_40ms = [activation_vol_400K_40ms, activation_vol_500K_40ms, activation_vol_600K_40ms, activation_vol_700K_40ms]
# mean_actv_40ms = np.average(Activation_Volumes_40ms)
#
# print('Activation Volume 400K, 40ms = ' + str(activation_vol_400K_40ms))
# print('Activation Volume 500K, 40ms = ' + str(activation_vol_500K_40ms))
# print('Activation Volume 600K, 40ms = ' + str(activation_vol_600K_40ms))
# print('Activation Volume 700K, 40ms = ' + str(activation_vol_700K_40ms))
# #
# ############ Plotting lnk vs 1000/T  #########################
#
# Temperatures = [400, 500, 600, 700]
# Inverse_Temperatures = np.array([1/x for x in Temperatures])
#
# trend1GPa_40ms = np.polyfit(Inverse_Temperatures, Log_Rates_1GPa_40ms, 1)
# trend2GPa_40ms = np.polyfit(Inverse_Temperatures, Log_Rates_2GPa_40ms, 1)
# trend3GPa_40ms = np.polyfit(Inverse_Temperatures, Log_Rates_3GPa_40ms, 1)
# trend4GPa_40ms = np.polyfit(Inverse_Temperatures, Log_Rates_4GPa_40ms, 1)
# trend5GPa_40ms = np.polyfit(Inverse_Temperatures, Log_Rates_5GPa_40ms, 1)
#
# fig1, ax1 = plt.subplots()
# ax1.set_title('Log of Dissociation Rates against Inverse of Temperatures, 40ms')
# ax1.set_xlabel('1000/T (K-1)')
# ax1.set_ylabel('ln(Rate) (ns-1)')
# ax1.scatter(Inverse_Temperatures, Log_Rates_1GPa_40ms)
# ax1.scatter(Inverse_Temperatures, Log_Rates_2GPa_40ms)
# ax1.scatter(Inverse_Temperatures, Log_Rates_3GPa_40ms)
# ax1.scatter(Inverse_Temperatures, Log_Rates_4GPa_40ms)
# ax1.scatter(Inverse_Temperatures, Log_Rates_5GPa_40ms)
#
# Fit1GPa_40ms = np.poly1d(trend1GPa_40ms)
# Fit2GPa_40ms = np.poly1d(trend2GPa_40ms)
# Fit3GPa_40ms = np.poly1d(trend3GPa_40ms)
# Fit4GPa_40ms = np.poly1d(trend4GPa_40ms)
# Fit5GPa_40ms = np.poly1d(trend5GPa_40ms)
#
# ax1.plot(Inverse_Temperatures, Fit1GPa_40ms(Inverse_Temperatures), label='1GPa')
# ax1.plot(Inverse_Temperatures, Fit2GPa_40ms(Inverse_Temperatures), label='2GPa')
# ax1.plot(Inverse_Temperatures, Fit3GPa_40ms(Inverse_Temperatures), label='3GPa')
# ax1.plot(Inverse_Temperatures, Fit4GPa_40ms(Inverse_Temperatures), label='4GPa')
# ax1.plot(Inverse_Temperatures, Fit5GPa_40ms(Inverse_Temperatures), label='5GPa')
# ax1.legend()
# plt.show()
#
# ActivationEnergy_1GPa_40ms = (((1 * 1e9 * (np.average(Average_Mu_List_1GPa_40ms) * mean_actv_40ms) * 1e-30) - 1.381 * trend1GPa_40ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
# ActivationEnergy_2GPa_40ms = (((2 * 1e9 * (np.average(Average_Mu_List_2GPa_40ms) * mean_actv_40ms) * 1e-30) - 1.381 * trend2GPa_40ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
# ActivationEnergy_3GPa_40ms = (((3 * 1e9 * (np.average(Average_Mu_List_3GPa_40ms) * mean_actv_40ms) * 1e-30) - 1.381 * trend3GPa_40ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
# ActivationEnergy_4GPa_40ms = (((4 * 1e9 * (np.average(Average_Mu_List_4GPa_40ms) * mean_actv_40ms) * 1e-30) - 1.381 * trend4GPa_40ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
# ActivationEnergy_5GPa_40ms = (((5 * 1e9 * (np.average(Average_Mu_List_5GPa_40ms) * mean_actv_40ms) * 1e-30) - 1.381 * trend5GPa_40ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
#
# print(np.average(Average_Mu_List_1GPa_40ms))
# print(np.average(Average_Mu_List_2GPa_40ms))
# print(np.average(Average_Mu_List_3GPa_40ms))
# print(np.average(Average_Mu_List_4GPa_40ms))
# print(np.average(Average_Mu_List_5GPa_40ms))
#
# lnA1GPa_40ms = np.log((np.exp(trend1GPa_40ms[1]) * (10 ** 9)))
# lnA2GPa_40ms = np.log((np.exp(trend2GPa_40ms[1]) * (10 ** 9)))
# lnA3GPa_40ms = np.log((np.exp(trend3GPa_40ms[1]) * (10 ** 9)))
# lnA4GPa_40ms = np.log((np.exp(trend4GPa_40ms[1]) * (10 ** 9)))
# lnA5GPa_40ms = np.log((np.exp(trend5GPa_40ms[1]) * (10 ** 9)))
#
# print(f"ln(A) at 1GPa, 40ms is {lnA1GPa_40ms}")
# print(f"ln(A) at 2GPa, 40ms is {lnA2GPa_40ms}")
# print(f"ln(A) at 3GPa, 40ms is {lnA3GPa_40ms}")
# print(f"ln(A) at 4GPa, 40ms is {lnA4GPa_40ms}")
# print(f"ln(A) at 5GPa, 40ms is {lnA5GPa_40ms}")
#
# print('Activation Energy 1GPa, 40ms =' + str(ActivationEnergy_1GPa_40ms))
# print('Activation Energy 2GPa, 40ms =' + str(ActivationEnergy_2GPa_40ms))
# print('Activation Energy 3GPa, 40ms =' + str(ActivationEnergy_3GPa_40ms))
# print('Activation Energy 4GPa, 40ms =' + str(ActivationEnergy_4GPa_40ms))
# print('Activation Energy 5GPa, 40ms =' + str(ActivationEnergy_5GPa_40ms))
#
# def linear(x, m, n):
#     return m * x + n
#
# coef_p_cf_1GPa_40ms, coef_p_pcov_1GPa_40ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_1GPa_40ms)
# sigma_A1GPa_40ms = coef_p_pcov_1GPa_40ms[1, 1] ** 0.5
# sigma_m1GPa_40ms = coef_p_pcov_1GPa_40ms[0, 0] ** 0.5
# dof_40ms = max(0, len(Log_Rates_1GPa_40ms) - len(coef_p_cf_1GPa_40ms))
#
# alpha = 0.05
# tval_40ms = t.ppf(1.0 - alpha / 2., 4.1532)
#
# sigma_40ms = np.std(Activation_Volumes_40ms)
# error_actv_40ms = sigma_40ms * tval_40ms / np.sqrt(len(Activation_Volumes_40ms))
# uncert_A1GPa_40ms = sigma_A1GPa_40ms * tval
# uncert_m1GPa_40ms = sigma_m1GPa_40ms * tval
#
# coef_p_cf_2GPa_40ms, coef_p_pcov_2GPa_40ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_2GPa_40ms)
# sigma_A2GPa_40ms = coef_p_pcov_2GPa_40ms[1, 1] ** 0.5
# sigma_m2GPa_40ms = coef_p_pcov_2GPa_40ms[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_2GPa_40ms) - len(coef_p_cf_2GPa_40ms))
# uncert_A2GPa_40ms = sigma_A2GPa_40ms * tval_40ms
# uncert_m2GPa_40ms = sigma_m2GPa_40ms * tval_40ms
#
# coef_p_cf_3GPa_40ms, coef_p_pcov_3GPa_40ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_3GPa_40ms)
# sigma_A3GPa_40ms = coef_p_pcov_3GPa_40ms[1, 1] ** 0.5
# sigma_m3GPa_40ms = coef_p_pcov_3GPa_40ms[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_3GPa_40ms) - len(coef_p_cf_3GPa_40ms))
# uncert_A3GPa_40ms = sigma_A3GPa_40ms * tval_40ms
# uncert_m3GPa_40ms = sigma_m3GPa_40ms * tval_40ms
#
# coef_p_cf_4GPa_40ms, coef_p_pcov_4GPa_40ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_4GPa_40ms)
# sigma_A4GPa_40ms = coef_p_pcov_4GPa_40ms[1, 1] ** 0.5
# sigma_m4GPa_40ms = coef_p_pcov_4GPa_40ms[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_4GPa_40ms) - len(coef_p_cf_4GPa_40ms))
# uncert_A4GPa_40ms = sigma_A4GPa_40ms * tval_40ms
# uncert_m4GPa_40ms = sigma_m4GPa_40ms * tval_40ms
#
# coef_p_cf_5GPa_40ms, coef_p_pcov_5GPa_40ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_5GPa_40ms)
# sigma_A5GPa_40ms = coef_p_pcov_5GPa_40ms[1, 1] ** 0.5
# sigma_m5GPa_40ms = coef_p_pcov_5GPa_40ms[0, 0] ** 0.5
# dof = max(0, len(Log_Rates_5GPa_40ms) - len(coef_p_cf_5GPa_40ms))
# uncert_A5GPa_40ms = sigma_A5GPa_40ms * tval_40ms
# uncert_m5GPa_40ms = sigma_m5GPa_40ms * tval_40ms
#
# err_E_1GPa_40ms = (np.sqrt(((1 * 1e9 * np.average(Average_Mu_List_1GPa_40ms) * error_actv_40ms * 1e-30) ** 2 + (1.381 * uncert_m1GPa_40ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
# err_E_2GPa_40ms = (np.sqrt(((2 * 1e9 * np.average(Average_Mu_List_2GPa_40ms) * error_actv_40ms * 1e-30) ** 2 + (1.381 * uncert_m2GPa_40ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
# err_E_3GPa_40ms = (np.sqrt(((3 * 1e9 * np.average(Average_Mu_List_3GPa_40ms) * error_actv_40ms * 1e-30) ** 2 + (1.381 * uncert_m3GPa_40ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
# err_E_4GPa_40ms = (np.sqrt(((4 * 1e9 * np.average(Average_Mu_List_4GPa_40ms) * error_actv_40ms * 1e-30) ** 2 + (1.381 * uncert_m4GPa_40ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
# err_E_5GPa_40ms = (np.sqrt(((5 * 1e9 * np.average(Average_Mu_List_5GPa_40ms) * error_actv_40ms * 1e-30) ** 2 + (1.381 * uncert_m5GPa_40ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
#
# print(f"Error for Activation Energy at 1GPa = {err_E_1GPa_40ms}")
# print(f"Error for Activation Energy at 2GPa = {err_E_2GPa_40ms}")
# print(f"Error for Activation Energy at 3GPa = {err_E_3GPa_40ms}")
# print(f"Error for Activation Energy at 4GPa = {err_E_4GPa_40ms}")
# print(f"Error for Activation Energy at 5GPa = {err_E_5GPa_40ms}")
#
# uncertlnA1GPa_40ms = np.log(uncert_A1GPa_40ms * (10 ** 9))
# uncertlnA2GPa_40ms = np.log(uncert_A2GPa_40ms * (10 ** 9))
# uncertlnA3GPa_40ms = np.log(uncert_A3GPa_40ms * (10 ** 9))
# uncertlnA4GPa_40ms = np.log(uncert_A4GPa_40ms * (10 ** 9))
# uncertlnA5GPa_40ms = np.log(uncert_A5GPa_40ms * (10 ** 9))
#
# print(f"Error for ln A at 1GPa = {uncert_A1GPa_40ms}")
# print(f"Error for ln A at 2GPa = {uncert_A2GPa_40ms}")
# print(f"Error for ln A at 3GPa = {uncert_A3GPa_40ms}")
# print(f"Error for ln A at 4GPa = {uncert_A4GPa_40ms}")
# print(f"Error for ln A at 5GPa = {uncert_A5GPa_40ms}")
#
# pressures = ['2GPa', '3GPa', '4GPa', '5GPa']
#
# ########################## Plotting ln(Rates) vs Normal Stress #####################################
# x = np.array([0, 1, 2, 3, 4, 5])
# params400K_40ms = np.polyfit(NormalStressMeans_400K_40ms, Log_Rates_400K_40ms, 1)
# params500K_40ms = np.polyfit(NormalStressMeans_500K_40ms, Log_Rates_500K_40ms, 1)
# params600K_40ms = np.polyfit(NormalStressMeans_600K_40ms, Log_Rates_600K_40ms, 1)
# params700K_40ms = np.polyfit(NormalStressMeans_700K_40ms, Log_Rates_700K_40ms, 1)
#
# RatesvsNormal, RvN3 = plt.subplots()
# RvN3.set_title('Log of Dissociation Rates vs Normal Stress - Alpha Fe, 40ms')
# RvN3.set_xlabel('Normal Stress(GPa)')
# RvN3.set_ylabel('Log of Dissociation Rate (ns-1)')
# RvN3.scatter(NormalStressMeans_400K_40ms, Log_Rates_400K_40ms)
# RvN3.scatter(NormalStressMeans_500K_40ms, Log_Rates_500K_40ms)
# RvN3.scatter(NormalStressMeans_600K_40ms, Log_Rates_600K_40ms)
# RvN3.scatter(NormalStressMeans_700K_40ms, Log_Rates_700K_40ms)
# RvN3.plot(x, params400K_40ms[0] * x + params400K_40ms[1], label='400K Fitted')
# RvN3.plot(x, params500K_40ms[0] * x + params500K_40ms[1], label='500K Fitted')
# RvN3.plot(x, params600K_40ms[0] * x + params600K_40ms[1], label='600K Fitted')
# RvN3.plot(x, params700K_40ms[0] * x + params700K_40ms[1], label='700K Fitted')
# #RvN3.set_xlim(1, 5)
# # RvN3.set_ylim(0, 25)
# RvN3.legend(loc='lower right')
# plt.show()
#
# def function(data, a, b, c):
#     x = data[0]
#     y = data[1]
#     return a * (x**b) * (y**c) #TODO change fitting function
#
# x_data = []
# y_data = []
# z_data = []
#
# data = [[400, Average_Shear_Stress_List_400K_40ms[0], LogRate_40ms_400K_1GPa], [400, Average_Shear_Stress_List_400K_40ms[1], LogRate_40ms_400K_2GPa], [400, Average_Shear_Stress_List_400K_40ms[2], LogRate_40ms_400K_3GPa], [400, Average_Shear_Stress_List_400K_40ms[3], LogRate_40ms_400K_4GPa], [400, Average_Shear_Stress_List_400K_40ms[4], LogRate_40ms_400K_5GPa],
#         [500, Average_Shear_Stress_List_500K_40ms[0], LogRate_40ms_500K_1GPa], [500, Average_Shear_Stress_List_500K_40ms[1], LogRate_40ms_500K_2GPa], [500, Average_Shear_Stress_List_500K_40ms[2], LogRate_40ms_500K_3GPa], [500, Average_Shear_Stress_List_500K_40ms[3], LogRate_40ms_500K_4GPa], [500, Average_Shear_Stress_List_500K_40ms[4], LogRate_40ms_500K_5GPa],
#         [600, Average_Shear_Stress_List_600K_40ms[0], LogRate_40ms_600K_1GPa], [600, Average_Shear_Stress_List_600K_40ms[1], LogRate_40ms_600K_2GPa], [600, Average_Shear_Stress_List_600K_40ms[2], LogRate_40ms_600K_3GPa], [600, Average_Shear_Stress_List_600K_40ms[3], LogRate_40ms_600K_4GPa], [600, Average_Shear_Stress_List_600K_40ms[4], LogRate_40ms_600K_5GPa],
#         [700, Average_Shear_Stress_List_700K_40ms[0], LogRate_40ms_700K_1GPa], [700, Average_Shear_Stress_List_700K_40ms[1], LogRate_40ms_700K_2GPa], [700, Average_Shear_Stress_List_700K_40ms[2], LogRate_40ms_700K_3GPa], [700, Average_Shear_Stress_List_700K_40ms[3], LogRate_40ms_700K_4GPa], [700, Average_Shear_Stress_List_700K_40ms[4], LogRate_40ms_700K_5GPa]]
#
# for item in data:
#     x_data.append(item[0])
#     y_data.append(item[1])
#     z_data.append(item[2])
# #
# #
# parameters, covariance = optimize.curve_fit(function, [x_data, y_data], z_data)
#
# # create surface function model
# # setup data points for calculating surface model
# model_x_data = np.linspace(min(x_data), max(x_data), 40)
# model_y_data = np.linspace(min(y_data), max(y_data), 40)
#
# # create coordinate arrays for vectorized evaluations
# X, Y = np.meshgrid(model_x_data, model_y_data)
# # calculate Z coordinate array
# Z = function(np.array([X, Y]), *parameters)
#
# z = []
# for row in Z:
#     row.sort()
#     z.append(row)
#
# zlogs = np.array(z)
# import matplotlib.cm
# cm = plt.get_cmap("jet")
# fig = plt.figure()
# ax4 = plt.axes(projection='3d')
# ax4.plot_surface(X, Y, Z, cmap=cm, alpha=0.5, edgecolor='black', linewidth=0.3)
# ax4.set_title('3D Plot - Variation in Log of Dissociation Rates - Alpha Fe, 40ms')
# ax4.scatter(x_data, y_data, z_data, color='black', alpha=1)
# ax4.set_xlabel('Temperature (K)')
# ax4.invert_xaxis()
# ax4.set_ylabel('Shear Stress (GPa)')
# ax4.set_zlabel('Log of Dissociation Rate (per ns)')
# plt.show()
# plt.close(fig)
#
# x_data = []
# y_data = []
# z_data = []
#
# data = [[400, Average_Shear_Stress_List_400K_40ms[0], Dissociation_Rate_40ms_400K_1GPa], [400, Average_Shear_Stress_List_400K_40ms[1], Dissociation_Rate_40ms_400K_2GPa], [400, Average_Shear_Stress_List_400K_40ms[2], Dissociation_Rate_40ms_400K_3GPa], [400, Average_Shear_Stress_List_400K_40ms[3], Dissociation_Rate_40ms_400K_4GPa], [400, Average_Shear_Stress_List_400K_40ms[4], Dissociation_Rate_40ms_400K_5GPa],
#         [500, Average_Shear_Stress_List_500K_40ms[0], Dissociation_Rate_40ms_500K_1GPa], [500, Average_Shear_Stress_List_500K_40ms[1], Dissociation_Rate_40ms_500K_2GPa], [500, Average_Shear_Stress_List_500K_40ms[2], Dissociation_Rate_40ms_500K_3GPa], [500, Average_Shear_Stress_List_500K_40ms[3], Dissociation_Rate_40ms_500K_4GPa], [500, Average_Shear_Stress_List_500K_40ms[4], Dissociation_Rate_40ms_500K_5GPa],
#         [600, Average_Shear_Stress_List_600K_40ms[0], Dissociation_Rate_40ms_600K_1GPa], [600, Average_Shear_Stress_List_600K_40ms[1], Dissociation_Rate_40ms_600K_2GPa], [600, Average_Shear_Stress_List_600K_40ms[2], Dissociation_Rate_40ms_600K_3GPa], [600, Average_Shear_Stress_List_600K_40ms[3], Dissociation_Rate_40ms_600K_4GPa], [600, Average_Shear_Stress_List_600K_40ms[4], Dissociation_Rate_40ms_600K_5GPa],
#         [700, Average_Shear_Stress_List_700K_40ms[0], Dissociation_Rate_40ms_700K_1GPa], [700, Average_Shear_Stress_List_700K_40ms[1], Dissociation_Rate_40ms_700K_2GPa], [700, Average_Shear_Stress_List_700K_40ms[2], Dissociation_Rate_40ms_700K_3GPa], [700, Average_Shear_Stress_List_700K_40ms[3], Dissociation_Rate_40ms_700K_4GPa], [700, Average_Shear_Stress_List_700K_40ms[4], Dissociation_Rate_40ms_700K_5GPa]]
#
# for item in data:
#     x_data.append(item[0])
#     y_data.append(item[1])
#     z_data.append(item[2])
#
# parameters, covariance = optimize.curve_fit(function, [x_data, y_data], z_data)
#
# # create surface function model
# # setup data points for calculating surface model
# model_x_data = np.linspace(min(x_data), max(x_data), 40)
# model_y_data = np.linspace(min(y_data), max(y_data), 40)
# # create coordinate arrays for vectorized evaluations
# X, Y = np.meshgrid(model_x_data, model_y_data)
# # calculate Z coordinate array
# Z = function(np.array([X, Y]), *parameters)
#
# z = []
# for row in Z:
#     row.sort()
#     z.append(row)
#
# z = np.array(z)
# # print(Z)
# # print('####################')
# # print(z)
#
# cm = plt.get_cmap("jet")
# fig = plt.figure()
# ax4 = plt.axes(projection='3d')
# ax4.plot_surface(X, Y, Z, cmap=cm, alpha=0.5, edgecolor= 'black', linewidth=0.3)
# ax4.set_title('3D Plot - Variation in Dissociation Rates - Alpha Fe, 40ms')
# ax4.scatter(x_data, y_data, z_data, color='black', alpha=1)
# ax4.set_xlabel('Temperature (K)')
# ax4.invert_xaxis()
# ax4.set_ylabel('Shear Stress (GPa)')
# ax4.set_zlabel('Dissociation Rate (per ns)')
# plt.show()
#
# # ########## CALCULATING THE ACTIVATION ENERGY, VOLUME AND PREFACTOR FROM 3D FIT ##########
# #
# logRates3D_40ms = zlogs
# shearstresses40ms = model_y_data
# Index = 0
# sigma = []
# params = []
# while Index < len(shearstresses40ms) - 1:
#     for logRates3Drow in logRates3D_40ms:
#         coef_sh_cf, coef_sh_pcov = optimize.curve_fit(linear, shearstresses40ms, logRates3Drow)
#         paramsrow = np.polyfit(shearstresses40ms, logRates3Drow, 1)
#         sigmarow = coef_sh_pcov[0, 0] ** 0.5
#         sigma.append(sigmarow)
#         params.append(paramsrow[0])
#         Index +=1
#
# #
# sigma = np.array(sigma)
# params = np.array(params)
# params = np.average(params)
# sigma = np.average(sigma)
#
# activation_vol_3D_40ms = (params) * (1.38065) * 500 * 1e-2
# print(f'3D Activation Volume, 40ms is {activation_vol_3D_40ms}')
#
# alpha = 0.05
# sigma = sigma ** 0.5
# dof = 2
# tval = t.ppf(1.0 - alpha / 2., dof)
# uncert = sigma * tval
# uncert_ActivationVolume_40ms = uncert * (1.38065) * 500 * 1e-2
# print(f'Activation Volume uncertainty for 3D fit at 40ms is {uncert_ActivationVolume_40ms}')
#
# logRates3D_40ms = zlogs
# temperatures = model_x_data
#
# inverse_temperatures = [1 / x for x in temperatures]
#
# Index = 0
# sigma = []
# SigmaA = []
# params = []
# interceptaverage = []
# while Index < len(inverse_temperatures) - 1:
#     for logRates3Drow in logRates3D_40ms:
#         coef_sh_cf, coef_sh_pcov = optimize.curve_fit(linear, inverse_temperatures, logRates3Drow)
#         paramsrow = np.polyfit(inverse_temperatures, logRates3Drow, 1)
#         sigmarow = coef_sh_pcov[0, 0] ** 0.5
#         sigma_A = coef_sh_pcov[1, 1] ** 0.5
#         sigma.append(sigmarow)
#         SigmaA.append(sigma_A)
#         params.append(paramsrow[0])
#         intercept = paramsrow[1]
#         interceptaverage.append((intercept))
#         Index +=1
#
# sigma = np.array(sigma)
# params = np.array(params)
# SigmaA = np.array(SigmaA)
# interceptaverage = np.array(interceptaverage)
# params_40ms = np.average(params)
# sigma = np.average(sigma)
# interceptaverage = np.average(interceptaverage)
# alpha = 0.05
# sigma = sigma ** 0.5
# tval = t.ppf(1.0 - alpha / 2., dof)
# uncert = sigma * tval
#
# Mu1GPa_40ms = np.average(Average_Mu_List_1GPa_40ms)
# Mu2GPa_40ms = np.average(Average_Mu_List_2GPa_40ms)
# Mu3GPa_40ms = np.average(Average_Mu_List_3GPa_40ms)
# Mu4GPa_40ms = np.average(Average_Mu_List_4GPa_40ms)
# Mu5GPa_40ms = np.average(Average_Mu_List_5GPa_40ms)
#
# MuAveragesDifferentPressures = np.array([Mu1GPa_40ms, Mu2GPa_40ms, Mu3GPa_40ms, Mu4GPa_40ms, Mu5GPa_40ms])
# AverageMu_40ms = np.average(MuAveragesDifferentPressures)
#
# ActivationEnergy_3D_40ms = (((3 * 1e9 * (AverageMu_40ms * activation_vol_3D_40ms) * 1e-30) - 1.381 * params_40ms * 1e-23) * 6.02214076 * (10**23)) / 1000
#
# print(f'Activation Energy for 3D fit at 40ms is {ActivationEnergy_3D_40ms}')
#
# error_3D_40ms = (np.sqrt((3 * 1e9 * np.average(AverageMu_40ms) * uncert_ActivationVolume_40ms * 1e-30) ** 2 + (1.381 * uncert * 1e-23) ** 2) * 6.02214076 * (10**23)) / 1000
# print(f"Activation_Energy Error at 40ms is {error_3D_40ms}")
#
# uncert_prefactor_3D_40ms = sigma_A * tval
# lnA_3D_40ms = np.log((np.exp(interceptaverage) * (10 ** 9)))
#
# print(f"ln(A) for 3D fit is {lnA_3D_40ms}")
# print(f"ln(A) uncertainty is {uncert_prefactor_3D_40ms}" + str(10 ** 9))
#
# ######### Plotting 3D fit vs 2D fit results along with their error margins ##########
# """
# Need to get a list with:
# - Values for Activation Volumes at different temperatures
# - Values for Activation Energies at different pressures
# - Values for ln(A) at different pressures
# - List of errors for each of the above quantities
# - Make a graph for each surface chemistry
# - Eventually will need to do the same for the different sliding speeds
#
# """
# Temperatures = [400, 500, 600, 700]
# Pressures = [1, 2, 3, 4, 5]
#
#
# Activation_Energies_40ms = [ActivationEnergy_1GPa_40ms, ActivationEnergy_2GPa_40ms, ActivationEnergy_3GPa_40ms, ActivationEnergy_4GPa_40ms, ActivationEnergy_5GPa_40ms]
# Activation_Energy_Errors_40ms = [err_E_1GPa_40ms, err_E_2GPa_40ms, err_E_3GPa_40ms, err_E_4GPa_40ms, err_E_5GPa_40ms]
# ActivationEnergy_3D_40ms = ActivationEnergy_3D_40ms
# ActivationEnergy_3D_error_40ms = error_3D_40ms
# ActivationEnergy_3D_error_UpperBound_Value_40ms = float(float(ActivationEnergy_3D_40ms) + float(ActivationEnergy_3D_error_40ms))
# ActivationEnergy_3D_error_LowerBound_Value_40ms = float(float(ActivationEnergy_3D_40ms) - float(ActivationEnergy_3D_error_40ms))
#
# Activation_Energy_Error_Plot, Ea2Dvs3d  = plt.subplots()
# Ea2Dvs3d.set_title('Comparison of Activation Energies from 2D and 3D Fits - AlphaFe, 40ms')
# Ea2Dvs3d.set_xlabel('Normal Stress (GPa)')
# Ea2Dvs3d.set_ylabel('Activation Energy')
# Ea2Dvs3d.scatter(Pressures, Activation_Energies_40ms)
# Ea2Dvs3d.errorbar(Pressures, Activation_Energies_40ms, yerr=Activation_Energy_Errors_40ms, linestyle="None", fmt='o', capsize=3)
# Ea2Dvs3d.axhline(y=ActivationEnergy_3D_40ms)
# Pressures = [0.5, 1, 2, 3, 4, 5, 5.5]
# Ea2Dvs3d.fill_between(Pressures, ActivationEnergy_3D_error_LowerBound_Value_40ms, ActivationEnergy_3D_error_UpperBound_Value_40ms, alpha=0.4)
# Ea2Dvs3d.set_xlim(0.5, 5.5)
# Ea2Dvs3d.set_ylim(0, 35)
# plt.show()
#
# ################ Activation Volume Errors #############################
#
# Activation_Volumes_40ms = [activation_vol_400K_40ms, activation_vol_500K_40ms, activation_vol_600K_40ms, activation_vol_700K_40ms]
# Activation_Volume_Errors_40ms = [uncert400_40ms, uncert500_40ms, uncert600_40ms, uncert700_40ms]
# Activation_Volume_3D_40ms = activation_vol_3D_40ms
# Activation_Volume_3D_Error_40ms = uncert_ActivationVolume_40ms
# ActivationVolume_3D_error_UpperBound_Value_40ms = float(float(Activation_Volume_3D_40ms) + float(Activation_Volume_3D_Error_40ms))
# ActivationVolume_3D_error_LowerBound_Value_40ms = float(float(Activation_Volume_3D_40ms) - float(Activation_Volume_3D_Error_40ms))
#
# Activation_Volume_Error_Plot, Av2Dvs3d  = plt.subplots()
# Av2Dvs3d.set_title('Comparison of Activation Volumes from 2D and 3D Fits - AlphaFe, 40ms')
# Av2Dvs3d.set_xlabel('Normal Stress(GPa)')
# Av2Dvs3d.set_ylabel('Activation Volume')
# Av2Dvs3d.scatter(Temperatures, Activation_Volumes_40ms)
# Av2Dvs3d.errorbar(Temperatures, Activation_Volumes_40ms, yerr=Activation_Volume_Errors_40ms, linestyle="None", fmt='o', capsize=3)
# Av2Dvs3d.axhline(y=Activation_Volume_3D_40ms)
# Temperatures = [350, 400, 500, 600, 700, 750]
# Av2Dvs3d.fill_between(Temperatures, ActivationVolume_3D_error_LowerBound_Value_40ms, ActivationVolume_3D_error_UpperBound_Value_40ms, alpha=0.4)
# Av2Dvs3d.set_xlim(350, 750)
# Av2Dvs3d.set_ylim(0, 30)
# plt.show()
#
# #################### Prefactor Errors #########################
# Pressures = [1, 2, 3, 4, 5]
# Prefactors_40ms = [lnA1GPa_40ms, lnA2GPa_40ms, lnA3GPa_40ms, lnA4GPa_40ms, lnA5GPa_40ms]
# PrefactorErrors_40ms = [uncert_A1GPa_40ms, uncert_A2GPa_40ms, uncert_A3GPa_40ms, uncert_A4GPa_40ms, uncert_A5GPa_40ms]
# lnA_3D_error_40ms = uncert_prefactor_3D_40ms
# lnA_3D_error_UpperBound_Value_40ms = float(float(lnA_3D_40ms) + float(lnA_3D_error_40ms))
# lnA_3D_error_LowerBound_Value_40ms = float(float(lnA_3D_40ms) - float(lnA_3D_error_40ms))
#
# Prefactor_Error_Plot, lnA2Dvs3d  = plt.subplots()
# lnA2Dvs3d.set_title('Comparison of Prefactors from 2D and 3D Fits - AlphaFe, 40ms')
# lnA2Dvs3d.set_xlabel('Normal Stress(GPa)')
# lnA2Dvs3d.set_ylabel('Prefactor')
# lnA2Dvs3d.scatter(Pressures, Prefactors_40ms)
# lnA2Dvs3d.errorbar(Pressures, Prefactors_40ms, yerr=PrefactorErrors_40ms, linestyle="None", fmt='o', capsize=3)
# lnA2Dvs3d.axhline(y=lnA_3D_40ms)
# Pressures = [0.5, 1, 2, 3, 4, 5, 5.5]
# lnA2Dvs3d.fill_between(Pressures, lnA_3D_error_LowerBound_Value_40ms, lnA_3D_error_UpperBound_Value_40ms, alpha=0.4)
# lnA2Dvs3d.set_xlim(0.5, 5.5)
# lnA2Dvs3d.set_ylim(20, 30)
# plt.show()
#
# ########## Checking for Kinetic Compensation Effect ##########
#
# KineticCompeEffectPlot, EavsLnA  = plt.subplots()
# EavsLnA.set_title('Activation Energy vs Prefactor - AlphaFe, 40ms')
# EavsLnA.set_xlabel('Ea')
# EavsLnA.set_ylabel('ln A')
# EavsLnA.scatter(Activation_Energies_40ms, Prefactors_40ms)
# EavsLnA.set_ylim(24, 27)
# plt.show()

##################### 50ms ###################################

Index = 0
Temperatures = ["500K", "600K", "700K"]
Pressures = ['2GPa', '3GPa', '4GPa', '5GPa']
Speeds = ['20ms', '30ms', '40ms', '50ms']

Dissociation_Rate_40ms_400K_1GPa, LogRate_40ms_400K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_400K_1GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_400K_2GPa, LogRate_40ms_400K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_400K_2GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_400K_3GPa, LogRate_40ms_400K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_400K_3GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_400K_4GPa, LogRate_40ms_400K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_400K_4GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_400K_5GPa, LogRate_40ms_400K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_400K_5GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rates_40ms_400K = [Dissociation_Rate_40ms_400K_1GPa, Dissociation_Rate_40ms_400K_2GPa, Dissociation_Rate_40ms_400K_3GPa, Dissociation_Rate_40ms_400K_4GPa, Dissociation_Rate_40ms_400K_5GPa]
Log_Rates_400K_40ms = [LogRate_40ms_400K_1GPa, LogRate_40ms_400K_2GPa, LogRate_40ms_400K_3GPa, LogRate_40ms_400K_4GPa, LogRate_40ms_400K_5GPa]

Dissociation_Rate_40ms_500K_1GPa, LogRate_40ms_500K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_500K_1GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_500K_2GPa, LogRate_40ms_500K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_500K_2GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_500K_3GPa, LogRate_40ms_500K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_500K_3GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_500K_4GPa, LogRate_40ms_500K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_500K_4GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_500K_5GPa, LogRate_40ms_500K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_500K_5GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rates_40ms_500K = [Dissociation_Rate_40ms_500K_1GPa, Dissociation_Rate_40ms_500K_2GPa, Dissociation_Rate_40ms_500K_3GPa, Dissociation_Rate_40ms_500K_4GPa, Dissociation_Rate_40ms_500K_5GPa]
Log_Rates_500K_40ms = [LogRate_40ms_500K_1GPa, LogRate_40ms_500K_2GPa, LogRate_40ms_500K_3GPa, LogRate_40ms_500K_4GPa, LogRate_40ms_500K_5GPa]

Dissociation_Rate_40ms_600K_1GPa, LogRate_40ms_600K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_600K_1GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_600K_2GPa, LogRate_40ms_600K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_600K_2GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_600K_3GPa, LogRate_40ms_600K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_600K_3GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_600K_4GPa, LogRate_40ms_600K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_600K_4GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_600K_5GPa, LogRate_40ms_600K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_600K_5GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rates_40ms_600K = [Dissociation_Rate_40ms_600K_1GPa, Dissociation_Rate_40ms_600K_2GPa, Dissociation_Rate_40ms_600K_3GPa, Dissociation_Rate_40ms_600K_4GPa, Dissociation_Rate_40ms_600K_5GPa]
Log_Rates_600K_40ms = [LogRate_40ms_600K_1GPa, LogRate_40ms_600K_2GPa, LogRate_40ms_600K_3GPa, LogRate_40ms_600K_4GPa, LogRate_40ms_600K_5GPa]

Dissociation_Rate_40ms_700K_1GPa, LogRate_40ms_700K_1GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_700K_1GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_700K_2GPa, LogRate_40ms_700K_2GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_700K_2GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_700K_3GPa, LogRate_40ms_700K_3GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_700K_3GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_700K_4GPa, LogRate_40ms_700K_4GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_700K_4GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rate_40ms_700K_5GPa, LogRate_40ms_700K_5GPa = get_MATLABFIT_dissociation_rates(Timestep, Dissociation_Rates.Dissociation_Rate_40ms_700K_5GPa_Coefficients[Index], Cutoff=Cutoff)
Dissociation_Rates_40ms_700K = [Dissociation_Rate_40ms_700K_1GPa, Dissociation_Rate_40ms_700K_2GPa, Dissociation_Rate_40ms_700K_3GPa, Dissociation_Rate_40ms_700K_4GPa, Dissociation_Rate_40ms_700K_5GPa]
Log_Rates_700K_40ms = [LogRate_40ms_700K_1GPa, LogRate_40ms_700K_2GPa, LogRate_40ms_700K_3GPa, LogRate_40ms_700K_4GPa, LogRate_40ms_700K_5GPa]

Dissociation_Rates_1GPa = [Dissociation_Rates_40ms_400K[0], Dissociation_Rates_40ms_500K[0], Dissociation_Rates_40ms_600K[0], Dissociation_Rates_40ms_700K[0]]
Dissociation_Rates_2GPa = [Dissociation_Rates_40ms_400K[1], Dissociation_Rates_40ms_500K[1], Dissociation_Rates_40ms_600K[1], Dissociation_Rates_40ms_700K[1]]
Dissociation_Rates_3GPa = [Dissociation_Rates_40ms_400K[2], Dissociation_Rates_40ms_500K[2], Dissociation_Rates_40ms_600K[2], Dissociation_Rates_40ms_700K[2]]
Dissociation_Rates_4GPa = [Dissociation_Rates_40ms_400K[3], Dissociation_Rates_40ms_500K[3], Dissociation_Rates_40ms_600K[3], Dissociation_Rates_40ms_700K[3]]
Dissociation_Rates_5GPa = [Dissociation_Rates_40ms_400K[4], Dissociation_Rates_40ms_500K[4], Dissociation_Rates_40ms_600K[4], Dissociation_Rates_40ms_700K[4]]

Log_Rates_1GPa_40ms = [LogRate_40ms_400K_1GPa, LogRate_40ms_500K_1GPa, LogRate_40ms_600K_1GPa, LogRate_40ms_700K_1GPa]
Log_Rates_2GPa_40ms = [LogRate_40ms_400K_2GPa, LogRate_40ms_500K_2GPa, LogRate_40ms_600K_2GPa, LogRate_40ms_700K_2GPa]
Log_Rates_3GPa_40ms = [LogRate_40ms_400K_3GPa, LogRate_40ms_500K_3GPa, LogRate_40ms_600K_3GPa, LogRate_40ms_700K_3GPa]
Log_Rates_4GPa_40ms = [LogRate_40ms_400K_4GPa, LogRate_40ms_500K_4GPa, LogRate_40ms_600K_4GPa, LogRate_40ms_700K_4GPa]
Log_Rates_5GPa_40ms = [LogRate_40ms_400K_5GPa, LogRate_40ms_500K_5GPa, LogRate_40ms_600K_5GPa, LogRate_40ms_700K_5GPa]

EquilibriumFactor = [30, 50] # How many rows (out of 99) to ignore before calculating shear stress/friction coefficient, as it won't stabilise until after a certain number of timesteps

def get_average_shear_normal_stress_and_average_mu_constant_temperature(Temperature, Pressures, EquilibriumFactor, Speed):
    Friction_Coefficient_Dataframe_Unnamed = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{Speed}/{Temperature}/1GPa/'
                                'fc_ave.dump'.format(Speed=Speed, Temperature=Temperature), sep=' ')
    Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe_Unnamed.rename(columns={'v_s_bot' : 'Shear Stress 1GPa', 'v_p_bot' : 'Normal Stress 1GPa'})

    for P in Pressures:
        Dataframe = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{Speed}/{Temperature}/{P}/'
                                'fc_ave.dump'.format(Speed=Speed, Temperature=Temperature, P=P), sep=' ')
        Big_DataframeP = Dataframe.rename(columns= {'Timestep': 'Timestep {}'.format(P),
                                                        'v_s_bot': 'Shear Stress {}'.format(P),
                                                        'v_p_bot': 'Normal Stress {}'.format(P)})

        Friction_Coefficient_Dataframe = pd.concat([Friction_Coefficient_Dataframe, Big_DataframeP], axis =1)
        Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe.dropna()


    #print(Friction_Coefficient_Dataframe)
    Mu_Final_Dataframe = Friction_Coefficient_Dataframe.iloc[:, [0, 1, 2, 4, 5, 7, 8, 10, 11, 13, 14]]
    Mu_Final_Dataframe = Mu_Final_Dataframe.iloc[EquilibriumFactor[0]:EquilibriumFactor[1], :]
    #print(Mu_Final_Dataframe)

    ShearStressMeans = Mu_Final_Dataframe[['Shear Stress 1GPa', 'Shear Stress 2GPa', 'Shear Stress 3GPa', 'Shear Stress 4GPa', 'Shear Stress 5GPa']].mean()
    Average_Shear_Stress_Dictionary = ShearStressMeans.to_dict()
    #print(ShearStressMeans)
    NormalStressMeans = Mu_Final_Dataframe[['Normal Stress 1GPa', 'Normal Stress 2GPa', 'Normal Stress 3GPa', 'Normal Stress 4GPa', 'Normal Stress 5GPa']].mean()
    NormalStressMeans = NormalStressMeans.to_dict()
    #print(NormalStressMeans)

    Average_Mu_Dictionary = {}

    NormalStressMeansList = list(NormalStressMeans.values())
    Normal_Stress = NormalStressMeans.get('Normal Stress 1GPa')
    Shear_Stress = ShearStressMeans.get('Shear Stress 1GPa')
    Average_Mu = Shear_Stress / Normal_Stress
    Average_Mu_Dictionary.update({'Average Mu 1GPa': Average_Mu})

    for P in Pressures:

        Normal_Stress = NormalStressMeans.get('Normal Stress {}'.format(P))
        Shear_Stress = ShearStressMeans.get('Shear Stress {}'.format(P))
        Average_Mu = Shear_Stress / Normal_Stress
        Average_Mu_Dictionary.update({'Average Mu {}'.format(P): Average_Mu})


    Average_Shear_Stress_List = list(Average_Shear_Stress_Dictionary.values())
    #print(Average_Shear_Stress_List)
    Average_Mu_List = list(Average_Mu_Dictionary.values())
    Average_Shear_Stress_List = [x / 10000 for x in Average_Shear_Stress_List] # Conversion to GPa
    NormalStressMeansList = [x/10000 for x in NormalStressMeansList]
    #print(Average_Shear_Stress_List)

    return Average_Shear_Stress_List, Average_Mu_List, NormalStressMeansList

######## Getting Average Shear Stress, Friction Coefficient and Normal Stress #################
Average_Shear_Stress_List_400K_40ms, Average_Mu_List_400K_40ms, NormalStressMeans_400K_40ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("400K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="40ms")
Average_Shear_Stress_List_500K_40ms, Average_Mu_List_500K_40ms, NormalStressMeans_500K_40ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("500K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="40ms")
Average_Shear_Stress_List_600K_40ms, Average_Mu_List_600K_40ms, NormalStressMeans_600K_40ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("600K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="40ms")
Average_Shear_Stress_List_700K_40ms, Average_Mu_List_700K_40ms, NormalStressMeans_700K_40ms = get_average_shear_normal_stress_and_average_mu_constant_temperature("700K", Pressures=Pressures, EquilibriumFactor=EquilibriumFactor, Speed="40ms")

Average_Mu_List_1GPa_40ms = [Average_Mu_List_400K_40ms[0], Average_Mu_List_500K_40ms[0], Average_Mu_List_600K_40ms[0], Average_Mu_List_700K_40ms[0]]
Average_Mu_List_2GPa_40ms = [Average_Mu_List_400K_40ms[1], Average_Mu_List_500K_40ms[1], Average_Mu_List_600K_40ms[1], Average_Mu_List_700K_40ms[1]]
Average_Mu_List_3GPa_40ms = [Average_Mu_List_400K_40ms[2], Average_Mu_List_500K_40ms[2], Average_Mu_List_600K_40ms[2], Average_Mu_List_700K_40ms[2]]
Average_Mu_List_4GPa_40ms = [Average_Mu_List_400K_40ms[3], Average_Mu_List_500K_40ms[3], Average_Mu_List_600K_40ms[3], Average_Mu_List_700K_40ms[3]]
Average_Mu_List_5GPa_40ms = [Average_Mu_List_400K_40ms[4], Average_Mu_List_500K_40ms[4], Average_Mu_List_600K_40ms[4], Average_Mu_List_700K_40ms[4]]

Average_Shear_Stress_List_1GPa_40ms = [Average_Shear_Stress_List_400K_40ms[0], Average_Shear_Stress_List_500K_40ms[0], Average_Shear_Stress_List_600K_40ms[0], Average_Shear_Stress_List_700K_40ms[0]]
Average_Shear_Stress_List_2GPa_40ms = [Average_Shear_Stress_List_400K_40ms[1], Average_Shear_Stress_List_500K_40ms[1], Average_Shear_Stress_List_600K_40ms[1], Average_Shear_Stress_List_700K_40ms[1]]
Average_Shear_Stress_List_3GPa_40ms = [Average_Shear_Stress_List_400K_40ms[2], Average_Shear_Stress_List_500K_40ms[2], Average_Shear_Stress_List_600K_40ms[2], Average_Shear_Stress_List_700K_40ms[2]]
Average_Shear_Stress_List_4GPa_40ms = [Average_Shear_Stress_List_400K_40ms[3], Average_Shear_Stress_List_500K_40ms[3], Average_Shear_Stress_List_600K_40ms[3], Average_Shear_Stress_List_700K_40ms[3]]
Average_Shear_Stress_List_5GPa_40ms = [Average_Shear_Stress_List_400K_40ms[4], Average_Shear_Stress_List_500K_40ms[4], Average_Shear_Stress_List_600K_40ms[4], Average_Shear_Stress_List_700K_40ms[4]]

plot_shear_stress_vs_normal_stress(Average_Shear_Stress_List_400K_40ms, Average_Shear_Stress_List_500K_40ms, Average_Shear_Stress_List_600K_40ms, Average_Shear_Stress_List_700K_40ms,
                                 "400K", "500K", "600K", "700K", Speed="40ms")

########################## Plotting ln(Rates) vs Shear Stress #####################################
x = np.array([0, 1, 2, 3, 4, 5])
params400K_40ms = np.polyfit(Average_Shear_Stress_List_400K_40ms, Log_Rates_400K_40ms, 1)
params500K_40ms = np.polyfit(Average_Shear_Stress_List_500K_40ms, Log_Rates_500K_40ms, 1)
params600K_40ms = np.polyfit(Average_Shear_Stress_List_600K_40ms, Log_Rates_600K_40ms, 1)
params700K_40ms = np.polyfit(Average_Shear_Stress_List_700K_40ms, Log_Rates_700K_40ms, 1)

RatesvsShear, RvS3 = plt.subplots()
RvS3.set_title('Log of Dissociation Rates vs Shear Stress, 40ms')
RvS3.set_xlabel('Shear Stress(GPa)')
RvS3.set_ylabel('Log of Dissociation Rate (per nanosecond)')
RvS3.scatter(Average_Shear_Stress_List_400K_40ms, Log_Rates_400K_40ms)
RvS3.scatter(Average_Shear_Stress_List_500K_40ms, Log_Rates_500K_40ms)
RvS3.scatter(Average_Shear_Stress_List_600K_40ms, Log_Rates_600K_40ms)
RvS3.scatter(Average_Shear_Stress_List_700K_40ms, Log_Rates_700K_40ms)
RvS3.plot(x, params400K_40ms[0] * x + params400K_40ms[1], label='400K Fitted')
RvS3.plot(x, params500K_40ms[0] * x + params500K_40ms[1], label='500K Fitted')
RvS3.plot(x, params600K_40ms[0] * x + params600K_40ms[1], label='600K Fitted')
RvS3.plot(x, params700K_40ms[0] * x + params700K_40ms[1], label='700K Fitted')
RvS3.set_xlim(0, 2)
RvS3.set_ylim(-0.5, 4)
RvS3.legend(loc='lower right')
plt.show()

####### Calculate Activation Volume, Using  Carlos' conversion to get in Angstrom^3 #################

activation_vol_400K_40ms = (params400K_40ms[0]) * (1.38065) * 400 * 1e-2
activation_vol_500K_40ms = (params500K_40ms[0]) * (1.38065) * 500 * 1e-2
activation_vol_600K_40ms = (params600K_40ms[0]) * (1.38065) * 600 * 1e-2
activation_vol_700K_40ms = (params700K_40ms[0]) * (1.38065) * 700 * 1e-2

print(params400K_40ms)
print(params500K_40ms)
print(params600K_40ms)
print(params700K_40ms)

alpha = 0.05
def linear(x, m, n):
    return m * x + n

coef_sh_cf, coef_sh_pcov = optimize.curve_fit(linear, Average_Shear_Stress_List_700K_40ms, Log_Rates_700K_40ms)
sigma = coef_sh_pcov[0, 0] ** 0.5
dof = max(0, len(Log_Rates_700K_40ms) - len(params700K_40ms))

tval = t.ppf(1.0 - alpha / 2., dof)
uncert = sigma * tval
uncert400_40ms = uncert * (1.38065) * 400 * 1e-2
uncert500_40ms = uncert * (1.38065) * 500 * 1e-2
uncert600_40ms = uncert * (1.38065) * 600 * 1e-2
uncert700_40ms = uncert * (1.38065) * 700 * 1e-2

print(f'Activation Volume uncertainty at 400K, 40ms is {uncert400_40ms}')
print(f'Activation Volume uncertainty at 500K, 40ms is {uncert500_40ms}')
print(f'Activation Volume uncertainty at 600K, 40ms is {uncert600_40ms}')
print(f'Activation Volume uncertainty at 700K, 40ms is {uncert700_40ms}')

Activation_Volumes_40ms = [activation_vol_400K_40ms, activation_vol_500K_40ms, activation_vol_600K_40ms, activation_vol_700K_40ms]
mean_actv_40ms = np.average(Activation_Volumes_40ms)



print('Activation Volume 400K, 40ms = ' + str(activation_vol_400K_40ms))
print('Activation Volume 500K, 40ms = ' + str(activation_vol_500K_40ms))
print('Activation Volume 600K, 40ms = ' + str(activation_vol_600K_40ms))
print('Activation Volume 700K, 40ms = ' + str(activation_vol_700K_40ms))
#
############ Plotting lnk vs 1000/T  #########################

Temperatures = [400, 500, 600, 700]
Inverse_Temperatures = np.array([1/x for x in Temperatures])

trend1GPa_40ms = np.polyfit(Inverse_Temperatures, Log_Rates_1GPa_40ms, 1)
trend2GPa_40ms = np.polyfit(Inverse_Temperatures, Log_Rates_2GPa_40ms, 1)
trend3GPa_40ms = np.polyfit(Inverse_Temperatures, Log_Rates_3GPa_40ms, 1)
trend4GPa_40ms = np.polyfit(Inverse_Temperatures, Log_Rates_4GPa_40ms, 1)
trend5GPa_40ms = np.polyfit(Inverse_Temperatures, Log_Rates_5GPa_40ms, 1)

fig1, ax1 = plt.subplots()
ax1.set_title('Log of Dissociation Rates against Inverse of Temperatures, 40ms')
ax1.set_xlabel('1000/T (K-1)')
ax1.set_ylabel('ln(Rate) (ns-1)')
ax1.scatter(Inverse_Temperatures, Log_Rates_1GPa_40ms)
ax1.scatter(Inverse_Temperatures, Log_Rates_2GPa_40ms)
ax1.scatter(Inverse_Temperatures, Log_Rates_3GPa_40ms)
ax1.scatter(Inverse_Temperatures, Log_Rates_4GPa_40ms)
ax1.scatter(Inverse_Temperatures, Log_Rates_5GPa_40ms)

Fit1GPa_40ms = np.poly1d(trend1GPa_40ms)
Fit2GPa_40ms = np.poly1d(trend2GPa_40ms)
Fit3GPa_40ms = np.poly1d(trend3GPa_40ms)
Fit4GPa_40ms = np.poly1d(trend4GPa_40ms)
Fit5GPa_40ms = np.poly1d(trend5GPa_40ms)

ax1.plot(Inverse_Temperatures, Fit1GPa_40ms(Inverse_Temperatures), label='1GPa')
ax1.plot(Inverse_Temperatures, Fit2GPa_40ms(Inverse_Temperatures), label='2GPa')
ax1.plot(Inverse_Temperatures, Fit3GPa_40ms(Inverse_Temperatures), label='3GPa')
ax1.plot(Inverse_Temperatures, Fit4GPa_40ms(Inverse_Temperatures), label='4GPa')
ax1.plot(Inverse_Temperatures, Fit5GPa_40ms(Inverse_Temperatures), label='5GPa')
ax1.legend()
plt.show()

Mu1GPa_40ms = np.average(Average_Mu_List_1GPa_40ms)
Mu2GPa_40ms = np.average(Average_Mu_List_2GPa_40ms)
Mu3GPa_40ms = np.average(Average_Mu_List_3GPa_40ms)
Mu4GPa_40ms = np.average(Average_Mu_List_4GPa_40ms)
Mu5GPa_40ms = np.average(Average_Mu_List_5GPa_40ms)

MuAveragesDifferentPressures = np.array([Mu1GPa_40ms, Mu2GPa_40ms, Mu3GPa_40ms, Mu4GPa_40ms, Mu5GPa_40ms])
AverageMu_40ms = np.average(MuAveragesDifferentPressures)

ActivationEnergy_1GPa_40ms = (((1 * 1e9 * (np.average(Average_Mu_List_1GPa_40ms) * mean_actv_40ms) * 1e-30) - 1.381 * trend1GPa_40ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
ActivationEnergy_2GPa_40ms = (((2 * 1e9 * (np.average(Average_Mu_List_2GPa_40ms) * mean_actv_40ms) * 1e-30) - 1.381 * trend2GPa_40ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
ActivationEnergy_3GPa_40ms = (((3 * 1e9 * (np.average(Average_Mu_List_3GPa_40ms) * mean_actv_40ms) * 1e-30) - 1.381 * trend3GPa_40ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
ActivationEnergy_4GPa_40ms = (((4 * 1e9 * (np.average(Average_Mu_List_4GPa_40ms) * mean_actv_40ms) * 1e-30) - 1.381 * trend4GPa_40ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000
ActivationEnergy_5GPa_40ms = (((5 * 1e9 * (np.average(Average_Mu_List_5GPa_40ms) * mean_actv_40ms) * 1e-30) - 1.381 * trend5GPa_40ms[0] * 1e-23) * 6.02214076 * (10**23)) / 1000

lnA1GPa_40ms = np.log((np.exp(trend1GPa_40ms[1]) * (10 ** 9)))
lnA2GPa_40ms = np.log((np.exp(trend2GPa_40ms[1]) * (10 ** 9)))
lnA3GPa_40ms = np.log((np.exp(trend3GPa_40ms[1]) * (10 ** 9)))
lnA4GPa_40ms = np.log((np.exp(trend4GPa_40ms[1]) * (10 ** 9)))
lnA5GPa_40ms = np.log((np.exp(trend5GPa_40ms[1]) * (10 ** 9)))

print(f"ln(A) at 1GPa, 40ms is {lnA1GPa_40ms}")
print(f"ln(A) at 2GPa, 40ms is {lnA2GPa_40ms}")
print(f"ln(A) at 3GPa, 40ms is {lnA3GPa_40ms}")
print(f"ln(A) at 4GPa, 40ms is {lnA4GPa_40ms}")
print(f"ln(A) at 5GPa, 40ms is {lnA5GPa_40ms}")

print('Activation Energy 1GPa, 40ms =' + str(ActivationEnergy_1GPa_40ms))
print('Activation Energy 2GPa, 40ms =' + str(ActivationEnergy_2GPa_40ms))
print('Activation Energy 3GPa, 40ms =' + str(ActivationEnergy_3GPa_40ms))
print('Activation Energy 4GPa, 40ms =' + str(ActivationEnergy_4GPa_40ms))
print('Activation Energy 5GPa, 40ms =' + str(ActivationEnergy_5GPa_40ms))

def linear(x, m, n):
    return m * x + n

coef_p_cf_1GPa_40ms, coef_p_pcov_1GPa_40ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_1GPa_40ms)
sigma_A1GPa_40ms = coef_p_pcov_1GPa_40ms[1, 1] ** 0.5
sigma_m1GPa_40ms = coef_p_pcov_1GPa_40ms[0, 0] ** 0.5
dof_40ms = max(0, len(Log_Rates_1GPa_40ms) - len(coef_p_cf_1GPa_40ms))

alpha = 0.05
tval_40ms = t.ppf(1.0 - alpha / 2., 4.1532)

sigma_40ms = np.std(Activation_Volumes_40ms)
error_actv_40ms = sigma_40ms * tval_40ms / np.sqrt(len(Activation_Volumes_40ms))
uncert_A1GPa_40ms = sigma_A1GPa_40ms * tval
uncert_m1GPa_40ms = sigma_m1GPa_40ms * tval

coef_p_cf_2GPa_40ms, coef_p_pcov_2GPa_40ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_2GPa_40ms)
sigma_A2GPa_40ms = coef_p_pcov_2GPa_40ms[1, 1] ** 0.5
sigma_m2GPa_40ms = coef_p_pcov_2GPa_40ms[0, 0] ** 0.5
dof = max(0, len(Log_Rates_2GPa_40ms) - len(coef_p_cf_2GPa_40ms))
uncert_A2GPa_40ms = sigma_A2GPa_40ms * tval_40ms
uncert_m2GPa_40ms = sigma_m2GPa_40ms * tval_40ms

coef_p_cf_3GPa_40ms, coef_p_pcov_3GPa_40ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_3GPa_40ms)
sigma_A3GPa_40ms = coef_p_pcov_3GPa_40ms[1, 1] ** 0.5
sigma_m3GPa_40ms = coef_p_pcov_3GPa_40ms[0, 0] ** 0.5
dof = max(0, len(Log_Rates_3GPa_40ms) - len(coef_p_cf_3GPa_40ms))
uncert_A3GPa_40ms = sigma_A3GPa_40ms * tval_40ms
uncert_m3GPa_40ms = sigma_m3GPa_40ms * tval_40ms

coef_p_cf_4GPa_40ms, coef_p_pcov_4GPa_40ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_4GPa_40ms)
sigma_A4GPa_40ms = coef_p_pcov_4GPa_40ms[1, 1] ** 0.5
sigma_m4GPa_40ms = coef_p_pcov_4GPa_40ms[0, 0] ** 0.5
dof = max(0, len(Log_Rates_4GPa_40ms) - len(coef_p_cf_4GPa_40ms))
uncert_A4GPa_40ms = sigma_A4GPa_40ms * tval_40ms
uncert_m4GPa_40ms = sigma_m4GPa_40ms * tval_40ms

coef_p_cf_5GPa_40ms, coef_p_pcov_5GPa_40ms = optimize.curve_fit(linear, Inverse_Temperatures, Log_Rates_5GPa_40ms)
sigma_A5GPa_40ms = coef_p_pcov_5GPa_40ms[1, 1] ** 0.5
sigma_m5GPa_40ms = coef_p_pcov_5GPa_40ms[0, 0] ** 0.5
dof = max(0, len(Log_Rates_5GPa_40ms) - len(coef_p_cf_5GPa_40ms))
uncert_A5GPa_40ms = sigma_A5GPa_40ms * tval_40ms
uncert_m5GPa_40ms = sigma_m5GPa_40ms * tval_40ms

err_E_1GPa_40ms = (np.sqrt(((1 * 1e9 * np.average(Average_Mu_List_1GPa_40ms) * error_actv_40ms * 1e-30) ** 2 + (1.381 * uncert_m1GPa_40ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
err_E_2GPa_40ms = (np.sqrt(((2 * 1e9 * np.average(Average_Mu_List_2GPa_40ms) * error_actv_40ms * 1e-30) ** 2 + (1.381 * uncert_m2GPa_40ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
err_E_3GPa_40ms = (np.sqrt(((3 * 1e9 * np.average(Average_Mu_List_3GPa_40ms) * error_actv_40ms * 1e-30) ** 2 + (1.381 * uncert_m3GPa_40ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
err_E_4GPa_40ms = (np.sqrt(((4 * 1e9 * np.average(Average_Mu_List_4GPa_40ms) * error_actv_40ms * 1e-30) ** 2 + (1.381 * uncert_m4GPa_40ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000
err_E_5GPa_40ms = (np.sqrt(((5 * 1e9 * np.average(Average_Mu_List_5GPa_40ms) * error_actv_40ms * 1e-30) ** 2 + (1.381 * uncert_m5GPa_40ms * 1e-23) ** 2)) * (6.02214076 * (10**23))) / 1000

print(f"Error for Activation Energy at 1GPa = {err_E_1GPa_40ms}")
print(f"Error for Activation Energy at 2GPa = {err_E_2GPa_40ms}")
print(f"Error for Activation Energy at 3GPa = {err_E_3GPa_40ms}")
print(f"Error for Activation Energy at 4GPa = {err_E_4GPa_40ms}")
print(f"Error for Activation Energy at 5GPa = {err_E_5GPa_40ms}")

uncertlnA1GPa_40ms = np.log(uncert_A1GPa_40ms * (10 ** 9))
uncertlnA2GPa_40ms = np.log(uncert_A2GPa_40ms * (10 ** 9))
uncertlnA3GPa_40ms = np.log(uncert_A3GPa_40ms * (10 ** 9))
uncertlnA4GPa_40ms = np.log(uncert_A4GPa_40ms * (10 ** 9))
uncertlnA5GPa_40ms = np.log(uncert_A5GPa_40ms * (10 ** 9))

print(f"Error for ln A at 1GPa = {uncert_A1GPa_40ms}")
print(f"Error for ln A at 2GPa = {uncert_A2GPa_40ms}")
print(f"Error for ln A at 3GPa = {uncert_A3GPa_40ms}")
print(f"Error for ln A at 4GPa = {uncert_A4GPa_40ms}")
print(f"Error for ln A at 5GPa = {uncert_A5GPa_40ms}")

pressures = ['2GPa', '3GPa', '4GPa', '5GPa']

########################## Plotting ln(Rates) vs Normal Stress #####################################
x = np.array([0, 1, 2, 3, 4, 5])
params400K_40ms = np.polyfit(NormalStressMeans_400K_40ms, Log_Rates_400K_40ms, 1)
params500K_40ms = np.polyfit(NormalStressMeans_500K_40ms, Log_Rates_500K_40ms, 1)
params600K_40ms = np.polyfit(NormalStressMeans_600K_40ms, Log_Rates_600K_40ms, 1)
params700K_40ms = np.polyfit(NormalStressMeans_700K_40ms, Log_Rates_700K_40ms, 1)

RatesvsNormal, RvN3 = plt.subplots()
RvN3.set_title('Log of Dissociation Rates vs Normal Stress - Alpha Fe, 40ms')
RvN3.set_xlabel('Normal Stress(GPa)')
RvN3.set_ylabel('Log of Dissociation Rate (ns-1)')
RvN3.scatter(NormalStressMeans_400K_40ms, Log_Rates_400K_40ms)
RvN3.scatter(NormalStressMeans_500K_40ms, Log_Rates_500K_40ms)
RvN3.scatter(NormalStressMeans_600K_40ms, Log_Rates_600K_40ms)
RvN3.scatter(NormalStressMeans_700K_40ms, Log_Rates_700K_40ms)
RvN3.plot(x, params400K_40ms[0] * x + params400K_40ms[1], label='400K Fitted')
RvN3.plot(x, params500K_40ms[0] * x + params500K_40ms[1], label='500K Fitted')
RvN3.plot(x, params600K_40ms[0] * x + params600K_40ms[1], label='600K Fitted')
RvN3.plot(x, params700K_40ms[0] * x + params700K_40ms[1], label='700K Fitted')
#RvN3.set_xlim(1, 5)
# RvN3.set_ylim(0, 25)
RvN3.legend(loc='lower right')
plt.show()

def function(data, a, b, c):
    x = data[0]
    y = data[1]
    return a * (x**b) * (y**c) #TODO change fitting function

x_data = []
y_data = []
z_data = []

data = [[400, Average_Shear_Stress_List_400K_40ms[0], LogRate_40ms_400K_1GPa], [400, Average_Shear_Stress_List_400K_40ms[1], LogRate_40ms_400K_2GPa], [400, Average_Shear_Stress_List_400K_40ms[2], LogRate_40ms_400K_3GPa], [400, Average_Shear_Stress_List_400K_40ms[3], LogRate_40ms_400K_4GPa], [400, Average_Shear_Stress_List_400K_40ms[4], LogRate_40ms_400K_5GPa],
        [500, Average_Shear_Stress_List_500K_40ms[0], LogRate_40ms_500K_1GPa], [500, Average_Shear_Stress_List_500K_40ms[1], LogRate_40ms_500K_2GPa], [500, Average_Shear_Stress_List_500K_40ms[2], LogRate_40ms_500K_3GPa], [500, Average_Shear_Stress_List_500K_40ms[3], LogRate_40ms_500K_4GPa], [500, Average_Shear_Stress_List_500K_40ms[4], LogRate_40ms_500K_5GPa],
        [600, Average_Shear_Stress_List_600K_40ms[0], LogRate_40ms_600K_1GPa], [600, Average_Shear_Stress_List_600K_40ms[1], LogRate_40ms_600K_2GPa], [600, Average_Shear_Stress_List_600K_40ms[2], LogRate_40ms_600K_3GPa], [600, Average_Shear_Stress_List_600K_40ms[3], LogRate_40ms_600K_4GPa], [600, Average_Shear_Stress_List_600K_40ms[4], LogRate_40ms_600K_5GPa],
        [700, Average_Shear_Stress_List_700K_40ms[0], LogRate_40ms_700K_1GPa], [700, Average_Shear_Stress_List_700K_40ms[1], LogRate_40ms_700K_2GPa], [700, Average_Shear_Stress_List_700K_40ms[2], LogRate_40ms_700K_3GPa], [700, Average_Shear_Stress_List_700K_40ms[3], LogRate_40ms_700K_4GPa], [700, Average_Shear_Stress_List_700K_40ms[4], LogRate_40ms_700K_5GPa]]

for item in data:
    x_data.append(item[0])
    y_data.append(item[1])
    z_data.append(item[2])
#
#
parameters, covariance = optimize.curve_fit(function, [x_data, y_data], z_data)

# create surface function model
# setup data points for calculating surface model
model_x_data = np.linspace(min(x_data), max(x_data), 30)
model_y_data = np.linspace(min(y_data), max(y_data), 30)

# create coordinate arrays for vectorized evaluations
X, Y = np.meshgrid(model_x_data, model_y_data)
# calculate Z coordinate array
Z = function(np.array([X, Y]), *parameters)

z = []
for row in Z:
    row.sort()
    z.append(row)

zlogs = np.array(z)
import matplotlib.cm
cm = plt.get_cmap("jet")
fig = plt.figure()
ax4 = plt.axes(projection='3d')
ax4.plot_surface(X, Y, Z, cmap=cm, alpha=0.5, edgecolor='black', linewidth=0.3)
ax4.set_title('3D Plot - Variation in Log of Dissociation Rates - Alpha Fe, 40ms')
ax4.scatter(x_data, y_data, z_data, color='black', alpha=1)
ax4.set_xlabel('Temperature (K)')
ax4.invert_xaxis()
ax4.set_ylabel('Shear Stress (GPa)')
ax4.set_zlabel('Log of Dissociation Rate (per ns)')
plt.show()
plt.close(fig)

x_data = []
y_data = []
z_data = []

data = [[400, Average_Shear_Stress_List_400K_40ms[0], Dissociation_Rate_40ms_400K_1GPa], [400, Average_Shear_Stress_List_400K_40ms[1], Dissociation_Rate_40ms_400K_2GPa], [400, Average_Shear_Stress_List_400K_40ms[2], Dissociation_Rate_40ms_400K_3GPa], [400, Average_Shear_Stress_List_400K_40ms[3], Dissociation_Rate_40ms_400K_4GPa], [400, Average_Shear_Stress_List_400K_40ms[4], Dissociation_Rate_40ms_400K_5GPa],
        [500, Average_Shear_Stress_List_500K_40ms[0], Dissociation_Rate_40ms_500K_1GPa], [500, Average_Shear_Stress_List_500K_40ms[1], Dissociation_Rate_40ms_500K_2GPa], [500, Average_Shear_Stress_List_500K_40ms[2], Dissociation_Rate_40ms_500K_3GPa], [500, Average_Shear_Stress_List_500K_40ms[3], Dissociation_Rate_40ms_500K_4GPa], [500, Average_Shear_Stress_List_500K_40ms[4], Dissociation_Rate_40ms_500K_5GPa],
        [600, Average_Shear_Stress_List_600K_40ms[0], Dissociation_Rate_40ms_600K_1GPa], [600, Average_Shear_Stress_List_600K_40ms[1], Dissociation_Rate_40ms_600K_2GPa], [600, Average_Shear_Stress_List_600K_40ms[2], Dissociation_Rate_40ms_600K_3GPa], [600, Average_Shear_Stress_List_600K_40ms[3], Dissociation_Rate_40ms_600K_4GPa], [600, Average_Shear_Stress_List_600K_40ms[4], Dissociation_Rate_40ms_600K_5GPa],
        [700, Average_Shear_Stress_List_700K_40ms[0], Dissociation_Rate_40ms_700K_1GPa], [700, Average_Shear_Stress_List_700K_40ms[1], Dissociation_Rate_40ms_700K_2GPa], [700, Average_Shear_Stress_List_700K_40ms[2], Dissociation_Rate_40ms_700K_3GPa], [700, Average_Shear_Stress_List_700K_40ms[3], Dissociation_Rate_40ms_700K_4GPa], [700, Average_Shear_Stress_List_700K_40ms[4], Dissociation_Rate_40ms_700K_5GPa]]

for item in data:
    x_data.append(item[0])
    y_data.append(item[1])
    z_data.append(item[2])

parameters, covariance = optimize.curve_fit(function, [x_data, y_data], z_data)

# create surface function model
# setup data points for calculating surface model
model_x_data = np.linspace(min(x_data), max(x_data), 30)
model_y_data = np.linspace(min(y_data), max(y_data), 30)
# create coordinate arrays for vectorized evaluations
X, Y = np.meshgrid(model_x_data, model_y_data)
# calculate Z coordinate array
Z = function(np.array([X, Y]), *parameters)

z = []
for row in Z:
    row.sort()
    z.append(row)

z = np.array(z)
# print(Z)
# print('####################')
# print(z)

cm = plt.get_cmap("jet")
fig = plt.figure()
ax4 = plt.axes(projection='3d')
ax4.plot_surface(X, Y, Z, cmap=cm, alpha=0.5, edgecolor= 'black', linewidth=0.3)
ax4.set_title('3D Plot - Variation in Dissociation Rates - Alpha Fe, 40ms')
ax4.scatter(x_data, y_data, z_data, color='black', alpha=1)
ax4.set_xlabel('Temperature (K)')
ax4.invert_xaxis()
ax4.set_ylabel('Shear Stress (GPa)')
ax4.set_zlabel('Dissociation Rate (per ns)')
plt.show()

# ########## CALCULATING THE ACTIVATION ENERGY, VOLUME AND PREFACTOR FROM 3D FIT ##########
#
logRates3D_40ms = zlogs
shearstresses40ms = model_y_data
Index = 0
sigma = []
params = []
while Index < len(shearstresses40ms) - 1:
    for logRates3Drow in logRates3D_40ms:
        coef_sh_cf, coef_sh_pcov = optimize.curve_fit(linear, shearstresses40ms, logRates3Drow)
        paramsrow = np.polyfit(shearstresses40ms, logRates3Drow, 1)
        sigmarow = coef_sh_pcov[0, 0] ** 0.5
        sigma.append(sigmarow)
        params.append(paramsrow)
        Index +=1

#
sigma = np.array(sigma)
params = np.array(params)
params = np.average(params)
sigma = np.average(sigma)

activation_vol_3D_40ms = (params) * (1.38065) * 550 * 1e-2
print(f'3D Activation Volume, 40ms is {activation_vol_3D_40ms}')

alpha = 0.05
sigma = sigma ** 0.5
dof = 2
tval = t.ppf(1.0 - alpha / 2., dof)
uncert = sigma * tval
uncert_ActivationVolume_40ms = uncert * (1.38065) * 550 * 1e-2
print(f'Activation Volume uncertainty for 3D fit at 40ms is {uncert_ActivationVolume_40ms}')

logRates3D_40ms = zlogs
temperatures = model_x_data

inverse_temperatures = [1 / x for x in temperatures]

Index = 0
sigma = []
SigmaA = []
params = []
interceptaverage = []
while Index < len(inverse_temperatures) - 1:
    for logRates3Drow in logRates3D_40ms:
        coef_sh_cf, coef_sh_pcov = optimize.curve_fit(linear, inverse_temperatures, logRates3Drow)
        paramsrow = np.polyfit(inverse_temperatures, logRates3Drow, 1)
        sigmarow = coef_sh_pcov[0, 0] ** 0.5
        sigma_A = coef_sh_pcov[1, 1] ** 0.5
        sigma.append(sigmarow)
        SigmaA.append(sigma_A)
        params.append(paramsrow[0])
        intercept = paramsrow[1]
        interceptaverage.append((intercept))
        Index +=1

sigma = np.array(sigma)
params = np.array(params)
SigmaA = np.array(SigmaA)
interceptaverage = np.array(interceptaverage)
params_40ms = np.average(params)
sigma = np.average(sigma)
interceptaverage = np.average(interceptaverage)
alpha = 0.05
sigma = sigma ** 0.5
tval = t.ppf(1.0 - alpha / 2., dof)
uncert = sigma * tval

Mu1GPa_40ms = np.average(Average_Mu_List_1GPa_40ms)
Mu2GPa_40ms = np.average(Average_Mu_List_2GPa_40ms)
Mu3GPa_40ms = np.average(Average_Mu_List_3GPa_40ms)
Mu4GPa_40ms = np.average(Average_Mu_List_4GPa_40ms)
Mu5GPa_40ms = np.average(Average_Mu_List_5GPa_40ms)

MuAveragesDifferentPressures = np.array([Mu1GPa_40ms, Mu2GPa_40ms, Mu3GPa_40ms, Mu4GPa_40ms, Mu5GPa_40ms])
AverageMu_40ms = np.average(MuAveragesDifferentPressures)

ActivationEnergy_3D_40ms = (((3 * 1e9 * (AverageMu_40ms * activation_vol_3D_40ms) * 1e-30) - 1.381 * params_40ms * 1e-23) * 6.02214076 * (10**23)) / 1000

print(f'Activation Energy for 3D fit at 40ms is {ActivationEnergy_3D_40ms}')

error_3D_40ms = (np.sqrt((3 * 1e9 * np.average(AverageMu_40ms) * uncert_ActivationVolume_40ms * 1e-30) ** 2 + (1.381 * uncert * 1e-23) ** 2) * 6.02214076 * (10**23)) / 1000
print(f"Activation_Energy Error at 40ms is {error_3D_40ms}")

uncert_prefactor_3D_40ms = sigma_A * tval
lnA_3D_40ms = np.log((np.exp(interceptaverage) * (10 ** 9)))

print(f"ln(A) for 3D fit is {lnA_3D_40ms}")
print(f"ln(A) uncertainty is {uncert_prefactor_3D_40ms}")

######### Plotting 3D fit vs 2D fit results along with their error margins ##########
"""
Need to get a list with:
- Values for Activation Volumes at different temperatures
- Values for Activation Energies at different pressures
- Values for ln(A) at different pressures
- List of errors for each of the above quantities
- Make a graph for each surface chemistry
- Eventually will need to do the same for the different sliding speeds

"""
Temperatures = [400, 500, 600, 700]
Pressures = [1, 2, 3, 4, 5]


Activation_Energies_40ms = [ActivationEnergy_1GPa_40ms, ActivationEnergy_2GPa_40ms, ActivationEnergy_3GPa_40ms, ActivationEnergy_4GPa_40ms, ActivationEnergy_5GPa_40ms]
Activation_Energy_Errors_40ms = [err_E_1GPa_40ms, err_E_2GPa_40ms, err_E_3GPa_40ms, err_E_4GPa_40ms, err_E_5GPa_40ms]
ActivationEnergy_3D_40ms = ActivationEnergy_3D_40ms
ActivationEnergy_3D_error_40ms = error_3D_40ms
ActivationEnergy_3D_error_UpperBound_Value_40ms = float(float(ActivationEnergy_3D_40ms) + float(ActivationEnergy_3D_error_40ms))
ActivationEnergy_3D_error_LowerBound_Value_40ms = float(float(ActivationEnergy_3D_40ms) - float(ActivationEnergy_3D_error_40ms))

Activation_Energy_Error_Plot, Ea2Dvs3d  = plt.subplots()
Ea2Dvs3d.set_title('Comparison of Activation Energies from 2D and 3D Fits - AlphaFe, 40ms')
Ea2Dvs3d.set_xlabel('Normal Stress (GPa)')
Ea2Dvs3d.set_ylabel('Activation Energy')
Ea2Dvs3d.scatter(Pressures, Activation_Energies_40ms)
Ea2Dvs3d.errorbar(Pressures, Activation_Energies_40ms, yerr=Activation_Energy_Errors_40ms, linestyle="None", fmt='o', capsize=3)
Ea2Dvs3d.axhline(y=ActivationEnergy_3D_40ms)
Pressures = [0.5, 1, 2, 3, 4, 5, 5.5]
Ea2Dvs3d.fill_between(Pressures, ActivationEnergy_3D_error_LowerBound_Value_40ms, ActivationEnergy_3D_error_UpperBound_Value_40ms, alpha=0.4)
Ea2Dvs3d.set_xlim(0.5, 5.5)
Ea2Dvs3d.set_ylim(0, 35)
plt.show()

################ Activation Volume Errors #############################

Activation_Volumes_40ms = [activation_vol_400K_40ms, activation_vol_500K_40ms, activation_vol_600K_40ms, activation_vol_700K_40ms]
Activation_Volume_Errors_40ms = [uncert400_40ms, uncert500_40ms, uncert600_40ms, uncert700_40ms]
Activation_Volume_3D_40ms = activation_vol_3D_40ms
Activation_Volume_3D_Error_40ms = uncert_ActivationVolume_40ms
ActivationVolume_3D_error_UpperBound_Value_40ms = float(float(Activation_Volume_3D_40ms) + float(Activation_Volume_3D_Error_40ms))
ActivationVolume_3D_error_LowerBound_Value_40ms = float(float(Activation_Volume_3D_40ms) - float(Activation_Volume_3D_Error_40ms))

Activation_Volume_Error_Plot, Av2Dvs3d  = plt.subplots()
Av2Dvs3d.set_title('Comparison of Activation Volumes from 2D and 3D Fits - AlphaFe, 40ms')
Av2Dvs3d.set_xlabel('Normal Stress(GPa)')
Av2Dvs3d.set_ylabel('Activation Volume')
Av2Dvs3d.scatter(Temperatures, Activation_Volumes_40ms)
Av2Dvs3d.errorbar(Temperatures, Activation_Volumes_40ms, yerr=Activation_Volume_Errors_40ms, linestyle="None", fmt='o', capsize=3)
Av2Dvs3d.axhline(y=Activation_Volume_3D_40ms)
Temperatures = [350, 400, 500, 600, 700, 750]
Av2Dvs3d.fill_between(Temperatures, ActivationVolume_3D_error_LowerBound_Value_40ms, ActivationVolume_3D_error_UpperBound_Value_40ms, alpha=0.4)
Av2Dvs3d.set_xlim(350, 750)
Av2Dvs3d.set_ylim(0, 30)
plt.show()

#################### Prefactor Errors #########################
Pressures = [1, 2, 3, 4, 5]
Prefactors_40ms = [lnA1GPa_40ms, lnA2GPa_40ms, lnA3GPa_40ms, lnA4GPa_40ms, lnA5GPa_40ms]
PrefactorErrors_40ms = [uncert_A1GPa_40ms, uncert_A2GPa_40ms, uncert_A3GPa_40ms, uncert_A4GPa_40ms, uncert_A5GPa_40ms]
lnA_3D_error_40ms = uncert_prefactor_3D_40ms
lnA_3D_error_UpperBound_Value_40ms = float(float(lnA_3D_40ms) + float(lnA_3D_error_40ms))
lnA_3D_error_LowerBound_Value_40ms = float(float(lnA_3D_40ms) - float(lnA_3D_error_40ms))

Prefactor_Error_Plot, lnA2Dvs3d  = plt.subplots()
lnA2Dvs3d.set_title('Comparison of Prefactors from 2D and 3D Fits - AlphaFe, 40ms')
lnA2Dvs3d.set_xlabel('Normal Stress(GPa)')
lnA2Dvs3d.set_ylabel('Prefactor')
lnA2Dvs3d.scatter(Pressures, Prefactors_40ms)
lnA2Dvs3d.errorbar(Pressures, Prefactors_40ms, yerr=PrefactorErrors_40ms, linestyle="None", fmt='o', capsize=3)
lnA2Dvs3d.axhline(y=lnA_3D_40ms)
Pressures = [0.5, 1, 2, 3, 4, 5, 5.5]
lnA2Dvs3d.fill_between(Pressures, lnA_3D_error_LowerBound_Value_40ms, lnA_3D_error_UpperBound_Value_40ms, alpha=0.4)
lnA2Dvs3d.set_xlim(0.5, 5.5)
lnA2Dvs3d.set_ylim(22, 28)
plt.show()

########## Checking for Kinetic Compensation Effect ##########

KineticCompeEffectPlot, EavsLnA  = plt.subplots()
EavsLnA.set_title('Activation Energy vs Prefactor - AlphaFe, 40ms')
EavsLnA.set_xlabel('Ea')
EavsLnA.set_ylabel('ln A')
EavsLnA.scatter(Activation_Energies_40ms, Prefactors_40ms)
EavsLnA.set_ylim(24, 26.5)
plt.show()

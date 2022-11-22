"""
Functions to help process TCP Data AT DIFFERENT SLIDING SPEEDS
Functions Present:
- Getting intact columns from processed experiments
- Plotting fitted exponential graphs of decay of intact molecules
- Get friction values from fc.ave.dump file
- Processing fc_ave.dump file to get average shear/normal stress and average mu

"""

import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy import optimize
import numpy as np
import statistics
import os
import os.path
import sys
import glob
import sys
from numpy import mean as m

def get_intact_columns_constant_temperature_and_speed(Temperature, Pressures, Speed):
    """
    speed: Speed directory that you ran the experiments at
    pressures: Different pressures that you ran the experiments at with constant speed
    :return: Dataframe with concatenated and renamed columns from each of the experiments at the input pressures and speed
    """

    #Makes a dataframe from the intact molecuels csv of the first speed/temperature/pressure experiment in your defined list
    All_Dataframe = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{Speed}/{Temperature}/1GPa/'
                                'processed/Fe2O3-200-iso-octane_IntactMols_0.3'.format(Temperature=Temperature, Speed=Speed), sep='\t')

    #Take the correct columns from the above datafram using iloc
    Big_Dataframe_NotProcessed = All_Dataframe.iloc[:, [0, 2]]

    #Rename the column to the pressure of the first experiment
    Big_Dataframe = Big_Dataframe_NotProcessed.rename(columns={'Intact_molecules_noomit' : '1GPa'})

    #Quick for loop to get the rest of the columns from the remaining

    Pressures = Pressures
    for P in Pressures:
        Dataframe = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{Speed}/{Temperature}/{P}/'
                                'processed/Fe2O3-200-iso-octane_IntactMols_0.3'.format(Temperature=Temperature, Speed=Speed, P=P), sep='\t')
        UnnamedBig_DataframeP = Dataframe.iloc[:, [0, 2]]
        #print(Big_DataframeP)
        Big_DataframeP = UnnamedBig_DataframeP.rename(columns= {'Timestep' : 'Timestep_{}'.format(P),
                                                                'Intact_molecules_noomit' : '{}'.format(P)})
        Big_Dataframe = pd.concat([Big_Dataframe, Big_DataframeP], axis =1)

    # Using .dropna from pandas library to get rid of rows with NaNa (which mean the simulation ran out of time
    Big_Dataframe = Big_Dataframe.dropna(axis=0)

    return Big_Dataframe

def get_intact_columns_constant_pressure_and_speed(Speed, Pressure, Temperatures):

    """
    speed: Speed directory that you ran the experiments at
    :pressures: Different pressures that you ran the experiments at with constant speed
    :return: Dataframe with concatenated and renamed columns from each of the experiments at the input pressures and speed
    """

    #Makes a dataframe from the intact molecuels csv of the first speed/temperature/pressure experiment in your defined list
    All_Dataframe = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{Speed}/400K/{Pressure}/'
                                'processed/Fe2O3-200-iso-octane_IntactMols_0.3'.format(Speed=Speed, Pressure=Pressure), sep='\t')

    #Take the correct columns from the above datafram using iloc
    Big_Dataframe_NotProcessed = All_Dataframe.iloc[:, [0, 2]]

    #Rename the column to the pressure of the first experiment
    Big_Dataframe = Big_Dataframe_NotProcessed.rename(columns={'Intact_molecules_noomit' : '400K'})

    #Quick for loop to get the rest of the columns from the remaining

    for T in Temperatures:
        Dataframe = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{Speed}/{T}/{Pressure}/'
                                'processed/Fe2O3-200-iso-octane_IntactMols_0.3'.format(Speed=Speed, T=T, Pressure=Pressure), sep='\t')
        UnnamedBig_DataframeP = Dataframe.iloc[:, [0, 2]]
        Big_DataframeP = UnnamedBig_DataframeP.rename(columns= {'Timestep' : 'Timestep_{}'.format(T),
                                                                'Intact_molecules_noomit' : '{}'.format(T)})
        Big_Dataframe = pd.concat([Big_Dataframe, Big_DataframeP], axis =1)

    # Using .dropna from pandas library to get rid of rows with NaNa (which mean the simulation ran out of time
    Big_Dataframe = Big_Dataframe.dropna(axis=0)

    return Big_Dataframe

def get_intact_columns_constant_pressure_and_temperature(Pressure, Speeds, Temperature):

    """
    speed: Speed directory that you ran the experiments at
    :pressures: Different pressures that you ran the experiments at with constant speed
    :return: Dataframe with concatenated and renamed columns from each of the experiments at the input temperatures
    """

    # Makes a dataframe from the intact molecuels csv of the first speed/temperature/pressure experiment in your defined list
    All_Dataframe = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/10ms/{Temperature}/{Pressure}/'
                                'processed/Fe2O3-200-iso-octane_IntactMols_0.3'.format(Pressure=Pressure, Temperature=Temperature), sep='\t')

    # Take the correct columns from the above datafram using iloc
    Big_Dataframe_NotProcessed = All_Dataframe.iloc[:, [0, 2]]

    # Rename the column to the pressure of the first experiment
    Big_Dataframe = Big_Dataframe_NotProcessed.rename(columns={'Intact_molecules_noomit': '10ms'})

    # Quick for loop to get the rest of the columns from the remaining

    for S in Speeds:
        Dataframe = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{S}/{Temperature}/{Pressure}/'
                                'processed/Fe2O3-200-iso-octane_IntactMols_0.3'.format(S=S, Temperature=Temperature, Pressure=Pressure), sep='\t')
        UnnamedBig_DataframeP = Dataframe.iloc[:, [0, 2]]
        # print(Big_DataframeP)
        Big_DataframeP = UnnamedBig_DataframeP.rename(columns={'Timestep': 'Timestep_{}'.format(S),
                                                               'Intact_molecules_noomit': '{}'.format(S)})
        Big_Dataframe = pd.concat([Big_Dataframe, Big_DataframeP], axis=1)

    # Using .dropna from pandas library to get rid of rows with NaNa (which mean the simulation ran out of time
    Big_Dataframe = Big_Dataframe.dropna(axis=0)

    return Big_Dataframe

def get_average_shear_normal_stress_and_average_mu_constant_temperature(Temperature, Pressures):
    Friction_Coefficient_Dataframe_Unnamed = pd.read_csv('C:/Users/eeo21/Documents/PhD/TCPDecompositionExperiments/Completed/AlphaFe/{}/1GPa/'
                                'fc_ave.dump'.format(Temperature), sep=' ')
    Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe_Unnamed.rename(columns={'v_s_bot' : 'Shear Stress 1GPa', 'v_p_bot' : 'Normal Stress 1GPa'})

    for P in Pressures:
        Dataframe = pd.read_csv('C:/Users/eeo21/Documents/PhD/TCPDecompositionExperiments/Completed/AlphaFe/{}/{}/'
                                'fc_ave.dump'.format(Temperature, P), sep=' ')
        Big_DataframeP = Dataframe.rename(columns= {'Timestep': 'Timestep {}'.format(P),
                                                        'v_s_bot': 'Shear Stress {}'.format(P),
                                                        'v_p_bot': 'Normal Stress {}'.format(P)})

        Friction_Coefficient_Dataframe = pd.concat([Friction_Coefficient_Dataframe, Big_DataframeP], axis =1)
        Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe.dropna()


    #print(Friction_Coefficient_Dataframe)
    Mu_Final_Dataframe = Friction_Coefficient_Dataframe.iloc[:, [0, 1, 2, 4, 5, 7, 8, 10, 11, 13, 14]]
    #print(Mu_Final_Dataframe)

    ShearStressMeans = Mu_Final_Dataframe[['Shear Stress 1GPa', 'Shear Stress 2GPa', 'Shear Stress 3GPa', 'Shear Stress 4GPa', 'Shear Stress 5GPa']].mean()
    Average_Shear_Stress_Dictionary = ShearStressMeans.to_dict()
    #print(ShearStressMeans)
    NormalStressMeans = Mu_Final_Dataframe[['Normal Stress 1GPa', 'Normal Stress 2GPa', 'Normal Stress 3GPa', 'Normal Stress 4GPa', 'Normal Stress 5GPa']].mean()
    NormalStressMeans = NormalStressMeans.to_dict()
    #print(NormalStressMeans)


    Average_Mu_Dictionary = {}

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
    Average_Mu_List = list(Average_Mu_Dictionary.values())
    Average_Shear_Stress_List = [x / 10000 for x in Average_Shear_Stress_List]

    return Average_Shear_Stress_List, Average_Mu_List, NormalStressMeans

def get_average_shear_normal_stress_and_average_mu_constant_pressure(Pressure, Temperatures):
    Friction_Coefficient_Dataframe_Unnamed = pd.read_csv('C:/Users/eeo21/Documents/PhD/TCPDecompositionExperiments/Completed/AlphaFe/300K/{}/'
                                'fc_ave.dump'.format(Pressure), sep=' ')
    Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe_Unnamed.rename(columns={'v_s_bot' : 'Shear Stress 300K', 'v_p_bot' : 'Normal Stress 300K'})

    for T in Temperatures:
        Dataframe = pd.read_csv('C:/Users/eeo21/Documents/PhD/TCPDecompositionExperiments/Completed/AlphaFe/{}/{}/'
                                'fc_ave.dump'.format(T, Pressure), sep=' ')
        Big_DataframeP = Dataframe.rename(columns= {'Timestep': 'Timestep {}'.format(T),
                                                        'v_s_bot': 'Shear Stress {}'.format(T),
                                                        'v_p_bot': 'Normal Stress {}'.format(T)})

        Friction_Coefficient_Dataframe = pd.concat([Friction_Coefficient_Dataframe, Big_DataframeP], axis =1)
        Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe.dropna()


    #print(Friction_Coefficient_Dataframe)
    Mu_Final_Dataframe = Friction_Coefficient_Dataframe.iloc[:, [0, 1, 2, 4, 5, 7, 8, 10, 11, 13, 14]]

    ShearStressMeans = Mu_Final_Dataframe[['Shear Stress 300K', 'Shear Stress 400K', 'Shear Stress 500K', 'Shear Stress 600K', 'Shear Stress 700K']].mean()
    Average_Shear_Stress_Dictionary = ShearStressMeans.to_dict()
    #print(ShearStressMeans)
    NormalStressMeans = Mu_Final_Dataframe[['Normal Stress 300K', 'Normal Stress 400K', 'Normal Stress 500K', 'Normal Stress 600K', 'Normal Stress 700K']].mean()
    NormalStressMeans = NormalStressMeans.to_dict()
    #print(NormalStressMeans)

    Average_Mu_Dictionary = {}

    Normal_Stress = NormalStressMeans.get('Normal Stress 300K')
    #print(Normal_Stress)
    Shear_Stress = ShearStressMeans.get('Shear Stress 300K')
    #print(Shear_Stress)
    Average_Mu = Shear_Stress / Normal_Stress
    Average_Mu_Dictionary.update({'Average Mu 300K': Average_Mu})

    for T in Temperatures:

        Normal_Stress = NormalStressMeans.get('Normal Stress {}'.format(T))
        Shear_Stress = ShearStressMeans.get('Shear Stress {}'.format(T))
        Average_Mu = Shear_Stress / Normal_Stress
        Average_Mu_Dictionary.update({'Average Mu {}'.format(T): Average_Mu})

    Average_Shear_Stress_List = list(Average_Shear_Stress_Dictionary.values())
    Average_Mu_List = list(Average_Mu_Dictionary.values())
    Average_Shear_Stress_List = [x / 10000 for x in Average_Shear_Stress_List]

    return Average_Shear_Stress_List, Average_Mu_List, NormalStressMeans

def plot_shear_stress_vs_normal_stress(Average_Shear_Stress_List_1, Average_Shear_Stress_List_2, Average_Shear_Stress_List_3, Average_Shear_Stress_List_4, Temp1, Temp2, Temp3, Temp4, Speed):
    x = np.array([1, 2, 3, 4, 5])
    a, b = np.polyfit(x, Average_Shear_Stress_List_1, 1)
    c, d = np.polyfit(x, Average_Shear_Stress_List_2, 1)
    e, f = np.polyfit(x, Average_Shear_Stress_List_3, 1)
    g, h = np.polyfit(x, Average_Shear_Stress_List_4, 1)

    fig1, ax2 = plt.subplots()
    ax2.set_title(f'Shear Stress vs Normal Stress at Different Temperatures {Speed}')
    ax2.set_xlabel('Normal Stress (GPa)')
    ax2.set_ylabel('Shear Stress (GPa)')
    ax2.scatter(x, Average_Shear_Stress_List_1)
    ax2.scatter(x, Average_Shear_Stress_List_2)
    ax2.scatter(x, Average_Shear_Stress_List_3)
    ax2.scatter(x, Average_Shear_Stress_List_4)
    ax2.plot(x, a * x + b, label=Temp1)
    ax2.plot(x, c * x + d, label=Temp2)
    ax2.plot(x, e * x + f, label=Temp3)
    ax2.plot(x, g * x + h, label=Temp4)
    ax2.legend()
    plt.show()

def plot_variation_in_mu(Average_Mu_List_1, Average_Mu_List_2, Average_Mu_List_3, Average_Mu_List_4, Average_Mu_List_5, temp1, temp2, temp3, temp4, temp5):
    x = np.array([1, 2, 3, 4, 5])
    a, b = np.polyfit(x, Average_Mu_List_1, 1)
    c, d = np.polyfit(x, Average_Mu_List_2, 1)
    e, f = np.polyfit(x, Average_Mu_List_3, 1)
    g, h = np.polyfit(x, Average_Mu_List_4, 1)
    i, j = np.polyfit(x, Average_Mu_List_5, 1)

    fig1, ax1 = plt.subplots()
    ax1.set_title('Normal Stress vs Mu at Different Speeds')
    ax1.set_xlabel('Normal Stress (GPa)')
    ax1.set_ylabel('Shear Stress (GPa)')
    ax1.scatter(x, Average_Mu_List_1)
    ax1.scatter(x, Average_Mu_List_2)
    ax1.scatter(x, Average_Mu_List_3)
    ax1.scatter(x, Average_Mu_List_4)
    ax1.scatter(x, Average_Mu_List_5)
    ax1.plot(x, a * x + b, label=temp1)
    ax1.plot(x, c * x + d, label=temp2)
    ax1.plot(x, e * x + f, label=temp3)
    ax1.plot(x, g * x + h, label=temp4)
    ax1.plot(x, i * x + j, label=temp5)
    ax1.legend()
    plt.show()

def get_dissociation_rates(Timestep_List, fitted_functionOneGPa, fitted_functionTwoGPa, fitted_functionThreeGPa,
                           fitted_functionFourGPa, fitted_functionFiveGPa, cutoff):
    Dissociation_rates = []

    Time = []
    for x in Timestep_List:
        Time.append(((x - 400000) / int(4 * 10 ** 6)))
    Value = list(range(0, len(Time)))

    Value_Delta_List_1GPa = [(fitted_functionOneGPa[x + 1] - fitted_functionOneGPa[x]) for x in Value if
                             fitted_functionOneGPa[x] > cutoff]
    Value_Delta_List_Mean_1GPa = statistics.fmean(Value_Delta_List_1GPa)
    Rate_Per_ns_1GPa = Value_Delta_List_Mean_1GPa * 1000
    Dissociation_rates.append(Rate_Per_ns_1GPa)

    Value_Delta_List_2GPa = [(fitted_functionTwoGPa[x + 1] - fitted_functionTwoGPa[x]) for x in Value if
                             fitted_functionTwoGPa[x] > cutoff]
    Value_Delta_List_Mean_2GPa = statistics.fmean(Value_Delta_List_2GPa)
    Rate_Per_ns_2GPa = Value_Delta_List_Mean_2GPa * 1000
    Dissociation_rates.append(Rate_Per_ns_2GPa)

    Value_Delta_List_3GPa = [(fitted_functionThreeGPa[x + 1] - fitted_functionThreeGPa[x]) for x in Value if
                             fitted_functionThreeGPa[x] > cutoff]
    Value_Delta_List_Mean_3GPa = statistics.fmean(Value_Delta_List_3GPa)
    Rate_Per_ns_3GPa = Value_Delta_List_Mean_3GPa * 1000
    Dissociation_rates.append(Rate_Per_ns_3GPa)

    Value_Delta_List_4GPa = [(fitted_functionFourGPa[x + 1] - fitted_functionFourGPa[x]) for x in Value if
                             fitted_functionFourGPa[x] > cutoff]
    Value_Delta_List_Mean_4GPa = statistics.fmean(Value_Delta_List_4GPa)
    Rate_Per_ns_4GPa = Value_Delta_List_Mean_4GPa * 1000
    Dissociation_rates.append(Rate_Per_ns_4GPa)

    Value_Delta_List_5GPa = [(fitted_functionFiveGPa[x + 1] - fitted_functionFiveGPa[x]) for x in Value if
                             fitted_functionFiveGPa[x] > cutoff]
    Value_Delta_List_Mean_5GPa = statistics.fmean(Value_Delta_List_5GPa)
    Rate_Per_ns_5GPa = Value_Delta_List_Mean_5GPa * 1000
    Dissociation_rates.append(Rate_Per_ns_5GPa)

    Dissociation_rates = [x * -1. for x in Dissociation_rates]

    return Dissociation_rates

def get_MATLABFIT_dissociation_rates(TimestepList, Coefficient, Cutoff):
    Time = []
    for x in TimestepList:
        Time.append(((x - 400000) / int(4 * 10 ** 6)))
    FittedCurveValues = []
    for x in Time:
        Value = 48 * np.exp(Coefficient*(x))
        FittedCurveValues.append(Value)

    NumberList = list(range(0, len(Time)))

    if Cutoff == None:
        NanosecondRate = Coefficient * -1
        LnRate = np.log(NanosecondRate)

    else:
        TimeCutoff = []
        [TimeCutoff.append(Time[x]) for x in NumberList if Time[x] <= Cutoff]
        fitted_function = []
        for x in TimeCutoff:
            fitted_function.append(48 * np.exp(Coefficient * x))

        ShortIntactMoleculeList = []
        [ShortIntactMoleculeList.append(fitted_function[x]) for x in NumberList[:len(TimeCutoff)]]
        UnextrapolatedRate = (np.log((ShortIntactMoleculeList[-1] / 48))) / TimeCutoff[-1]

        Extrapolation_Constant = 1001 / len(TimeCutoff)
        NanosecondRate = UnextrapolatedRate * -1 * Extrapolation_Constant
        LnRate = np.log(NanosecondRate)

    return NanosecondRate, LnRate

def plot_variation_in_shear_stress_constanttemp(temperature, Pressures):
    Friction_Coefficient_Dataframe_Unnamed = pd.read_csv('C:/Users/eeo21/Documents/PhD/TCPDecompositionExperiments/Completed/AlphaFe/{}/1GPa/'
                                'fc_ave.dump'.format(temperature), sep=' ')
    Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe_Unnamed.rename(columns={'v_s_bot' : 'Shear Stress 1GPa', 'v_p_bot' : 'Normal Stress 1GPa'})

    for P in Pressures:
        Dataframe = pd.read_csv('C:/Users/eeo21/Documents/PhD/TCPDecompositionExperiments/Completed/AlphaFe/{}/{}/'
                                'fc_ave.dump'.format(temperature, P), sep=' ')
        Big_DataframeP = Dataframe.rename(columns= {'TimeStep': 'Timestep {}'.format(P),
                                                        'v_s_bot': 'Shear Stress {}'.format(P),
                                                        'v_p_bot': 'Normal Stress {}'.format(P)})

        Friction_Coefficient_Dataframe = pd.concat([Friction_Coefficient_Dataframe, Big_DataframeP], axis =1)
        Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe.dropna()

    Timestep = Friction_Coefficient_Dataframe.TimeStep.tolist()

    Shear_Stress_1GPa = Friction_Coefficient_Dataframe['Shear Stress 1GPa'].tolist()
    Shear_Stress_2GPa = Friction_Coefficient_Dataframe['Shear Stress 2GPa'].tolist()
    Shear_Stress_3GPa = Friction_Coefficient_Dataframe['Shear Stress 3GPa'].tolist()
    Shear_Stress_4GPa = Friction_Coefficient_Dataframe['Shear Stress 4GPa'].tolist()
    Shear_Stress_5GPa = Friction_Coefficient_Dataframe['Shear Stress 5GPa'].tolist()

    Shear_Stress_1GPa = [x / 10000 for x in Shear_Stress_1GPa]
    Shear_Stress_2GPa = [x / 10000 for x in Shear_Stress_2GPa]
    Shear_Stress_3GPa = [x / 10000 for x in Shear_Stress_3GPa]
    Shear_Stress_4GPa = [x / 10000 for x in Shear_Stress_4GPa]
    Shear_Stress_5GPa = [x / 10000 for x in Shear_Stress_5GPa]

    Time = []
    for x in Timestep:
        Time.append(((x - 400000) / int(4 * 10 ** 6)))

    fig3, ax3 = plt.subplots()
    ax3.set_title('Variation in Shear Stress at {}'.format(temperature))
    ax3.set_xlabel('Time (ns)')
    ax3.set_ylabel('Shear Stress (GPa)')
    ax3.plot(Time, Shear_Stress_1GPa, label='1GPa')
    ax3.plot(Time, Shear_Stress_2GPa, label='2GPa')
    ax3.plot(Time, Shear_Stress_3GPa, label='3GPa')
    ax3.plot(Time, Shear_Stress_4GPa, label='4GPa')
    ax3.plot(Time, Shear_Stress_5GPa, label='5GPa')
    ax3.legend()
    plt.show()

def plot_shear_stress_vs_normal_stress_different_sliding_speeds(Average_Shear_Stress_List_1, Average_Shear_Stress_List_2,
                                       Average_Shear_Stress_List_3, Average_Shear_Stress_List_4,
                                       Average_Shear_Stress_List_5, Average_Shear_Stress_List_6, Speed1, Speed2, Speed3, Speed4, Speed5, Speed6):
    x = np.array([1, 2, 3, 4, 5])
    a, b = np.polyfit(x, Average_Shear_Stress_List_1, 1)
    c, d = np.polyfit(x, Average_Shear_Stress_List_2, 1)
    e, f = np.polyfit(x, Average_Shear_Stress_List_3, 1)
    g, h = np.polyfit(x, Average_Shear_Stress_List_4, 1)
    i, j = np.polyfit(x, Average_Shear_Stress_List_5, 1)
    k, l = np.polyfit(x, Average_Shear_Stress_List_6, 1)

    fig1, ax2 = plt.subplots()
    ax2.set_title('Shear Stress vs Normal Stress at Different Sliding Speeds')
    ax2.set_xlabel('Normal Stress (GPa)')
    ax2.set_ylabel('Shear Stress (GPa)')
    ax2.scatter(x, Average_Shear_Stress_List_1)
    ax2.scatter(x, Average_Shear_Stress_List_2)
    ax2.scatter(x, Average_Shear_Stress_List_3)
    ax2.scatter(x, Average_Shear_Stress_List_4)
    ax2.scatter(x, Average_Shear_Stress_List_5)
    ax2.scatter(x, Average_Shear_Stress_List_6)
    ax2.plot(x, a * x + b, label=Speed1)
    ax2.plot(x, c * x + d, label=Speed2)
    ax2.plot(x, e * x + f, label=Speed3)
    ax2.plot(x, g * x + h, label=Speed4)
    ax2.plot(x, i * x + j, label=Speed5)
    ax2.plot(x, k * x + l, label=Speed6)
    ax2.legend()
    plt.show()

def plot_variation_in_shear_stress_constant_speed(speed, Pressures):
    Friction_Coefficient_Dataframe_Unnamed = pd.read_csv(
        'D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{}/400K/1GPa/'
        'fc_ave.dump'.format(speed), sep=' ')
    Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe_Unnamed.rename(
        columns={'v_s_bot': 'Shear Stress 1GPa', 'v_p_bot': 'Normal Stress 1GPa'})

    for P in Pressures:
        Dataframe = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{}/400K/{}/'
                                'fc_ave.dump'.format(speed, P), sep=' ')
        Big_DataframeP = Dataframe.rename(columns={'TimeStep': 'Timestep {}'.format(P),
                                                   'v_s_bot': 'Shear Stress {}'.format(P),
                                                   'v_p_bot': 'Normal Stress {}'.format(P)})

        Friction_Coefficient_Dataframe = pd.concat([Friction_Coefficient_Dataframe, Big_DataframeP], axis=1)
        Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe.dropna()


    Timestep = Friction_Coefficient_Dataframe.TimeStep.tolist()

    Shear_Stress_1GPa = Friction_Coefficient_Dataframe['Shear Stress 1GPa'].tolist()
    Shear_Stress_2GPa = Friction_Coefficient_Dataframe['Shear Stress 2GPa'].tolist()
    Shear_Stress_3GPa = Friction_Coefficient_Dataframe['Shear Stress 3GPa'].tolist()
    Shear_Stress_4GPa = Friction_Coefficient_Dataframe['Shear Stress 4GPa'].tolist()
    Shear_Stress_5GPa = Friction_Coefficient_Dataframe['Shear Stress 5GPa'].tolist()

    Shear_Stress_1GPa = [x / 10000 for x in Shear_Stress_1GPa]
    Shear_Stress_2GPa = [x / 10000 for x in Shear_Stress_2GPa]
    Shear_Stress_3GPa = [x / 10000 for x in Shear_Stress_3GPa]
    Shear_Stress_4GPa = [x / 10000 for x in Shear_Stress_4GPa]
    Shear_Stress_5GPa = [x / 10000 for x in Shear_Stress_5GPa]

    Time = []
    for x in Timestep:
        Time.append(((x - 400000) / int(4 * 10 ** 6)))

    fig3, ax3 = plt.subplots()
    ax3.set_title('Variation in Shear Stress at {}'.format(speed))
    ax3.set_xlabel('Time (ns)')
    ax3.set_ylabel('Shear Stress (GPa)')
    ax3.plot(Time, Shear_Stress_1GPa, label='1GPa')
    ax3.plot(Time, Shear_Stress_2GPa, label='2GPa')
    ax3.plot(Time, Shear_Stress_3GPa, label='3GPa')
    ax3.plot(Time, Shear_Stress_4GPa, label='4GPa')
    ax3.plot(Time, Shear_Stress_5GPa, label='5GPa')
    ax3.legend()
    plt.show()

def get_average_shear_normal_stress_and_average_mu_constant_speed(Speed, Pressures, EquilibriumFactor):
    Friction_Coefficient_Dataframe_Unnamed = pd.read_csv(
        'D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{}/400K/1GPa/'
        'fc_ave.dump'.format(Speed), sep=' ')
    Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe_Unnamed.rename(
        columns={'v_s_bot': 'Shear Stress 1GPa', 'v_p_bot': 'Normal Stress 1GPa'})

    for P in Pressures:
        Dataframe = pd.read_csv('D:/PhD/TCPDecompositionExperiments/Completed/DifferentSlidingSpeeds/{}/400K/{}/'
                                'fc_ave.dump'.format(Speed, P), sep=' ')
        Big_DataframeP = Dataframe.rename(columns={'Timestep': 'Timestep {}'.format(P),
                                                   'v_s_bot': 'Shear Stress {}'.format(P),
                                                   'v_p_bot': 'Normal Stress {}'.format(P)})

        Friction_Coefficient_Dataframe = pd.concat([Friction_Coefficient_Dataframe, Big_DataframeP], axis=1)
        Friction_Coefficient_Dataframe = Friction_Coefficient_Dataframe.dropna()

    # print(Friction_Coefficient_Dataframe)
    Mu_Final_Dataframe = Friction_Coefficient_Dataframe.iloc[:, [0, 1, 2, 4, 5, 7, 8, 10, 11, 13, 14]]
    Mu_Final_Dataframe = Mu_Final_Dataframe.iloc[EquilibriumFactor:, :]
    # print(Mu_Final_Dataframe)

    ShearStressMeans = Mu_Final_Dataframe[
        ['Shear Stress 1GPa', 'Shear Stress 2GPa', 'Shear Stress 3GPa', 'Shear Stress 4GPa',
         'Shear Stress 5GPa']].mean()
    Average_Shear_Stress_Dictionary = ShearStressMeans.to_dict()
    # print(ShearStressMeans)
    NormalStressMeans = Mu_Final_Dataframe[
        ['Normal Stress 1GPa', 'Normal Stress 2GPa', 'Normal Stress 3GPa', 'Normal Stress 4GPa',
         'Normal Stress 5GPa']].mean()
    NormalStressMeans = NormalStressMeans.to_dict()
    # print(NormalStressMeans)

    Average_Mu_Dictionary = {}

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
    # print(Average_Shear_Stress_List)
    Average_Mu_List = list(Average_Mu_Dictionary.values())
    Average_Shear_Stress_List = [x / 10000 for x in Average_Shear_Stress_List]  # Conversion to GPa
    # print(Average_Shear_Stress_List)

    return Average_Shear_Stress_List, Average_Mu_List, NormalStressMeans


from datafunctions import LQModel, LQModelHighLET
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

pd.set_option('mode.chained_assignment',None) #Ignoring errors from the pandas package 


#Setting Default font
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'cm'
mpl.rcParams['font.sans-serif'] = 'Times New Roman'
mpl.rcParams['font.family'] = 'Times New Roman'

#Declare colors
colors = ['maroon','orangered','orange','green','midnightblue']
colors4 = ['maroon','orange','green','midnightblue']
colors3 = ['maroon','orange','midnightblue']


def Plotting_one_surviving_curve(Dataframe, Alpha, Beta, Filepath, Title):
    #Getting information to set y and x limits for plot
    MaxDose = Dataframe['Dose [Gy]'].values.max()
    MinSurv = Dataframe['SurvWMulti'].values.min()

    #Calculating LQ-model
    TestingDose = np.linspace(0,int(MaxDose+2),int(MaxDose+3))
    LQFitted = LQModel(TestingDose,Alpha,Beta)
    #Plotting
    fig, ax = plt.subplots()
    plt.plot(Dataframe['Dose [Gy]'], Dataframe['SurvWMulti'],'o', markersize = 2, color = 'midnightblue', label='Average \nexperimental \ndata')
    plt.errorbar(Dataframe['Dose [Gy]'], Dataframe['SurvWMulti'], yerr=Dataframe['SurvWMulti SE'], fmt='none', elinewidth = 1, ecolor='midnightblue',capthick=1, capsize=2, label='')
    plt.plot(TestingDose, LQFitted,'--',color = 'midnightblue', linewidth = 1, label= 'LQ-model')
    
    #Formatting the plot
    plt.legend(loc='best', fontsize = 13)
    plt.grid(which='major',linestyle='--', linewidth='0.5', color='lightgray')
    plt.xlim(-0.05,MaxDose*1.05)
    plt.title(Title, fontsize = 15)
    plt.xlabel('Dose [Gy]', fontsize = 14)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.ylim(MinSurv*0.2,2)
    plt.ylabel('Surviving fraction', fontsize = 14)
    plt.yscale('log')
    plt.savefig('{}'.format(Filepath))
    plt.show()
    
def Plotting_proton(AverageSurv, Pos1A, Pos1B, Pos5A, Filepath, Title):
    #Getting information to set y and x limits for plot
    MaxDose = AverageSurv['Dose [Gy]'].values.max()
    MinSurv = AverageSurv['SurvWMulti'].values.min()
    Pos1Data = AverageSurv.loc[AverageSurv['Position']==1]
    Pos5Data = AverageSurv.loc[AverageSurv['Position']==5]
    C = AverageSurv.loc[AverageSurv['Dose [Gy]']==0]
    C['Corrected Dose [Gy]'] = 0
    Pos1 = pd.concat([C,Pos1Data])
    Pos5 = pd.concat([C,Pos5Data])

    TestingDosePos1 = np.linspace(0,int(MaxDose+1),int(MaxDose+2))
    TestingDosePos5 = np.linspace(0,15,16)

    LQPos1 = LQModel(TestingDosePos1,Pos1A,Pos1B)
    LQPos5 = LQModelHighLET(TestingDosePos5,Pos5A)

    fig, ax = plt.subplots()
    #Pos 1
    plt.plot(Pos1['Corrected Dose [Gy]'], Pos1['SurvWMulti'],'o', markersize = 3, color = 'midnightblue', label='Position 1')
    # plt.errorbar(Pos1['Corrected Dose [Gy]'], Pos1['SurvWMulti'], yerr=Pos1['SurvWMulti SE'], xerr = Pos1['Corrected Dose SE'],fmt='none', elinewidth = 1, ecolor='midnightblue',capthick=1, capsize=2, label='')
    plt.errorbar(Pos1['Corrected Dose [Gy]'], Pos1['SurvWMulti'], yerr=Pos1['SurvWMulti SE'],fmt='none', elinewidth = 1, ecolor='midnightblue',capthick=1, capsize=2, label='')
    plt.plot(TestingDosePos1, LQPos1, '--',color = 'midnightblue', linewidth = 1, label= 'LQ-model position 1')
    #Pos 5
    plt.plot(Pos5['Corrected Dose [Gy]'], Pos5['SurvWMulti'],'o', markersize = 3, color = 'maroon', label='Position 5')
    # plt.errorbar(Pos5['Corrected Dose [Gy]'], Pos5['SurvWMulti'], yerr=Pos5['SurvWMulti SE'], xerr = Pos5['Corrected Dose SE'],fmt='none', elinewidth = 1, ecolor='maroon',capthick=1, capsize=2, label='')
    plt.errorbar(Pos5['Corrected Dose [Gy]'], Pos5['SurvWMulti'], yerr=Pos5['SurvWMulti SE'],fmt='none', elinewidth = 1, ecolor='maroon',capthick=1, capsize=2, label='')
    plt.plot(TestingDosePos5, LQPos5, '--',color = 'maroon', linewidth = 1, label= 'LQ-model position 5')
    #Formatting the plot
    plt.legend(loc='best',fontsize=14)
    plt.grid(which='major',linestyle='--', linewidth='0.5', color='lightgray')
    plt.xlim(-0.05,MaxDose*1.08)
    plt.title(Title, fontsize=16)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.xlabel('Dose [Gy]',fontsize=14)
    plt.ylim(MinSurv*0.1,1.5)
    plt.ylabel('Surviving fraction',fontsize=14)
    plt.yscale('log')
    plt.savefig('%s'  % (Filepath ))
    # plt.show()

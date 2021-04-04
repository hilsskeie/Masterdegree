from statistics import stdev

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

pd.set_option('mode.chained_assignment',None) #Ignoring errors from the pandas package 




def Average_signal(Dataframe, irradiation):
    if 'xray' in irradiation:
        number_of_replicates = Dataframe.groupby(['Time', 'Dose']).size()
        average_signal = Dataframe.groupby(['Time', 'Dose']).mean()
        std_average_signal = Dataframe.groupby(['Time', 'Dose']).std()
    elif 'proton' in irradiation:
        number_of_replicates = Dataframe.groupby(['Position','Time', 'Dose']).size()
        average_signal = Dataframe.groupby(['Position','Time', 'Dose']).mean()
        std_average_signal = Dataframe.groupby(['Position','Time', 'Dose']).std()


    #Create a new column with error squared
    Dataframe['Ratio error full^2']=Dataframe['Ratio error full']**2
    Dataframe['Ratio error G1^2'] = Dataframe['Ratio error G1']**2
    Dataframe['Ratio error G2^2'] = Dataframe['Ratio error G2']**2
    Dataframe['Ratio error S^2'] = Dataframe['Ratio error S']**2
    Dataframe['Ratio error G1G2^2'] = Dataframe['Ratio error G1G2']**2

    #Sum error for each time and dose point
    if 'xray' in irradiation:
        sum_dataframe = Dataframe.groupby(['Time', 'Dose']).sum()
    elif 'proton' in irradiation:
        sum_dataframe = Dataframe.groupby(['Position','Time', 'Dose']).sum()
    
    #Calculate the total error from each of the median errors
    se_median_full = np.sqrt(sum_dataframe['Ratio error full^2']/number_of_replicates**2)
    se_median_G1G2 = np.sqrt(sum_dataframe['Ratio error G1G2^2']/number_of_replicates**2)
    se_median_G1 = np.sqrt(sum_dataframe['Ratio error G1^2']/number_of_replicates**2)
    se_median_G2 = np.sqrt(sum_dataframe['Ratio error G2^2']/number_of_replicates**2)
    se_median_S = np.sqrt(sum_dataframe['Ratio error S^2']/number_of_replicates**2)

    #Calculate error for the mean values
    se_average_ratio_full = std_average_signal['Ratio full']/np.sqrt(number_of_replicates) 
    se_average_ratio_G1G2 = std_average_signal['Ratio G1G2']/np.sqrt(number_of_replicates)
    se_average_ratio_G1 = std_average_signal['Ratio G1']/np.sqrt(number_of_replicates)
    se_average_ratio_G2 = std_average_signal['Ratio G2']/np.sqrt(number_of_replicates)
    se_average_ratio_S = std_average_signal['Ratio S']/np.sqrt(number_of_replicates)

    average_signal['Mean error full'] = np.sqrt(se_median_full**2 + se_average_ratio_full**2)
    average_signal['Mean error G1G2'] = np.sqrt(se_median_G1G2**2 + se_average_ratio_G1G2**2)
    average_signal['Mean error G1'] = np.sqrt(se_median_G1**2 + se_average_ratio_G1**2)
    average_signal['Mean error G2'] = np.sqrt(se_median_G2**2 + se_average_ratio_G2**2)
    average_signal['Mean error S'] = np.sqrt(se_median_S**2 + se_average_ratio_S**2)
    
    #Get index where there is only one replicate
    One_replicate_index = number_of_replicates.index[number_of_replicates <2].tolist()
    
    if len(One_replicate_index)>0:
        #Getting the values for Time and Dose from the index
        Time = One_replicate_index[0][0]
        Dose = One_replicate_index[0][1]

        #Getting the ratio error from the original dataframe for the given replicate
        Ratio_error_full = Dataframe.loc[(Dataframe.Time == Time) & (Dataframe.Dose == Dose),'Ratio error full']
        Ratio_error_G1G2 = Dataframe.loc[(Dataframe.Time == Time) & (Dataframe.Dose == Dose),'Ratio error G1G2']     
        Ratio_error_G1 = Dataframe.loc[(Dataframe.Time == Time) & (Dataframe.Dose == Dose),'Ratio error G1']     
        Ratio_error_G2 = Dataframe.loc[(Dataframe.Time == Time) & (Dataframe.Dose == Dose),'Ratio error G2']     
        Ratio_error_S = Dataframe.loc[(Dataframe.Time == Time) & (Dataframe.Dose == Dose),'Ratio error S']     
        
        #Setting the Ratio error as Mean error when only one replicate
        average_signal.loc[One_replicate_index, 'Mean error full'] = Ratio_error_full.values[0]
        average_signal.loc[One_replicate_index, 'Mean error G1G2'] = Ratio_error_G1G2.values[0]
        average_signal.loc[One_replicate_index, 'Mean error G1'] = Ratio_error_G1.values[0]
        average_signal.loc[One_replicate_index, 'Mean error G2'] = Ratio_error_G2.values[0]
        average_signal.loc[One_replicate_index, 'Mean error S'] = Ratio_error_S.values[0]
    else:
        pass
    return average_signal


def Normalize_to_Controls(filepath):
    filepathinformation = filepath.split('/') #Do not want the folders before file
    fileinformation = filepathinformation[-1].split('_')
    date = fileinformation[0]
    cell_line = fileinformation[1]
    irradiation = fileinformation[2]
    if 'proton' in irradiation.lower():
        position = fileinformation[3]
        
    else:
        pass

    #read the data that we want to plot
    data = pd.read_csv(filepath, index_col = [0]).drop(columns = ['Replicate','Fractions'])
    samples = data.loc[data.Dose!=0]
    control = data.loc[data.Dose==0]
        
    #Calculate average values based on irradiation, groupby position as well if proton
    average_samples = Average_signal(samples, irradiation)
    average_control = Average_signal(control,irradiation)
    
    #Get length of indexes to check how to normalize to control
    SampleIndex = average_samples.index.get_level_values(0).unique()
    ControlIndex = average_control.index.get_level_values(0).unique()
    
    if len(SampleIndex.values) != len(ControlIndex.values) and len(ControlIndex.values)==1:
        #Do not have a set of controls for each timepoint
        control_corr_full = average_samples['Ratio full']/average_control['Ratio full'].values[0]
        control_corr_full_se = control_corr_full*np.sqrt((average_samples['Mean error full']/average_samples['Ratio full'])**2 + (average_control['Mean error full'].values[0]/average_control['Ratio full'].values[0])**2)
        control_corr_G1G2 = average_samples['Ratio G1G2']/average_control['Ratio G1G2'].values[0]
        control_corr_G1G2_se = control_corr_G1G2*np.sqrt((average_samples['Mean error G1G2']/average_samples['Ratio G1G2'])**2 + (average_control['Mean error G1G2'].values[0]/average_control['Ratio G1G2'].values[0])**2)
        control_corr_G1 = average_samples['Ratio G1']/average_control['Ratio G1'].values[0]
        control_corr_G1_se = control_corr_G1*np.sqrt((average_samples['Mean error G1']/average_samples['Ratio G1'])**2 + (average_control['Mean error G1'].values[0]/average_control['Ratio G1'].values[0])**2)
        control_corr_G2 = average_samples['Ratio G2']/average_control['Ratio G2'].values[0]
        control_corr_G2_se = control_corr_G2*np.sqrt((average_samples['Mean error G2']/average_samples['Ratio G2'])**2 + (average_control['Mean error G2'].values[0]/average_control['Ratio G2'].values[0])**2)
        control_corr_S = average_samples['Ratio S']/average_control['Ratio S'].values[0]
        control_corr_S_se = control_corr_S*np.sqrt((average_samples['Mean error S']/average_samples['Ratio S'])**2 + (average_control['Mean error S'].values[0]/average_control['Ratio S'].values[0])**2)

    elif len(SampleIndex.values) == len(ControlIndex.values):
        #Have a set of controls for each timepoint
        #Adds time as column and remove as index for control
        average_samples['Time'] = average_samples.index.get_level_values(0)
        average_control.index = average_control.index.droplevel(1)

        #Calculate control corrected values and error
        control_corr_full = average_samples['Ratio full']/average_samples['Time'].map(average_control['Ratio full'])
        control_corr_full_se = control_corr_full*np.sqrt((average_samples['Mean error full']/average_samples['Ratio full'])**2 + average_samples['Time'].map(average_control['Mean error full']/average_control['Ratio full'])**2)
        control_corr_G1G2 = average_samples['Ratio G1G2']/average_samples['Time'].map(average_control['Ratio G1G2'])
        control_corr_G1G2_se = control_corr_G1G2*np.sqrt((average_samples['Mean error G1G2']/average_samples['Ratio G1G2'])**2 + average_samples['Time'].map(average_control['Mean error G1G2']/average_control['Ratio G1G2'])**2)
        control_corr_G1 = average_samples['Ratio G1']/average_samples['Time'].map(average_control['Ratio G1'])
        control_corr_G1_se = control_corr_G1*np.sqrt((average_samples['Mean error G1']/average_samples['Ratio G1'])**2 + average_samples['Time'].map(average_control['Mean error G1']/average_control['Ratio G1'])**2)
        control_corr_G2 = average_samples['Ratio G2']/average_samples['Time'].map(average_control['Ratio G2'])
        control_corr_G2_se = control_corr_G2*np.sqrt((average_samples['Mean error G2']/average_samples['Ratio G2'])**2 + average_samples['Time'].map(average_control['Mean error G2']/average_control['Ratio G2'])**2)
        control_corr_S = average_samples['Ratio S']/average_samples['Time'].map(average_control['Ratio S'])
        control_corr_S_se = control_corr_S*np.sqrt((average_samples['Mean error S']/average_samples['Ratio S'])**2 + average_samples['Time'].map(average_control['Mean error S']/average_control['Ratio S'])**2)
    else:
        print('Did not normalize to control')
    
    average_ratio_control_corrected = pd.DataFrame({'Ratio full':control_corr_full, 'Ratio G1G2': control_corr_G1G2, 'Ratio G1': control_corr_G1, 'Ratio G2': control_corr_G2, 'Ratio S': control_corr_S, 'Error full': control_corr_full_se, 'Error G1G2': control_corr_G1G2_se, 'Error G1': control_corr_G1_se, 'Error S': control_corr_S_se, 'Error G2': control_corr_G2_se}, index = average_samples.index)
    
    sorted_average = average_ratio_control_corrected.sort_values(by= 'Dose', ascending = False)
    return sorted_average


def CellCycle(Dataframe,G1Values,G2Values):
    #Dividing the dataframe into two separate based on cell cycle
    G1min = Dataframe.loc[Dataframe['FL2-A']>=G1Values[0]]
    G1 = G1min.loc[G1min['FL2-A']<=G1Values[1]] 
    G2min = Dataframe.loc[Dataframe['FL2-A']>=G2Values[0]]
    G2 = G2min.loc[G2min['FL2-A']<=G2Values[1]]
    S = Dataframe.loc[(Dataframe['FL2-A'] > G1Values[1]) & (Dataframe['FL2-A'] < G2Values[0])]
    return G1, G2, S

def CellCycleFiltering(Dataframe,CellCyclelimits):
    #Filtering out the cells outside the cell cycle
    Min = Dataframe.loc[Dataframe['FL2-A']>=CellCyclelimits[0][0]]
    Max = Min.loc[Min['FL2-A']<=CellCyclelimits[1][0]]
    return Max

def DNAContentCorrection(Dataframe,BCG1Median):
    #Correct each H2AX siganl by the PI signal (FL1-A/FL2-A) and finding min,max,median and stanard error
    H2AXSignal = Dataframe.loc[:,'FL1-A']
    CellContent = Dataframe.loc[:,'FL2-A']/BCG1Median

    Dataframe.loc[:,'H2AX DNA Corr'] = H2AXSignal/CellContent
    H2AXmin = Dataframe['H2AX DNA Corr'].min()
    H2AXmax = Dataframe['H2AX DNA Corr'].max()
    H2AXmedian = Dataframe['H2AX DNA Corr'].median()
    H2AXSE = stdev(Dataframe['H2AX DNA Corr'])/np.sqrt(len(Dataframe))
    SortedDataframe = Dataframe.sort_values(by='H2AX DNA Corr')
    NewDF = SortedDataframe.reset_index()
    return NewDF,H2AXmedian, H2AXSE

def NumberOfCells(Dataframe,Limits):
    #Finding the amount of elements that is within the provided limits
    Dataframemin = Dataframe.loc[Dataframe['FL2-A']>=Limits[0]]
    NewDf = Dataframemin.loc[Dataframe['FL2-A']<=Limits[1]]
    Amountwithinlimits = len(NewDf)/len(Dataframe)*100
    return NewDf, Amountwithinlimits

def Phaseanalyze(Dataframe, Type, Irradiation):
    #Finding the average amount and calculating the standard error
    if Irradiation == 'xray':
        Average_dist = Dataframe.groupby(['Time', 'Dose']).mean()
        Number_of_replicates = Dataframe.groupby(['Time', 'Dose']).size()
        Std_dist = Dataframe.groupby(['Time', 'Dose']).std()
    elif Irradiation =='proton':
        Average_dist = Dataframe.groupby(['Position','Time', 'Dose']).mean()
        Number_of_replicates = Dataframe.groupby(['Position','Time', 'Dose']).size()
        Std_dist = Dataframe.groupby(['Position','Time', 'Dose']).std()
    
    #Calculating the standard error
    SE_dist_G1 = Std_dist[Type+' G1 amount']/np.sqrt(Number_of_replicates)
    SE_dist_G2 = Std_dist[Type+' G2 amount']/np.sqrt(Number_of_replicates)
    SE_dist_S = Std_dist[Type+' S amount']/np.sqrt(Number_of_replicates)
    #Adding error to dataframe
    G1 = pd.concat([Average_dist[Type+' G1 amount'],SE_dist_G1], axis = 1).rename(columns ={Type+' G1 amount':'Amount',0:'Error'})
    G2 = pd.concat([Average_dist[Type+' G2 amount'],SE_dist_G2], axis = 1).rename(columns ={Type+' G2 amount':'Amount',0:'Error'})
    S = pd.concat([Average_dist[Type+' S amount'],SE_dist_S], axis = 1).rename(columns ={Type+' S amount':'Amount',0:'Error'})
    
    return G1,S,G2



def H2AX_Find_Median(Sample,SampleLimits,BC,BCLimits):
    #Getting the signal in G1 and G2
    BCG1, BCG2, BCS = CellCycle(BC,BCLimits[0],BCLimits[1])
    SampleG1, SampleG2, SampleS = CellCycle(Sample,SampleLimits[0],SampleLimits[1])

    BCG1Median = BCG1['FL2-A'].median()

    #Joining G1 and G2 data
    BCG1G2 = pd.DataFrame.append(BCG1,BCG2, ignore_index=True)
    SampleG1G2 = pd.DataFrame.append(SampleG1,SampleG2, ignore_index=True)
    
    """
    Correct for FL2-A (DNA Content)
    """
    #Whole cell cycle
    SampleDNACorrected,SampleDNACorrMedian, SampleDNACorrSE = DNAContentCorrection(Sample, BCG1Median)
    BCDNACorrected, BCDNACorrMedian, BCDNACorrSE = DNAContentCorrection(BC,BCG1Median)
    #G1 
    SampleG1DNACorrected,SampleG1DNACorrMedian, SampleG1DNACorrSE = DNAContentCorrection(SampleG1, BCG1Median)
    BCG1DNACorrected, BCG1DNACorrMedian, BCG1DNACorrSE = DNAContentCorrection(BCG1,BCG1Median)
    #G2
    SampleG2DNACorrected,SampleG2DNACorrMedian, SampleG2DNACorrSE = DNAContentCorrection(SampleG2, BCG1Median)
    BCG2DNACorrected, BCG2DNACorrMedian, BCG2DNACorrSE = DNAContentCorrection(BCG2,BCG1Median)
    #S
    SampleSDNACorrected,SampleSDNACorrMedian, SampleSDNACorrSE = DNAContentCorrection(SampleS, BCG1Median)
    BCSDNACorrected, BCSDNACorrMedian, BCSDNACorrSE = DNAContentCorrection(BCS,BCG1Median)
    
    #G1 and G2
    SampleG1G2DNACorrected, SampleG1G2DNACorrMedian,SampleG1G2DNACorrSE = DNAContentCorrection(SampleG1G2, BCG1Median)
    BCG1G2DNACorrected, BCG1G2DNACorrMedian, BCG1G2DNACorrSE = DNAContentCorrection(BCG1G2,BCG1Median)


    #Caluclating median ratios with errors
    Full_ratio = SampleDNACorrMedian/BCDNACorrMedian
    Full_se = Full_ratio*np.sqrt((SampleDNACorrSE/SampleDNACorrMedian)**2 + (BCDNACorrSE/BCDNACorrMedian)**2)

    G1_ratio = SampleG1DNACorrMedian/BCG1DNACorrMedian
    G1_se = G1_ratio*np.sqrt((SampleG1DNACorrSE/SampleG1DNACorrMedian)**2 + (BCG1DNACorrSE/BCG1DNACorrMedian)**2)
    G2_ratio = SampleG2DNACorrMedian/BCG2DNACorrMedian
    G2_se = G2_ratio*np.sqrt((SampleG2DNACorrSE/SampleG2DNACorrMedian)**2 + (BCG2DNACorrSE/BCG2DNACorrMedian)**2)
    S_ratio = SampleSDNACorrMedian/BCSDNACorrMedian
    S_se = S_ratio*np.sqrt((SampleSDNACorrSE/SampleSDNACorrMedian)**2 + (BCSDNACorrSE/BCSDNACorrMedian)**2)

    G1G2_ratio = SampleG1G2DNACorrMedian/BCG1G2DNACorrMedian
    G1G2_se = G1G2_ratio*np.sqrt((SampleG1G2DNACorrSE/SampleG1G2DNACorrMedian)**2 + (BCG1G2DNACorrSE/BCG1G2DNACorrMedian)**2)

    Sample_Length = len(Sample)
    Sample_Length_G1 = len(SampleG1)/Sample_Length*100
    Sample_Length_G2 = len(SampleG2)/Sample_Length*100
    Sample_Length_S = len(SampleS)/Sample_Length*100
    BC_Length = len(BC)
    BC_Length_G1 = len(BCG1)/BC_Length*100
    BC_Length_G2 = len(BCG2)/BC_Length*100
    BC_Length_S = len(BCS)/BC_Length*100


    Median_Data = pd.DataFrame({'Sample median full' : SampleDNACorrMedian, 'Sample error full' : SampleDNACorrSE, \
        'BC median full' : BCDNACorrMedian, 'BC error full': BCDNACorrSE,\
        'Ratio full': Full_ratio,'Ratio error full': Full_se,\
        'Sample median G1' : SampleG1DNACorrMedian, 'Sample error G1': SampleG1DNACorrSE, \
        'Sample median G2' : SampleG2DNACorrMedian, 'Sample error G2': SampleG2DNACorrSE,\
        'Sample median S' : SampleSDNACorrMedian, 'Sample error S': SampleSDNACorrSE,\
        'BC median G1' : BCG1DNACorrMedian, 'BC error G1': BCG1DNACorrSE, \
        'BC median G2' : BCG2DNACorrMedian, 'BC error G2': BCG2DNACorrSE,\
        'BC median S' : BCSDNACorrMedian, 'BC error S': BCSDNACorrSE,\
        'Ratio G1': G1_ratio, 'Ratio error G1': G1_se, \
        'Ratio G2': G2_ratio, 'Ratio error G2': G2_se,\
        'Ratio S': S_ratio, 'Ratio error S': S_se,\
        'Sample G1 amount': Sample_Length_G1, 'Sample G2 amount': Sample_Length_G2, 'Sample S amount': Sample_Length_S,
        'BC G1 amount': BC_Length_G1,'BC G2 amount': BC_Length_G2,'BC S amount': BC_Length_S,
        'Ratio G1G2': G1G2_ratio, 'Ratio error G1G2': G2_se}, index = [0])

    return Median_Data
"""
Possible input:
Filepath of specific file - Results in a file and plots of the results for the given experiment
Filepath of Results folder - Results in the average result of all experiments of the given irradiation type
Dose correction factor - As option number 3, if film has been analyzed pass the correction factor. If not given it will be set to 1.0
"""
from readfile import ReadCount,ReadSurv, ReadFolder
from datafunctions import *
from plotting import *

import pandas as pd
import numpy as np
import sys
import os

Input = sys.argv[1]
try:
    DoseCorrectionFactor = float(sys.argv[2])
except IndexError:
    DoseCorrectionFactor = 1.0

try:
    if 'Results' in Input:
        #Get filenames
        AverageSurvFilenames = ReadFolder(Input+'Average/')
        AllSurvFilenames = ReadFolder(Input+'All/')
        AverageSurv,AllSurv = ReadSurv(AverageSurvFilenames,AllSurvFilenames,Input)
        if 'xray' in Input:
            #Calculate average for plotting
            Average = AverageResults(AverageSurv,'xray')
            #Want the dose as a column
            DatatoPlot = Average.reset_index()
            #Fit to LQ-model without weights
            Alpha,Beta,RSquared = LQModelFit(AverageSurv['Dose [Gy]'],AverageSurv['SurvWMulti'],AverageSurv['SurvWMulti SE'],AllSurv['Dose [Gy]'],AllSurv['SurvWMulti'],AllSurv['SurvWMulti SE'],'no')
            #Fit to LQ-model with weights
            # Alpha,Beta,RSquared = LQModelFit(AverageSurv['Dose [Gy]'],AverageSurv['SurvWMulti'],AverageSurv['SurvWMulti SE'],AllSurv['Dose [Gy]'],AllSurv['SurvWMulti'],AllSurv['SurvWMulti SE'],'yes')
            print('Alpha: {:.3f} +\- {:.3f}'.format(Alpha[0],Alpha[1]))
            print('Beta: {:.4f} +\- {:.4f}'.format(Beta[0], Beta[1]))
            print('R^2: {:.4f}'.format(RSquared))

            # #Plotting
            Plotting_one_surviving_curve(DatatoPlot,Alpha[0],Beta[0],Input+'Average_surv.pdf', 'Colony formation results for A549 cells after X-ray irradiation')
            # Plotting_one_surviving_curve(AllSurv,Alpha[0],Beta[0],Input+'All_surv.pdf', 'Colony formation results for A549 cells after X-ray irradiation')

        elif 'proton' in Input:
            #Calculating average for plotting
            Average = AverageResults(AverageSurv,'proton')
            #Want the dose as a column
            DatatoPlot = Average.reset_index()

            #Dividing into positions for LQ-fit
            Pos1 = AverageSurv.loc[AverageSurv['Position']==1] 
            Pos5 = AverageSurv.loc[AverageSurv['Position']==5]
            #Getting the data to test the fit 
            AllPos1 = AllSurv.loc[AllSurv['Position']==1]
            AllPos5Fit = AllSurv.loc[(AllSurv['Position']==5)&(AllSurv['Dose [Gy]']<6)]

            #Excluding the SF-plateau
            FitPos5 = Pos5.loc[Pos5['Dose [Gy]']<6]

            #Fitting to models
            AlphaPos1,BetaPos1,RSquaredPos1 = LQModelFit(Pos1['Corrected Dose [Gy]'],Pos1['SurvWMulti'],Pos1['SurvWMulti SE'],AllPos1['Corrected Dose [Gy]'],AllPos1['SurvWMulti'],AllPos1['SurvWMulti SE'],'no')
            AlphaPos5, RSquaredPos5 = LQModelFitHighLET(FitPos5['Corrected Dose [Gy]'],FitPos5['SurvWMulti'],FitPos5['SurvWMulti SE'],AllPos5Fit['Corrected Dose [Gy]'],AllPos5Fit['SurvWMulti'],AllPos5Fit['SurvWMulti SE'],'no')
            
            #With Weights
            # AlphaPos1,BetaPos1,RSquaredPos1 = LQModelFit(Pos1['Corrected Dose [Gy]'],Pos1['SurvWMulti'],Pos1['SurvWMulti SE'],AllPos1['Corrected Dose [Gy]'],AllPos1['SurvWMulti'],AllPos1['SurvWMulti SE'],'yes')
            # AlphaPos5, RSquaredPos5 = LQModelFitHighLET(FitPos5['Corrected Dose [Gy]'],FitPos5['SurvWMulti'],FitPos5['SurvWMulti SE'],AllPos1['Corrected Dose [Gy]'],AllPos1['SurvWMulti'],AllPos1['SurvWMulti SE'],'yes')

            Alpha = [AlphaPos1,AlphaPos5]
            print('Alpha pos 1: {:.2f} +\- {:.2f}'.format(Alpha[0][0],Alpha[0][1]))
            print('Beta pos 1: {:.3f} +\- {:.3f}'.format(BetaPos1[0], BetaPos1[1]))
            print('R^2 pos 1: {:.3f} \n'.format(RSquaredPos1))
            print('Alpha pos 5: {:.2f} +\- {:.2f}'.format(Alpha[1][0],Alpha[1][1]))
            print('R^2 pos5 : {:.3f} \n'.format(RSquaredPos5))

            #Setting dose as column not index
            DosetoColumn = Average.reset_index()
     
            #Plotting
            Plotting_proton(DosetoColumn,Alpha[0][0],BetaPos1[0],Alpha[1][0],Input+'Average_surv.pdf','Colony formation results for A549 cells after proton irradiation')
            # Plotting_proton(AllSurv,Alpha[0][0],BetaPos1[0],Alpha[1][0],Input+'All_surv.pdf','Colony formation results for A549 cells after proton irradiation')
        else:
            print('Something went wrong')

        #Saving to file
        Average.to_csv(Input+'Average_results.csv')
        AllSurv.to_csv(Input+'All_results.csv')

    else:
        Filepathinfo = Input.split('/')
        Resultfilepath = '/'.join(Filepathinfo[:-1])+'/Results/'

        #Create folder if it does not exists
        try:
            CreateFolderResults = os.mkdir(Resultfilepath)
            CreateFolderAverage = os.mkdir(Resultfilepath+'Average/')
            CreateFolderAll = os.mkdir(Resultfilepath+'All/')
        except FileExistsError:
            pass
        Filename = Filepathinfo[-1]

        #Reading the data
        Samples,Control,Multiplicity = ReadCount(Input)
        
        if 'proton' in Input:
            #Correct for dose given 
            DoseCorrectedSamples = DoseCorrection(Samples,DoseCorrectionFactor)
            #Correct for the counting ring
            SamplesCorrectedCellNumbers = AreaCorrection(DoseCorrectedSamples,5.2)
            ControlCorrectedCellNumbers = AreaCorrection(Control,5.2)
            
            #Calculate Average control and PE
            MeanControl,PlatingEfficiency = AverageCountControl(ControlCorrectedCellNumbers)
            MeanControl['Corrected Dose [Gy]'] = 0.0
            
            MeanSample = SamplesCorrectedCellNumbers.groupby(['Position','Dose [Gy]']).mean()
            StdCountSample = SamplesCorrectedCellNumbers.groupby(['Position','Dose [Gy]']).std()
            LengthSampels = SamplesCorrectedCellNumbers.groupby(['Position','Dose [Gy]']).size()
            MeanSample['SE Count'] = StdCountSample['Count']/np.sqrt(LengthSampels)
            MeanSample['Corrected Dose SE'] = StdCountSample['Corrected Dose [Gy]']/np.sqrt(LengthSampels)
            
            #Calculate the Survivingfraction
            AverageSurvData = AverageSurvivalWMulti(MeanSample,MeanControl,Multiplicity,PlatingEfficiency)
            AllSurvData = SurvWMultiAllSamples(SamplesCorrectedCellNumbers,MeanControl,Multiplicity['Count'].values[0])
            #Adding Control values to all data
            AllSurvData.loc[-1] = MeanControl.drop(['Cells/dish','Count'])
            if len(LengthSampels.loc[LengthSampels == 1])>0:
                #If only one sample the SE is taken from the given sample, not average
                Position = LengthSampels.loc[LengthSampels== 1].index.get_level_values(0)
                Dose = LengthSampels.loc[LengthSampels== 1].index.get_level_values(1)
                Data = AllSurvData.loc[(AllSurvData['Position']==Position[0]) & (AllSurvData['Dose [Gy]']==Dose[0])]

                AverageSurvData.loc[(AverageSurvData['Position']==Position[0]) & (AverageSurvData['Dose [Gy]']==Dose[0]),'Surv SE'] = Data['Surv SE'].values[0]
                AverageSurvData.loc[(AverageSurvData['Position']==Position[0]) & (AverageSurvData['Dose [Gy]']==Dose[0]),'SurvWMulti SE'] = Data['SurvWMulti SE'].values[0]
                
            #Finding the alpha- and beta values for pos 1 for this experiment based on the results and the LQ model
            FittingDatapos1 = AllSurvData.loc[AllSurvData['Surv']>0]
            LQFittedparametersPos1, LQFittedpcovPos1 = optimize.curve_fit(LogLQModel,FittingDatapos1['Dose [Gy]'].loc[FittingDatapos1['Position']==1],np.log(FittingDatapos1['SurvWMulti'].loc[FittingDatapos1['Position']==1]))
            AlphaPos1 = LQFittedparametersPos1[0]
            BetaPos1 = LQFittedparametersPos1[1]
            ParametersErrorPos1 = np.sqrt(np.diag(LQFittedpcovPos1))
        
            #Finding the alpha-value for pos 5 for this experiment based on the results and the LQ model, filtering on doses up to 7 Gy due to the tail
            FittingDatapos5 = AllSurvData.loc[(AllSurvData['Dose [Gy]']<6) & (AllSurvData['Position']==5) & (AllSurvData['Surv']>0)]
            LQFittedparametersPos5, LQFittedpcovPos5 = optimize.curve_fit(LogLQModelHighLET,FittingDatapos5['Dose [Gy]'],np.log(FittingDatapos5['SurvWMulti']), method = 'trf')
            AlphaPos5 = LQFittedparametersPos5[0]
            BetaPos5 = 0
            ParametersErrorPos5 = np.sqrt(LQFittedpcovPos5[0][0])

            #Plot
            Plotfilepath = Resultfilepath+Filename[:-5]+'.pdf'
            Plotting_proton(AverageSurvData,AlphaPos1,BetaPos1,AlphaPos5,Plotfilepath,'Surviving fraction for A549 cells after proton irradiation')

            #Remove unnecessary columns
            AverageSurvData = AverageSurvData.drop(columns=['Factor','MU prescribed','MU delivered'])
            AllSurvData = AllSurvData.drop(columns=['Cells/dish','Count','Factor','MU prescribed','MU delivered'])


        elif 'xray' in Input:
            #Find Average control counts and PE
            MeanControl,PlatingEfficiency = AverageCountControl(Control)

            #Calculating average and standard error
            MeanSample = Samples.groupby(['Dose [Gy]']).mean()
            StdCountSample = Samples.groupby(['Dose [Gy]']).std()
            LengthSampels = Samples.groupby(['Dose [Gy]']).size()
            MeanSample['SE Count'] = StdCountSample['Count']/np.sqrt(LengthSampels)
            
            #Calculate the data
            AverageSurvData = AverageSurvivalWMulti(MeanSample,MeanControl,Multiplicity,PlatingEfficiency)
            AllSurvData = SurvWMultiAllSamples(Samples,MeanControl,Multiplicity['Count'].values[0]) 
            #Adding Control values to all data
            AllSurvData.loc[-1] = MeanControl.drop(['Cells/dish','Count'])
            
            #Finding the alpha- and beta values for this experiment based on the results and the LQ model
            Alpha,Beta,RSquared = LQModelFit(AllSurvData['Dose [Gy]'],AllSurvData['SurvWMulti'],AllSurvData['SurvWMulti SE'],AllSurvData['Dose [Gy]'],AllSurvData['SurvWMulti'],AllSurvData['SurvWMulti SE'], 'no')
            print('R^2: {:.4f}'.format(RSquared))

            #Plot
            PlotFilePathAndName = Resultfilepath+Filename[:-5] + '.pdf'
            Plotting_one_surviving_curve(AverageSurvData, Alpha[0], Beta[0], PlotFilePathAndName, 'Surviving fraction for A549 cells after x-ray irradiation')

        else:
            print('Did not find the correct way to analyze the colony data')

        #Write to file
        AverageSurvData.to_csv('%sAverage/average_%s_surv_pe_%.3f.csv' % (Resultfilepath,Filename[:-5],PlatingEfficiency))
        AllSurvData.sort_values(by=['Dose [Gy]']).to_csv('%sAll/all_%s_surv_pe_%.3f.csv' % (Resultfilepath,Filename[:-5],PlatingEfficiency))
        
    
except FileNotFoundError:
    print('Did not find any file, maybe missing a \'/\'?')
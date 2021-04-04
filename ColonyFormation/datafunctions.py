import numpy as np
import pandas as pd
import sklearn.metrics as sk
from statistics import stdev
from scipy import optimize


   
def AreaCorrection(Samples,DiameterCellDish):
   #Samples is a dataframe, function to correct for the fact that we don't count the entire dish, so we can compare to x-ray data.
   CellNumber = Samples['Cells/dish']
   DiameterCountingRing = 0.4
   AreaDish = np.pi*(DiameterCellDish/2)**2
   AreaRing = np.pi*(DiameterCountingRing/2)**2
   CorrectionFactor = AreaRing/AreaDish
   NewCellNumber = CellNumber - CellNumber*CorrectionFactor
   Samples['Cells/dish'] = NewCellNumber
   return Samples

def AverageCountControl(Control):
   #Calculating mean and error values for control
   MeanControl = Control.mean()
   PlatingEfficiency = MeanControl['Count']/MeanControl['Cells/dish']
   MeanControl['SE Count'] = Control['Count'].std()/np.sqrt(len(Control))
   #Setting Surv to 1.0, (eqaual to doing the calculations)
   MeanControl['Surv'] = 1.0
   MeanControl['Surv SE'] = np.sqrt((MeanControl['Count']*(MeanControl['SE Count']/MeanControl['Cells/dish'])/(MeanControl['Cells/dish']*PlatingEfficiency**2))**2+(MeanControl['SE Count']/(PlatingEfficiency*MeanControl['Cells/dish']))**2)
   return MeanControl, PlatingEfficiency

def AverageResults(AverageSurv,Type):
   #Square the average data to later calculate the standard error
   SquaredData = AverageSurv.pow(2)
   if Type =='xray':
      Average = AverageSurv.groupby(['Dose [Gy]']).mean()
      STD = AverageSurv.groupby(['Dose [Gy]']).std()
      Size = AverageSurv.groupby(['Dose [Gy]']).size()
      SE = STD['SurvWMulti']/np.sqrt(Size)
      #Calculate sum 
      Sum = SquaredData.groupby(['Dose [Gy]']).sum()
      Sum.index = Size.index
      SumSE = np.sqrt((Sum['SurvWMulti SE']/Size**2))
      #Calculate total SE
      TotSE = np.sqrt(SumSE**2 + SE**2)
      #If only one replicate, use that SE
      if len(Size.loc[Size ==1])>1:
         Dose = Size.loc[Size==1].index.values
         for i in range(len(Dose)):
            TotSE.loc[Dose[i]]=AverageSurv.loc[AverageSurv['Dose [Gy]']==Dose[i],'SurvWMulti SE'].values[0]
      else:
         pass
      
   elif Type =='proton':
      #Calculate error and remove columns Factor, MU prescribed and delivered
      Average = AverageSurv.groupby(['Position','Dose [Gy]']).mean()
      STD = AverageSurv.groupby(['Position','Dose [Gy]']).std()
      Size = AverageSurv.groupby(['Position','Dose [Gy]']).size()
      #Calculate the sum of variances
      Sum = SquaredData.groupby(['Position','Dose [Gy]']).sum()
      Sum.index = Size.index
      SumSE = np.sqrt(Sum['SurvWMulti SE']/Size**2)
      SumDoseSE = np.sqrt(Sum['Corrected Dose SE']/Size**2)
      SE = STD['SurvWMulti']/np.sqrt(Size)
      
      Average['Corrected Dose SE']= np.sqrt(SumDoseSE**2+(STD['Corrected Dose [Gy]']/np.sqrt(Size))**2)
      
      #Calculate total SE
      TotSE= np.sqrt(SumSE**2 + SE**2)

      #If only one replicate, use that SE
      if len(Size.loc[Size ==1]):
         Indexes = Size.loc[Size ==1].index
         for i in range(len(Indexes)):

            TotSE.loc[Indexes[i]]=AverageSurv.loc[(AverageSurv['Position']==Indexes[i][0])&(AverageSurv['Dose [Gy]']==Indexes[i][1]),'SurvWMulti SE'].values[0]
            Average.loc[Indexes[i],'Corrected Dose SE']=AverageSurv.loc[(AverageSurv['Position']==Indexes[i][0])&(AverageSurv['Dose [Gy]']==Indexes[i][1]),'Corrected Dose SE'].values[0]
            
      else:
         pass
   else:
      print('Something went wrong')
   
   #Getting the total Standard error 
   Average['SurvWMulti SE'] = TotSE

   return Average

def AverageSurvivalWMulti(MeanSample,MeanControl,Multiplicity,PlatingEfficiency):
   #Calculating Surviving Fraction
   MeanSample = SurvivingFraction(MeanSample,PlatingEfficiency, MeanControl)
   #Get the dose as a column, not index
   MeanSample = MeanSample.reset_index()
   #Calculating Surviving Fraction with Multiplicity correction
   MeanSample = SurvCorrMulti(MeanSample,Multiplicity['Count'].values[0])
   #Set the same labels as for multi-corrected samples, but using surv values for control
   MeanControl['SurvWMulti'] = MeanControl['Surv']
   MeanControl['SurvWMulti SE'] = MeanControl['Surv SE']
   #Appedning controls to irradiated samples
   Data = MeanSample.append(MeanControl, ignore_index = True)
   return Data.sort_values(by = ['Dose [Gy]']).drop(columns=['Cells/dish','Count','SE Count'])

def DoseCorrection(Samples,DoseCorrectionFactor):
   #Samples is a dataframe, function which returns the dataframe with the corrected Dose values. 1.411 is the calibration factor for the ionchamber.
   #DoseCorrectionFactor is the correction factor found by analyzing the GafChromic films
   Factor = Samples['Factor']
   MU = Samples['MU delivered']
   CorrectedDose = MU*1.411*Factor*DoseCorrectionFactor
   Samples['Corrected Dose [Gy]'] = CorrectedDose 
   return Samples

def LQModel(Dose, alpha, beta):
   return np.exp(-alpha*Dose - beta*Dose**2)

def LQModelHighLET(Dose, alpha):
   return np.exp(-alpha*Dose)

def LogLQModel(Dose, alpha, beta):
   return (-alpha*Dose-beta*Dose**2)

def LQModelFit(Dose,Surv,SurvSE,AllDoses,AllSurv,AllSurvSE,weighted):
   #Finding the alpha and beta values from the experimental data
   if weighted == 'yes':
      Fitparam,FitPCov = optimize.curve_fit(LogLQModel,Dose,np.log(Surv),sigma = np.log(SurvSE))
      Pred = LogLQModel(AllDoses,Fitparam[0],Fitparam[1])
      Real = np.log(AllSurv)
      RSquared = sk.r2_score(Real,Pred, sample_weight=np.log(AllSurvSE))
   elif weighted == 'no':
      Fitparam,FitPCov = optimize.curve_fit(LogLQModel,Dose,np.log(Surv))
      Pred = LogLQModel(AllDoses,Fitparam[0],Fitparam[1])
      Real = np.log(AllSurv)
      RSquared = sk.r2_score(Real,Pred)

   else:
      print('Did not pass option for weighting of fit')

   Alpha = [Fitparam[0],np.sqrt(FitPCov[0][0])]
   Beta = [Fitparam[1],np.sqrt(FitPCov[1][1])]

   return Alpha,Beta, RSquared

def LogLQModelHighLET(Dose, alpha):
   return (-alpha*Dose)

def LQModelFitHighLET(Dose,Surv,SurvSE,AllDoses,AllSurv,AllSurvSE,weighted):
   #Finding the alpha and beta values from the experimental data
   if weighted == 'yes':
      Fitparam,FitPCov = optimize.curve_fit(LogLQModelHighLET,Dose,np.log(Surv),sigma = np.log(SurvSE))
      Pred = LogLQModelHighLET(AllDoses,Fitparam[0])
      Real = np.log(AllSurv)
      RSquared = sk.r2_score(Real,Pred, sample_weight = np.log(AllSurvSE))

   elif weighted == 'no':
      Fitparam,FitPCov = optimize.curve_fit(LogLQModelHighLET,Dose,np.log(Surv))
      Pred = LogLQModelHighLET(AllDoses,Fitparam[0])
      Real = np.log(AllSurv)
      
      RSquared = sk.r2_score(Real,Pred)
   else:
      print('Did not pass option for weighting of fit')

   Alpha = [Fitparam[0],np.sqrt(FitPCov[0][0])]
   return Alpha,RSquared

def SESurvivalWithMulti(M,Sample):
   #Calculating SE for Survival with Multiplicity correction
   a = (Sample['Surv SE'])/np.sqrt(M**2 - 4*(M - 1)*Sample['Surv'])
   b =  (M - 2*(M-1)*Sample['Surv'])/np.sqrt(M**2 - 4*(M - 1)*Sample['Surv'])
   c = 2*(M-1)**2
   DeltaM = 0.03 #Ref Sollund's Master
   return np.sqrt(a**2 +(((-1+b)/c)*DeltaM)**2)

def SurvCorrMulti(Sample,M):
   #Calcualting Survivalfraction with Multiplicity correction
   #But using Surv values if Surv >=1
   for i in range(len(Sample)):
      if Sample.loc[i,'Surv']>=1:
         Sample.loc[i,'SurvWMulti'] = Sample.loc[i,'Surv']
         Sample.loc[i,'SurvWMulti SE'] = Sample.loc[i,'Surv SE']
      else:
         Sample.loc[i,'SurvWMulti'] = (M-np.sqrt(M**2-4*(M-1)*Sample.loc[i,'Surv']))/(2*(M-1))
         Sample.loc[i,'SurvWMulti SE'] = SESurvivalWithMulti(M,Sample.loc[i,:])
   return Sample

def SurvivingFraction(Sample,PE,MeanControl):
   #Calculate by sending in average sample values
   Sample['Surv'] = Sample['Count']/(Sample['Cells/dish']*PE)
   DeltaPE = MeanControl['SE Count']/MeanControl['Cells/dish']
   Sample['Surv SE'] = np.sqrt((-Sample['Count']/(Sample['Cells/dish']*PE**2)*DeltaPE)**2 + (Sample['SE Count']/(Sample['Cells/dish']*PE))**2)
   return Sample


def SurvWMultiAllSamples(Samples,MeanControl, M):
   #Calculate by sending in average sample values
   PE = MeanControl['Count']/MeanControl['Cells/dish']
   Samples['Surv'] = Samples['Count']/(Samples['Cells/dish']*PE)
   DeltaPE = MeanControl['SE Count']/MeanControl['Cells/dish']
   #Different SE than when using average values for count
   Samples['Surv SE'] = Samples['Count']/(Samples['Cells/dish']*PE**2)*DeltaPE
   SurWMultiAll = SurvCorrMulti(Samples,M)
   return SurWMultiAll.sort_values(by=['Dose [Gy]','SurvWMulti'])#.drop(columns = ['Cells/dish','Count'])

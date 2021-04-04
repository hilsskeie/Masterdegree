"""
Reading and sorting count from file into lits. Needs filepath and filename in one string as input
"""
import pandas as pd
import os

def ReadCount(filepath):
    Datafile = pd.read_excel(filepath,['Samples','Control'])
    Samples = Datafile['Samples']
    ControlData = Datafile['Control']
    Multiplicity = ControlData[ControlData.Type=='Multiplicity']
    Control = ControlData[ControlData.Type=='Control']
    return(Samples,Control,Multiplicity)


def ReadFolder(Filepath):
    #Get the files within the director
    Files = os.listdir(Filepath)
    AnalyzeFiles=[]

    #Skipping files that is not CSV
    n = len(Files)
    for i in range(n):
        if Files[i].endswith('.csv'):
            Survdata = Files[i] 
            AnalyzeFiles.append(Survdata)
        else:
            pass
       
    return AnalyzeFiles

def ReadSurv(AverageSurvFilenames,AllSurvFilenames,Filepath):
    #Create empty Dataframe to put all the results in
    if 'proton' in Filepath:
        AverageSurv = pd.DataFrame(columns=['Position','Dose [Gy]','Corrected Dose [Gy]','Corrected Dose SE','Surv','Surv SE','SurvWMulti','SurvWMulti SE'])
        AllSurv = pd.DataFrame(columns=['Position','Dose [Gy]','Corrected Dose [Gy]','Surv','Surv SE','SurvWMulti','SurvWMulti SE'])
    else:
        AverageSurv = pd.DataFrame(columns=['Dose [Gy]','Surv','Surv SE','SurvWMulti','SurvWMulti SE'])
        AllSurv = pd.DataFrame(columns=['Dose [Gy]','Surv','SurvWMulti'])
    
    for i in range(len(AllSurvFilenames)):
        #Read file
        AverageData = pd.read_csv(Filepath+'Average/'+AverageSurvFilenames[i], index_col=[0]) 
        AllData = pd.read_csv(Filepath+'All/'+AllSurvFilenames[i], index_col=[0]) 
        #Add data from file to Dataframe
        AverageSurv = AverageSurv.append(AverageData, ignore_index= True, sort = False)
        AllSurv = AllSurv.append(AllData, ignore_index=True, sort = False)
    return AverageSurv,AllSurv

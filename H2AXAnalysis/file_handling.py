import pandas as pd
import os

def Write_to_file(Dataframe, filepath):
    filepathinformation = filepath.split('/') 
    fileinformation = filepathinformation[-1].split('_')
    date = fileinformation[0]

    try:
        #Try to open file
        data_in_file = pd.read_csv(filepath, index_col = 0)

        #Check if replicate already in file
        if 'xrays' in fileinformation:
            columns_to_check = ['Time', 'Dose', 'Replicate']
            appended_data_in_file = data_in_file.append(Dataframe)
            new_file = appended_data_in_file.drop_duplicates(subset =columns_to_check, keep ='first')
            new_file.to_csv(filepath)
            
        elif 'protons' in fileinformation:
            columns_to_check = ['Position','Time', 'Dose', 'Replicate']
            appended_data_in_file = data_in_file.append(Dataframe)
            new_file = appended_data_in_file.drop_duplicates(subset =columns_to_check, keep ='first')
            new_file.to_csv(filepath)
 
    except FileNotFoundError:
        #If file does not exists, create it
        Dataframe.to_csv(filepath)

def Info_From_Filename(Filename,Filepath_Result):
    #Getting information from filename
    filename = Filename
    fileinformation = filename.split('_')
    date = fileinformation[0]
   
    cell_line = fileinformation[1]
    irradiation = fileinformation[2]
    
    if 'proton' in irradiation:
        position = fileinformation[3]
        dose = fileinformation[4]
        number_of_fractions = fileinformation[5]
        time = fileinformation[6]
        replicate = fileinformation[7]

    else:
        dose = fileinformation[3]
        number_of_fractions = fileinformation[4]
        time = fileinformation[5]
        replicate = fileinformation[6][:1]

    #Want to format half hour
    if '05' in time:
        time = '0.5H'
    else:
        pass

    #Create Dataframe with info from file
    if 'xray' in irradiation:
        Dataframe = pd.DataFrame({'Time': float(time[:-1]), 'Dose': dose[:-2], 'Fractions':number_of_fractions, 'Replicate': int(replicate)}, index = [0])
    elif 'proton' in irradiation:
        Dataframe = pd.DataFrame({'Position': int(position[-1]), 'Time': float(time[:-1]), 'Dose': dose[:-2], 'Fractions':number_of_fractions, 'Replicate': int(replicate)}, index = [0])
    else:
        print('Not able to create file to save data')
    #Creating new filepath and name for resultfile
    New_file_name = Filepath_Result +date +'_'+ cell_line +'_'+irradiation+'_results.csv'

    return Dataframe,New_file_name


def Read_folder(Filepath):
    #Get the files within the director
    Files = os.listdir(Filepath)
    AnalyzeFiles=[]

    #Skipping files that is not CSV and saving the one that endes with _sample
    n = len(Files)
    for i in range(n):
        if Files[i].endswith('_sample.csv'):
            Replicate = Files[i].replace('_sample.csv','')         
            AnalyzeFiles.append(Replicate)
        else:
            pass
       
    return AnalyzeFiles

def Read_Sample_and_BC(FilePath,Filename):
    Fileinformation = Filename.split('_')
    FilenameSample = FilePath+Filename+'_sample.csv'
    FilenameBC = FilePath+Filename+'_barcode.csv'
    Sample = pd.read_csv(FilenameSample).drop(columns=['FL3-A','FL3-H','FL1-H','FL2-H','FL4-H','Width', 'Time'])
    BC = pd.read_csv(FilenameBC).drop(columns=['FL3-A','FL3-H','FL1-H','FL2-H','FL4-H','Width', 'Time'])

    if 0 in Sample['FL2-A'].values:
        print('Skipped zeros from dataframe from sample file: ', Filename)
        print('Skipped line: \n', Sample.loc[Sample['FL2-A']==0])

    elif 0 in BC['FL2-A'].values:
        print('Skipped zeros from dataframe from barcode file: ', Filename)
        print('Skipped line: \n', BC.loc[BC['FL2-A']==0])
    else:
        pass
    Samples_without_zeros = Sample.loc[Sample['FL2-A']!=0]
    BC_witout_zeros = BC.loc[BC['FL2-A']!=0]

    return Samples_without_zeros,BC_witout_zeros,Fileinformation

def Read_Cell_cycle_dist(Filepath):
    Result_data = pd.read_csv(Filepath, index_col = 0)
    Distribution_data = Result_data.drop(columns = ['Sample median full', 'Sample error full', 'BC median full', 'BC error full', 'Ratio full', 'Ratio error full','Ratio G1G2','Ratio error G1G2','Ratio G1','Ratio error G1','Ratio G2', 'Ratio error G2','Ratio S', 'Ratio error S', 'Sample median G1', 'Sample error G1', 'BC median G1', 'BC error G1', 'Sample median G2', 'Sample error G2', 'BC median G2', 'BC error G2', 'Sample median S', 'Sample error S', 'BC median S', 'BC error S', 'Replicate'])
    Distribution_Sample = Distribution_data.loc[Distribution_data.Dose!=0].drop(columns = ['BC G2 amount', 'BC G1 amount', 'BC S amount'])
    Distribution_BC = Distribution_data.drop(columns = ['Sample G2 amount', 'Sample G1 amount', 'Sample S amount'])
    Distribution_Control = Distribution_data.loc[Distribution_data.Dose==0]
    return Distribution_Sample, Distribution_BC, Distribution_Control


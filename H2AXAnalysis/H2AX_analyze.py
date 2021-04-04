from H2AX_functions import *
from file_handling import *
from plotting import Mark_cell_cycle_limits, Plotting_Results,Plotting_Results_Proton, Plotting_Cycle_phase, Plotting_Cycle_phase_proton, Mark_cell_cycle

import pandas as pd
import sys


"""
README
The first argument is the filepath of the choosen file to analyze. Then there are X options to add after the filepath:
filter  - Filter the cell cycle, meaning you will mark the start and end of the cell cycle and get filtered .csv files.
cycle - You will get a plot with the different phases and the distribution of cells. This will also be saved as a .csv file.

If you do not add anything after the filepath the full analyze will run, inlcuding markin the phases of the cell cycle and the plotting of median and cycle results.
XX insert file format note XX
"""
#Get filepath as input from command line
Filepath = sys.argv[1]

#Check if other commands are given
try:
    Argument = sys.argv[2]
except IndexError:
    Argument = ' '


if Filepath[-1]=='/':    
    #Get filenames
    ListOfFilenames = Read_folder(Filepath)
    count = len(ListOfFilenames)-1

    if Argument =='filter':
        for Filename in ListOfFilenames:
            Fileinformation = Filename.split('_')
            #Read files
            Sample = pd.read_csv(Filepath+Filename+'_sample.csv')
            BC = pd.read_csv(Filepath+Filename+'_barcode.csv')

            #Plot files and set Cell cycle limits
            TitleBC = Fileinformation[-4]+ ' - ' + Fileinformation[-2] + ' - Replicate ' + Fileinformation[-1] + ' - Barcode'
            TitleSample = Fileinformation[-4]+ ' - ' + Fileinformation[-2] + ' - Replicate ' + Fileinformation[-1] + ' - Sample'
            
            #Create new filenames, do not wnat to overwrite the old files
            new_filename_Sample = Filepath + 'Filtered/' + Filename + '_filtered_sample.csv'
            new_filename_BC = Filepath +'Filtered/' +Filename + '_filtered_barcode.csv'
            
            #Mark Start and Stop of cell cycle
            Sample_Cell_Cycle_limits = Mark_cell_cycle(Sample['FL2-A'], TitleSample)
            plt.close()
            BCCell_Cycle_limits = Mark_cell_cycle(BC['FL2-A'],TitleBC)
            plt.close()
            print(count, ' more replicates to set limits')
            count -= 1

            #Filter out cells that is not in the marked cycle
            Sample_filtered = CellCycleFiltering(Sample,Sample_Cell_Cycle_limits)
            BC_filtered = CellCycleFiltering(BC,BCCell_Cycle_limits)
            
            #Save the filtered samples
            try:
                Sample_filtered.to_csv(new_filename_Sample)
                BC_filtered.to_csv(new_filename_BC)
            except FileNotFoundError:
                os.mkdir(Filepath+'Filtered/')
                Sample_filtered.to_csv(new_filename_Sample)
                BC_filtered.to_csv(new_filename_BC)
    else:
        for Filename in ListOfFilenames:        
            #Read files
            Sample,BC,Fileinformation = Read_Sample_and_BC(Filepath,Filename)
            if 'filtered' in Fileinformation:
                if 'xrays' in Fileinformation:
                    TitleBC = Fileinformation[-5]+ ' - ' + Fileinformation[-3] + ' - Replicate ' + Fileinformation[-2] + ' - Barcode'
                    TitleSample = Fileinformation[-5]+ ' - ' + Fileinformation[-3] + ' - Replicate ' + Fileinformation[-2] + ' - Sample'
                else:
                    TitleBC = Fileinformation[-6]+ ' - '+ Fileinformation[-5]+ ' - ' + Fileinformation[-3] + ' - Replicate ' + Fileinformation[-2] + ' - Barcode'
                    TitleSample = Fileinformation[-6]+ ' - ' + Fileinformation[-5]+ ' - ' + Fileinformation[-3] + ' - Replicate ' + Fileinformation[-2] + ' - Sample'
            else:
                if 'xrays' in Fileinformation:
                    TitleBC = Fileinformation[-4]+ ' - ' + Fileinformation[-2] + ' - Replicate ' + Fileinformation[-1] + ' - Barcode'
                    TitleSample = Fileinformation[-4]+ ' - ' + Fileinformation[-2] + ' - Replicate ' + Fileinformation[-1] + ' - Sample'
                else:
                    TitleBC = Fileinformation[-5]+ ' - ' + Fileinformation[-4]+ ' - ' + Fileinformation[-2] + ' - Replicate ' + Fileinformation[-1] + ' - Barcode'
                    TitleSample = Fileinformation[-5]+ ' - ' + Fileinformation[-4]+ ' - ' + Fileinformation[-2] + ' - Replicate ' + Fileinformation[-1] + ' - Sample'
            
            #Plot files and set G1 and G2 limits
            Sample_Cell_Cycle_limits = Mark_cell_cycle_limits(Sample['FL2-A'], TitleSample)
            plt.close()
            BCCell_Cycle_limits = Mark_cell_cycle_limits(BC['FL2-A'],TitleBC)
            plt.close()
            print(count, ' more replicates to set limits')
            count -= 1

            #Calculating Median data
            Median_data = H2AX_Find_Median(Sample,Sample_Cell_Cycle_limits,BC,BCCell_Cycle_limits)

            #Create Results folder
            Filepathwithresults = Filepath + '/Results/'
            try:
                os.mkdir(Filepathwithresults)
            except FileExistsError:
                pass
            #Get info from filename
            Info_from_filenameDF, Newfilename = Info_From_Filename(Filename,Filepathwithresults)
            #Merge median data and info from file name
            ResultDF = pd.concat([Info_from_filenameDF,Median_data], axis=1, join='inner')

            #Write to file
            Write_to_file(ResultDF, Newfilename)
    

elif 'results' in Filepath:
    Average_ratio = Normalize_to_Controls(Filepath)
   
    FilepathInfo= Filepath.split('/')
    Filename = FilepathInfo[-1]
    NewDirectory = '/'.join(FilepathInfo[:-2])+'/Results/'
    title = 'PI fluorescence intensity - '+Filename[:6]

    if Argument == 'cycle':
        Distribution_Sample, Distribution_BC, Distribution_Control = Read_Cell_cycle_dist(Filepath)

        if 'xrays' in Filepath:
            #Cell cycle plots
            G1_sample,S_sample,G2_sample = Phaseanalyze(Distribution_Sample, 'Sample', 'xray')
            G1_control,S_control,G2_control = Phaseanalyze(Distribution_Control, 'Sample', 'xray')
            Plotting_Cycle_phase(G1_sample,G2_sample,S_sample,NewDirectory,title)
        
        elif 'protons' in Filepath:
            G1_sample,S_sample,G2_sample = Phaseanalyze(Distribution_Sample, 'Sample', 'proton')
            G1_control,S_control,G2_control = Phaseanalyze(Distribution_Control, 'Sample', 'proton')

            # Plotting_Cycle_phase_proton(G1_sample,G2_sample,S_sample,NewDirectory,title)
            G1 = pd.concat([G1_control,G1_sample], axis = 0)
            S = pd.concat([S_control,S_sample], axis = 0)
            G2 = pd.concat([G2_control,G2_sample], axis =0)
            Plotting_Cycle_phase_proton(G1,G2,S,NewDirectory,title)
            
        else:
            print('Did not analyze files, unknown irradiation type')
            
        G1_data = pd.concat([G1_control,G1_sample], axis = 0)
        S_data = pd.concat([S_control,S_sample], axis = 0)
        G2_data = pd.concat([G2_control,G2_sample], axis =0)
        
        #Save cycle distribtuion to file
        G1_data.to_csv(NewDirectory+'G1_data.csv')
        S_data.to_csv(NewDirectory+'S_data.csv')
        G2_data.to_csv(NewDirectory+'G2_data.csv')


    else:
        #Median plot
        if 'xrays' in Filepath:
            Plotting_Results(Average_ratio,NewDirectory,'G1', 'X-ray irradiation - G1 cells' )
            Plotting_Results(Average_ratio,NewDirectory,'G2','X-ray irradiation - G2 cells')
            Plotting_Results(Average_ratio,NewDirectory,'S','X-ray irradiation - S-phase cells')
            Plotting_Results(Average_ratio,NewDirectory,'full','X-ray irradiation - All cells')
        elif 'protons' in Filepath:
            Plotting_Results_Proton(Average_ratio,NewDirectory,'G1', 'Proton irradiation - G1 cells' )
            Plotting_Results_Proton(Average_ratio,NewDirectory,'G2','Proton irradiation - G2 cells')
            Plotting_Results_Proton(Average_ratio,NewDirectory,'S','Proton irradiation - S-phase cells')
            Plotting_Results_Proton(Average_ratio,NewDirectory,'full','Proton irradiation - All cells')
        else:
            pass

        Average_ratio.to_csv(NewDirectory+'H2AX_results.csv')

else:
    print('Did not analyze, wrong input - maybe missing a /? ')
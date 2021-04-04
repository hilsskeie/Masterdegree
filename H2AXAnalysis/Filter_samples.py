from file_handling import Read_Sample_and_BC, Read_folder
from H2AX_functions import CellCycleFiltering
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os


def Mark_thershold(Dataframe, title):
        """
        Give input for limits by clicking the plot
 
        The keyboard can also be used to select points in case your mouse does not have one or more of the buttons. 
        The delete and backspace keys act like right clicking (i.e., remove last point), the enter key terminates 
        input and any other key (not already used by the window manager) selects a point.
        """
        plt.figure()
        plt.scatter(Sample['FL2-A'],Sample['FL1-A'], alpha=0.5, s=0.3, color = 'midnightblue', label = 'Sample')
        plt.yscale('log')
        plt.title(title)
        marker_input = plt.ginput(2, timeout = -1)
        return marker_input
        plt.show()

Filepath = sys.argv[1]
#Get filenames
ListOfFilenames = Read_folder(Filepath)
count = len(ListOfFilenames)-1

for Filename in ListOfFilenames:
    Fileinformation = Filename.split('_')
    #Read files
    Sample = pd.read_csv(Filepath+Filename+'_sample.csv')
    BC = pd.read_csv(Filepath+Filename+'_barcode.csv')

    #Plot files and set Cell cycle limits
    
    TitleSample = Fileinformation[-4]+ ' - ' + Fileinformation[-2] + ' - Replicate ' + Fileinformation[-1] + ' - Sample'
    Sample_thershold = Mark_thershold(Sample, TitleSample)
    plt.close()
    FilteredSampleFL1 = Sample.loc[Sample['FL1-A']>=Sample_thershold[0][1]]
    FilteredSampleFL2 = FilteredSampleFL1.loc[Sample['FL2-A']>=Sample_thershold[0][0]]
    SampleFL2 = FilteredSampleFL2.loc[FilteredSampleFL2['FL2-A']<=Sample_thershold[1][0]]


    new_filename_Sample = Filepath + 'Thershold/' + Filename + '_thershold_sample.csv'
    new_filename_BC = Filepath +'Thershold/' +Filename + '_thershold_barcode.csv'

    try:
        SampleFL2.to_csv(new_filename_Sample)
        BC.to_csv(new_filename_BC)
    except FileNotFoundError:
        os.mkdir(Filepath+'Thershold/')
        SampleFL2.to_csv(new_filename_Sample)
        BC.to_csv(new_filename_BC)
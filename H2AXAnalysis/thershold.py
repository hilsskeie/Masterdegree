from file_handling import Read_folder,Read_Sample_and_BC
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import sys 

#Setting font
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'cm'
mpl.rcParams['font.sans-serif'] = 'Times New Roman'
mpl.rcParams['font.family'] = 'Times New Roman'

Filepath = sys.argv[1]


Files = Read_folder(Filepath)
for File in Files:
  Sample,BC,fileinfo = Read_Sample_and_BC(Filepath,File)
  plt.figure(figsize=(7,7))
  #Plot FL1 and FL2 for BC and Sample
  ax1 = plt.subplot(221)
  ax1.scatter(Sample['FL2-A'],Sample['FL1-A'], alpha=0.3, s=0.3, color = 'midnightblue', label = 'Sample')
  ax1.scatter(BC['FL2-A'],BC['FL1-A'], alpha=0.3, s=0.3, color = 'maroon', label = 'BC')
  ax1.set_yscale('log')
  ax1.set_title('Sample and BC')
  ax1.set_ylabel('FL1-A')
  ax1.set_xlabel('FL2-A')
  leg = plt.legend(loc='upper left')
  leg.legendHandles[0]._sizes = [8]
  leg.legendHandles[1]._sizes = [8]
  for lh in leg.legendHandles: 
    lh.set_alpha(1)

  #plot Cell cycle
  ax2 = plt.subplot(222)
  ax2.hist(Sample['FL2-A'], bins = 500, color = 'midnightblue')
  ax2.set_title('Cell cycle histogram Sample')
  ax2.set_xlabel('FL2-A')
  print(File)
  Title = fileinfo[0] + ' - ' + fileinfo[5] + ' - ' +fileinfo[3] + ' - '+ fileinfo[-4] + ' - '+ fileinfo[-2]
  plt.suptitle(Title)
 
  #Plot FL1 and FL2 for just sample
  ax3 = plt.subplot(212)
  ax3.scatter(Sample['FL2-A'],Sample['FL1-A'], alpha=0.5, s=0.3, color = 'midnightblue', label = 'Sample')
  ax3.set_yscale('log')
  ax3.set_title('Sample')
  ax3.set_ylabel('FL1-A')
  ax3.set_xlabel('FL2-A')
  # point = plt.ginput(1, timeout = -1)
  # plt.close()
  # FilteredSample = Sample.loc[(Sample['FL2-A']>=point[0][0]) & (Sample['FL1-A']>=point[0][1])]
  # plt.figure()
  # plt.scatter(FilteredSample['FL2-A'],FilteredSample['FL1-A'], alpha=0.5, s=0.3, color = 'midnightblue', label = 'Sample')
  # ax3.set_yscale('log')
  # ax3.set_title('Sample')
  # ax3.set_ylabel('FL1-A')
  # ax3.set_xlabel('FL2-A')
  NewFilePath = Filepath+'/Threshold/'
  plt.savefig(NewFilePath+File+'_thershold.pdf')
 
  #FilteredSample = Sample.loc[(Sample['FL1-A']<x0) and (Sample['FL1-A']>x1),(Sample['FL2-A']<y0) and(Sample['FL2-A']>y1)]

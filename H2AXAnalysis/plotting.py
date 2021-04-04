import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import time
import pandas as pd

#list of colors to plot
colors1 = ['maroon']
colors2 = ['maroon','midnightblue']
colors3 = ['maroon','green','midnightblue']
colors4 = ['maroon','orange','green','midnightblue']
colors5 = ['maroon','orangered','orange','green','midnightblue']


#Setting font
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'cm'
mpl.rcParams['font.sans-serif'] = 'Times New Roman'
mpl.rcParams['font.family'] = 'Times New Roman'



def Mark_cell_cycle(Dataframe, title):
        """
        Give input for limits by clicking the plot

        The keyboard can also be used to select points in case your mouse does not have one or more of the buttons. 
        The delete and backspace keys act like right clicking (i.e., remove last point), the enter key terminates 
        input and any other key (not already used by the window manager) selects a point.
        """
        plt.figure()
        plt.hist(Dataframe,bins=500,histtype='step')
        plt.title(title)
        marker_input = plt.ginput(2, timeout = -1)
        return marker_input
        plt.show()

def Mark_cell_cycle_limits(Dataframe, title):
        """
        Give input for limits by clicking the plot

        The keyboard can also be used to select points in case your mouse does not have one or more of the buttons. 
        The delete and backspace keys act like right clicking (i.e., remove last point), the enter key terminates 
        input and any other key (not already used by the window manager) selects a point.
        """
        minvalue = Dataframe.min()
        maxvalue = Dataframe.max()
        plt.figure()
        plt.hist(Dataframe,bins=250,histtype='step')
        plt.xlim(minvalue*0.6, maxvalue*1.2)
        plt.title(title)
        marker_input = plt.ginput(4, timeout = -1)

        #Only x-component is relevant
        G1Limits = [marker_input[0][0],marker_input[1][0]]
        G2Limits = [marker_input[2][0],marker_input[3][0]]
        return G1Limits, G2Limits
        plt.show()



def Plotting_Results(average_ratio,filepath,Phase,title):
        Ratio = 'Ratio '+ Phase
        Error = 'Error '+ Phase
        #Get positions, timepoints and dose. This assumes one dose
        timepoints = average_ratio.index.get_level_values(0).unique()
        dose = average_ratio.index.get_level_values(1).unique()
        if len(dose) == 1:
                colors = colors1
        elif len(dose) == 2:
                colors = colors2
        elif len(dose) == 3:
                colors = colors3
        elif len(dose) == 4:
                colors = colors4
        #Create figure
        fig,ax = plt.subplots(figsize=(6.5,5))
        box = ax.get_position()
        ax.set_position([box.x0,box.y0,box.width*0.9,box.height])
        for i in range(len(dose)):
                doses_to_plot = average_ratio.loc[average_ratio.index.get_level_values(level=1)==dose[i]]
                doses_to_plot.index = doses_to_plot.index.droplevel(1)
                plt.plot(doses_to_plot.index.get_level_values(level=0), doses_to_plot[Ratio],'o--', color = colors[i], linewidth=0.8, label= dose[i])
                plt.errorbar(doses_to_plot.index.get_level_values(level=0), doses_to_plot[Ratio], yerr= doses_to_plot[Error], fmt='none', elinewidth = 0.8, ecolor=colors[i],capthick=1, capsize=2, label ='')
                plt.title(title, fontsize = 14)
                plt.ylabel(r'Relative $\gamma-$H2AX fluorescence intensity', fontsize = 13)
                plt.xlabel('Time [h]', fontsize = 13)
                plt.xticks(timepoints)
                plt.grid(which='major',linestyle='--', linewidth='0.5', color='lightgray')
                leg = plt.legend(bbox_to_anchor = (1.24,0.62),fontsize = 12)
                leg.set_title('Dose [Gy]', prop = {'size':'large'})
        plt.savefig(filepath+'H2AX_'+ Phase+'_intensity.pdf')
        # plt.show()


def Plotting_Results_Proton(average_ratio,filepath,Phase,title):
        Ratio = 'Ratio '+ Phase
        Error = 'Error '+ Phase
        #Get positions, timepoints and dose. This assumes one dose
        positions = average_ratio.index.get_level_values(0).unique()
        timepoints = average_ratio.index.get_level_values(1).unique()
        dose = average_ratio.index.get_level_values(2).unique()
        
        #Finding the number of positions and timepoints
        number_of_positions = len(positions.values)
        number_of_timepoints = len(timepoints)

        #Create figure
        fig,ax = plt.subplots(figsize=(6,4.7))
        box = ax.get_position()
        ax.set_position([box.x0,box.y0,box.width*0.9,box.height])
        # axins = ax.inset_axes([0.6,0.25,0.35,0.35])
        for i in range(number_of_timepoints):
                positions_to_plot = average_ratio.loc[average_ratio.index.get_level_values(level=1)==timepoints[i]]
                positions_to_plot.index = positions_to_plot.index.droplevel(1)

                plt.plot(positions_to_plot.index.get_level_values(level=0), positions_to_plot[Ratio],'o--', color = colors2[i], linewidth=0.8, label= timepoints[i])
                plt.errorbar(positions_to_plot.index.get_level_values(level=0), positions_to_plot[Ratio], yerr= positions_to_plot[Error], fmt='none', elinewidth = 0.8, ecolor=colors2[i],capthick=1, capsize=2, label ='')
                
                # axins.plot(positions_to_plot.index.get_level_values(level=0), positions_to_plot[Ratio],'o--', color = colors2[i], linewidth=0.8, label= timepoints[i])
                # axins.errorbar(positions_to_plot.index.get_level_values(level=0), positions_to_plot[Ratio], yerr= positions_to_plot[Error], fmt='none', elinewidth = 0.8, ecolor=colors2[i],capthick=1, capsize=2, label ='')
                # axins.set_ylim([1,2.5])
                # axins.set_xticks([1,2,3,4,5])
                # axins.grid(which='major',linestyle='--', linewidth='0.5', color='lightgray')
                
                plt.title(title + ' - '+str(dose.values[0])+' Gy', fontsize = 16)
                plt.ylabel(r'Relative $\gamma$-H2AX fluorescence intensity', fontsize = 14)
                plt.xlabel('Position', fontsize = 14)
                plt.xticks([1,2,3,4,5],fontsize=13)
                plt.yticks(fontsize=13)
                plt.grid(which='major',linestyle='--', linewidth='0.5', color='lightgray')
                leg = plt.legend(bbox_to_anchor = (1.26,0.62),fontsize = 13)
                leg.set_title('Time [h]', prop = {'size':'x-large'})
        plt.savefig(filepath+'H2AX_'+ Phase+'_intensity.pdf')
        plt.show()


def Plotting_Cycle_phase(G1,G2,S,Filepath, Title):
        Time = G1.index.get_level_values(0).unique().tolist()
        Dose = G1.index.get_level_values(1).unique().tolist()

        #Setting colors based on number of doses
        if len(Dose) == 1:
                colors = colors1
        elif len(Dose) == 2:
                colors = colors2
        elif len(Dose) == 3:
                colors = colors3
        elif len(Dose) == 4:
                colors = colors4
        elif len(Dose) == 5:
                colors = colors5
        else:
                print('Do not know which colors to use')
        newfilepath = Filepath + 'Phase_plot.pdf'
                
        plt.figure(figsize = (9.4,6))
        grid = gridspec.GridSpec(nrows = 1, ncols = 3, left = 0.065, right = 0.99)
        grid.update(wspace = 0.06)
        
        ax1 = plt.subplot(grid[0])
        ax2 = plt.subplot(grid[1])
        ax3 = plt.subplot(grid[2])
        box1 = ax1.get_position()
        box2 = ax2.get_position()
        box3 = ax3.get_position()
        ax1.set_position([box1.x0+0.01,box1.y0,box1.width,box1.height])
        ax2.set_position([box2.x0+0.005,box2.y0,box2.width,box2.height])
        ax3.set_position([box3.x0+0.005,box3.y0,box3.width,box3.height])
        for i in range(len(Dose)):
                G1_to_plot = G1.loc[G1.index.get_level_values(level=1)==Dose[i]]
                G1_to_plot.index = G1_to_plot.index.droplevel(1)
                G2_to_plot = G2.loc[G2.index.get_level_values(level=1)==Dose[i]]
                G2_to_plot.index = G2_to_plot.index.droplevel(1)
                S_to_plot = S.loc[S.index.get_level_values(level=1)==Dose[i]]
                S_to_plot.index = S_to_plot.index.droplevel(1)
                if Dose[i] == 13:
                        Dose[i] = 12
                elif Dose[i] == 10:
                        Dose[i] = 9
                else: 
                        pass
                ax1.plot(G1_to_plot.index, G1_to_plot['Amount'],'o--',label = Dose[i], color = colors[i], linewidth = 1)
                ax1.errorbar(G1_to_plot.index, G1_to_plot['Amount'], yerr= G1_to_plot['Error'], fmt='none', elinewidth = 0.8, ecolor=colors[i],capthick=1, capsize=2, label ='')
                ax2.plot(S_to_plot.index, S_to_plot['Amount'],'o--',label = Dose[i], color = colors[i], linewidth = 1)
                ax2.errorbar(S_to_plot.index, S_to_plot['Amount'], yerr= S_to_plot['Error'], fmt='none', elinewidth = 0.8, ecolor=colors[i],capthick=1, capsize=2, label ='')
                ax3.plot(G2_to_plot.index, G2_to_plot['Amount'],'o--',label = Dose[i], color = colors[i], linewidth = 1)
                ax3.errorbar(G2_to_plot.index, G2_to_plot['Amount'], yerr= G2_to_plot['Error'], fmt='none', elinewidth = 0.8, ecolor=colors[i],capthick=1, capsize=2, label ='')
                ax1.set_title('G1', fontsize = 20)
                ax1.set_ylabel('$\%$ of total cells', fontsize = 16)
                ax1.set_xlabel('Time [h]', fontsize = 16)
                ax2.set_xlabel('Time [h]', fontsize = 16)
                ax3.set_xlabel('Time [h]', fontsize = 16)
                ax2.set_title('S - phase', fontsize = 20)
                ax3.set_title('G2', fontsize = 20)
                ax1.set_xticks(Time)
                ax2.set_xticks(Time)
                ax3.set_xticks(Time)
                ax1.set_xticklabels(['0.5','24','48','72'],fontsize=16)
                ax2.set_xticklabels(['0.5','24','48','72'],fontsize=16)
                ax3.set_xticklabels(['0.5','24','48','72'],fontsize=16)
                ax1.grid(which='major',linestyle='--', linewidth='0.5', color='lightgray')
                ax2.grid(which='major',linestyle='--', linewidth='0.5', color='lightgray')
                ax3.grid(which='major',linestyle='--', linewidth='0.5', color='lightgray')
                ax1.set_yticks([20,40,60,80,100])
                ax1.set_yticklabels(['20','40','60','80','100'],fontsize=16)
                ax2.set_yticks([20,40,60,80,100])
                ax2.set_yticklabels(['','','','',''],fontsize=15)
                ax3.set_yticks([20,40,60,80,100])
                ax3.set_yticklabels(['','','','',''],fontsize=15)
                ax1.set_ylim([-0.1, 100.1])
                ax2.set_ylim([-0.1, 100.1])
                ax3.set_ylim([-0.1, 100.1])

                plt.suptitle(Title, fontsize = 22)

        leg = plt.legend(bbox_to_anchor = (0.91,0.55),fontsize = 16)
        leg.set_title('Dose [Gy]', prop = {'size':'x-large'})
        plt.savefig(newfilepath)
        plt.show()

def Plotting_Cycle_phase_proton(G1,G2,S,Filepath, Title):
        Position = G1.index.get_level_values(0).unique().tolist()
        Time = G1.index.get_level_values(1).unique().tolist()
        if len(Time) == 2:
                colors = colors2
        elif len(Time) == 3:
                colors = colors3
        elif len(Time) == 4:
                colors = colors4
        else:
                print('Do not know which colors to use')

        plt.figure(figsize = (9,6))
        grid = gridspec.GridSpec(nrows = 1, ncols = 3, left = 0.065, right = 0.99)
        grid.update(wspace = 0.06)
        
        ax1 = plt.subplot(grid[0])
        ax2 = plt.subplot(grid[1])
        ax3 = plt.subplot(grid[2])

        newfilepath = Filepath + 'Phase_plot.pdf'
        Labels = ['C',1,2,3,4,5]
        box = ax1.get_position()
        box2 = ax2.get_position()
        ax1.set_position([box.x0+0.02,box.y0,box.width,box.height])
        ax2.set_position([box2.x0+0.01,box2.y0,box2.width,box2.height])
        for i in range(len(Time)):
                G1_to_plot = G1.loc[G1.index.get_level_values(1)==Time[i]]
                G1_to_plot.index = G1_to_plot.index.droplevel([1,2])
                G2_to_plot = G2.loc[G2.index.get_level_values(1)==Time[i]]
                G2_to_plot.index = G2_to_plot.index.droplevel([1,2])
                S_to_plot = S.loc[S.index.get_level_values(1)==Time[i]]
                S_to_plot.index = S_to_plot.index.droplevel([1,2])

                ax1.plot(G1_to_plot.index, G1_to_plot['Amount'],'o--',label =Time[i], color = colors[i], linewidth = 0.8)
                ax1.errorbar(G1_to_plot.index, G1_to_plot['Amount'], yerr=G1_to_plot['Error'], fmt='none', elinewidth = 1, ecolor=colors[i],capthick=1, capsize=2, label ='')
                ax2.plot(S_to_plot.index, S_to_plot['Amount'],'o--',label =Time[i], color = colors[i], linewidth = 0.8)
                ax2.errorbar(S_to_plot.index, S_to_plot['Amount'], yerr=S_to_plot['Error'], fmt='none', elinewidth = 1, ecolor=colors[i],capthick=1, capsize=2, label ='')
                ax3.plot(G2_to_plot.index, G2_to_plot['Amount'],'o--',label =Time[i], color = colors[i], linewidth = 0.8)
                ax3.errorbar(G2_to_plot.index, G2_to_plot['Amount'], yerr=G2_to_plot['Error'], fmt='none', elinewidth = 1, ecolor=colors[i],capthick=1, capsize=2, label ='')
                ax1.set_title('G1', fontsize = 20)
                ax1.set_ylabel('$\%$ of total cells', fontsize = 16)
                ax1.set_xlabel('Position', fontsize = 16)
                ax2.set_xlabel('Position', fontsize = 16)
                ax3.set_xlabel('Position', fontsize = 16)
                ax1.set_xticks([0,1,2,3,4,5])
                ax1.set_xticklabels(['C','1','2','3','4','5'], fontsize=16)
                ax2.set_xticks([0,1,2,3,4,5])
                ax2.set_xticklabels(['C','1','2','3','4','5'], fontsize=16)
                ax3.set_xticks([0,1,2,3,4,5])
                ax3.set_xticklabels(['C','1','2','3','4','5'], fontsize=16)
                ax2.set_title('S - phase', fontsize = 20)
                ax3.set_title('G2', fontsize = 20)
                ax1.grid(which='major',linestyle='--', linewidth='0.5', color='lightgray')
                ax2.grid(which='major',linestyle='--', linewidth='0.5', color='lightgray')
                ax3.grid(which='major',linestyle='--', linewidth='0.5', color='lightgray')
                ax1.set_yticks([20,40,60,80,100])
                ax1.set_yticklabels([20,40,60,80,100],fontsize=16)
                ax2.set_yticks([20,40,60,80,100])
                ax2.set_yticklabels(['','','','',''])
                ax3.set_yticks([20,40,60,80,100])
                ax3.set_yticklabels(['','','','',''])
                ax1.set_ylim([-0.1, 100.1])
                ax2.set_ylim([-0.1, 100.1])
                ax3.set_ylim([-0.1, 100.1])

                plt.suptitle(Title, fontsize = 22)
                leg = plt.legend(fontsize = 15)
                leg.set_title('Time [h]', prop = {'size':'x-large'})
        plt.savefig(newfilepath)
        plt.show()


def Plotting_average_xray(AverageData,Filepath, Phase, title):
        Ratio = 'Ratio '+ Phase
        Error = 'Error '+ Phase
        Doses = AverageData.index.get_level_values(1).unique()

        fig,ax = plt.subplots()
        axins = ax.inset_axes([0.43,0.42,0.5,0.5])
        box = ax.get_position()
        ax.set_position([box.x0,box.y0,box.width*0.95,box.height])
        for i in range(len(Doses)):
                Data_to_plot = AverageData.loc[AverageData.index.get_level_values(level=1)==Doses[i]]
                plt.plot(Data_to_plot.index.get_level_values(level=0), Data_to_plot[Ratio],'o--', color = colors4[i], linewidth=0.8, label= Doses[i])
                plt.errorbar(Data_to_plot.index.get_level_values(level=0), Data_to_plot[Ratio], yerr= Data_to_plot[Error], fmt='none', elinewidth = 0.8, ecolor=colors4[i],capthick=1, capsize=2, label ='')
                plt.title(title, fontsize=16)
                plt.ylabel(r'Relative $\gamma$-H2AX fluorescence intensity', fontsize=14)
                plt.xlabel('Time [h]', fontsize=14)
                leg = plt.legend(bbox_to_anchor = (1.2,0.65),fontsize=12)
                leg.set_title('Dose [Gy]', prop = {'size':'large'})
                axins.plot(Data_to_plot.index.get_level_values(level=0), Data_to_plot[Ratio],'o--', color = colors4[i], linewidth=0.8, label= Doses[i])
                axins.errorbar(Data_to_plot.index.get_level_values(level=0), Data_to_plot[Ratio], yerr= Data_to_plot[Error], fmt='none', elinewidth = 0.8, ecolor=colors4[i],capthick=1, capsize=2, label ='')
                axins.set_ylim([0.8,3.1])
                axins.set_xticks([0.5,24,48,72])
                axins.grid(which='major',linestyle='--', linewidth='0.5', color='lightgray')
                plt.xticks([0.5,24,48,72], fontsize=13)
                plt.yticks(fontsize=13)
                plt.grid(which='major',linestyle='--', linewidth='0.5', color='lightgray')
        plt.savefig(Filepath+'H2AX_'+Phase+'_intensity.pdf')
        plt.show()

def Plotting_average_proton(AverageData,Filepath, Phase, title):
        Ratio = 'Ratio '+ Phase
        Error = 'Error '+ Phase 
        positions = AverageData.index.get_level_values(0).unique()
        timepoints = AverageData.index.get_level_values(1).unique()
        dose  = AverageData.index.get_level_values(2).unique()
        number_of_timepoints = len(timepoints.values)
   
        #Remove dose from index
        AverageData.index = AverageData.index.droplevel(2)


        #Create figure with subplots matching the number of doses
        fig,ax = plt.subplots(figsize=(6,5))
        box = ax.get_position()
        ax.set_position([box.x0,box.y0,box.width*0.9,box.height])

        for i in range(number_of_timepoints):
                positions_to_plot = AverageData.loc[AverageData.index.get_level_values(level=1)==timepoints[i]]
                positions_to_plot.index = positions_to_plot.index.droplevel(1)

                plt.plot(positions_to_plot.index.get_level_values(level=0), positions_to_plot[Ratio],'o--', color = colors2[i], linewidth=0.8, label= timepoints[i])
                plt.errorbar(positions_to_plot.index.get_level_values(level=0), positions_to_plot[Ratio], yerr= positions_to_plot[Error], fmt='none', elinewidth = 0.8, ecolor=colors2[i],capthick=1, capsize=2, label ='')
                plt.title(title, fontsize = 16)
                plt.ylabel(r'Relative $\gamma$-H2AX fluorescence intensity', fontsize = 14)
                plt.xlabel('Position', fontsize = 14)
                plt.xticks([1,2,3,4,5])
                plt.grid(which='major',linestyle='--', linewidth='0.5', color='lightgray')
                leg = plt.legend(bbox_to_anchor = (1.24,0.62),fontsize = 12)
                leg.set_title('Dose [Gy]', prop = {'size':'large'})
        plt.savefig(Filepath+'H2AX_'+ Phase+'_intensity.pdf')
        plt.show()


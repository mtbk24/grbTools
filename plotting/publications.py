from __future__ import division
import matplotlib.pyplot as plt
import numpy as np


def journal():
    '''
    Appears nearly square. I like this best.
    'figure.figsize': [3.1, 2.6]
    
    Use: plt.tight_layout(pad=0.1, w_pad=0.0, h_pad=0.0)
    plt.xlabel('$E_{iso}$ $(erg)$',labelpad=-1)  
    plt.ylabel('$E^*_{pk}$ $(keV)$',labelpad=-2)
    
    '''
    params = {'backend': 'pdf',
              'axes.labelsize':  10,
              'font.size':       10,
              'legend.fontsize': 8,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'xtick.direction': 'in',
              'ytick.direction': 'in',
              'text.usetex':     True,
              'figure.figsize': [3.1, 2.6],
              'font.family': 'serif',}
    plt.rcParams.update(params)


def journal2():
    '''
    Wider than it is tall.
    'figure.figsize': [3.4, 2.45]
    
    Use: plt.tight_layout(pad=0.1, w_pad=0.0, h_pad=0.0)
    plt.xlabel('$E_{iso}$ $(erg)$',labelpad=-1)  
    plt.ylabel('$E^*_{pk}$ $(keV)$',labelpad=-2)
    
    '''
    params = {'backend': 'pdf',
              'axes.labelsize':  10,
              'font.size':       10,
              'legend.fontsize': 8,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'xtick.direction': 'in',
              'ytick.direction': 'in',
              'text.usetex':     True,
              'figure.figsize': [3.4, 2.45],
              'font.family': 'serif',}
    plt.rcParams.update(params)


def set_journal(fig_width=245.26653, height_factor=1.):
    """
    Sets the matplotlib rc params to generate plots that are the proper
    size for two column journals. Simple specify the column width
    :param fig_width: Get this from LaTeX using \showthe\columnwidth
    """
    

    fig_width_pt =fig_width 
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean * height_factor      # height in inches
    fig_size =  [fig_width,fig_height]
    params = {'backend': 'pdf',
              'axes.labelsize': 10,
              'font.size': 10,
              'legend.fontsize': 10,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'xtick.direction': 'in',
              'ytick.direction': 'in',
              'text.usetex': True,
              'figure.figsize': fig_size,
              'font.family': 'serif'}
    plt.rcParams.update(params)


def set_journal_2(fig_width=245.26653, height_factor=1.):
    """
    Sets the matplotlib rc params to generate plots that are the proper
    size for two column journals. Simple specify the column width
    :param fig_width: Get this from LaTeX using \showthe\columnwidth
    """
    

#    fig_width_pt = fig_width
#    inches_per_pt = 1.0/72.27               # Convert pt to inch
#    golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
#    fig_width = fig_width_pt*inches_per_pt  # width in inches
#    fig_height = fig_width*golden_mean * height_factor      # height in inches
    fig_size =  [3.39375300954753, 2.45]
    params = {'backend': 'pdf',
              'axes.labelsize': 10,
              'font.size': 10,
              'legend.fontsize': 10,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'xtick.direction': 'in',
              'ytick.direction': 'in',
              'text.usetex': True,
              'figure.figsize': fig_size,
              'font.family': 'serif',}
    plt.rcParams.update(params)






def set_journal_3(fig_width=245.26653, height_factor=1.):
    """
    Sets the matplotlib rc params to generate plots that are the proper
    size for two column journals. Simple specify the column width
    :param fig_width: Get this from LaTeX using \showthe\columnwidth
    """
    fig_size =  [3.4, 2.45]
    params = {'backend': 'pdf',
              'axes.labelsize': 10,
              'font.size': 10,
              'legend.fontsize': 10,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'xtick.direction': 'in',
              'ytick.direction': 'in',
              'text.usetex': True,
              'figure.figsize': fig_size,
              'font.family': 'serif',}
    plt.rcParams.update(params)






def pltreset():
    """
    Reset the plotting settings
    """
    
    plt.rcdefaults()
    #optional for notebook
    #%matplotlib inline
    #%matplotlib notebook




def print_vertical_subplots():
    toprint = '''
        Pars1 = {'x':           d1.eiso,
                 'y':           d1.epeak,
                 'xerr':        [d1.eiso_err.tolist(), d1.eiso_err.tolist()],
                 'yerr':        [d1.epeak_err_low.tolist(), d1.epeak_err_up.tolist()],
                 'fmt':          '.', 
                  'ms':           4,         
                 'color':        'darkslategrey',
                 'capsize':      0,
                 'mew':          0,               
                 'mec':          'darkslategrey',    
                 'lw':           1,                
                 'alpha':        0.65,
                 'label':        'Amati (38)',
                    }
                    
        Pars2 = {'x':           d2.eiso,
                 'y':           d2.epeak,
                 'xerr':        [d2.eiso_err.tolist(), d2.eiso_err.tolist()],
                 'yerr':        [d2.epeak_err_low.tolist(), d2.epeak_err_up.tolist()],
                 'fmt':          '.',          
                 'ms':           8,            
                 'color':        'blue',
                 'capsize':      0,
                 'mew':          0.25,                  
                 'mec':          'k',    
                 'lw':           1,                 
                 'alpha':        1,
                 'label':        'Ours, G+L (8)',
                    }


        # FOR COEFFICIENTS
        clrs   = ['k','k','dimgrey','darkslategrey','blue',
                  'firebrick','mediumturquoise']
                  
        locs   = [[0.72, 0.50], 
                  [0.72, 0.48],
                  [0.72, 0.43],
                  [0.72, 0.38],
                  [0.72, 0.33],
                  [0.72, 0.28],
                  [0.72, 0.23]]  
                  
        props  = dict(boxstyle='round', facecolor='white', alpha=0.0, lw=0)



        yLims       = [20, 0.6E4]
        xLims       = [3.0E51, 0.5E55]

        # CLEAR ALL PLOTS
        plt.clf()
        plt.close()
        #
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True,
                                       figsize=(3.5, 5.5))
        # SETTINGS FOR AMATI RELATION
        trendArgs   = dict(color='grey', alpha=0.7, lw=2)
        fillArgs    = dict(color='grey', alpha=0.08)
        for ax in (ax1, ax2):
            ax.grid(which='both', axis='both', alpha=0.1)
            PLT_Amati_Relation_2013_Fill(ax = ax,
                                        trendArgs = trendArgs,
                                        fillArgs  = fillArgs)
            ax.minorticks_on()    
            ax.loglog()
            ax.set_xlim(xLims[0], xLims[1])
            ax.set_ylim(yLims[0], yLims[1])
            ax.set_ylabel('$E^{*}_{pk}$ ($keV$)', labelpad=-2)


        ax1.text(0.35, 0.10, 'Full overlap', transform=ax1.transAxes, fontsize=12, fontstyle='italic')
        ax1.errorbar(**Pars1)
        ax1.errorbar(**Pars2)
        ax1.errorbar(**Pars3)

        ax2.text(0.35, 0.10, 'Constrained only', transform=ax2.transAxes, fontsize=12, fontstyle='italic')
        ax2.errorbar(**Pars4)
        ax2.errorbar(**Pars5)
        ax2.errorbar(**Pars6)


        for X,Y,DX,DY in zip(data_FULL.eiso_MINE, data_FULL.epeak_MINE, data_FULL.eiso_AMATI, data_FULL.epeak_AMATI):
            ax1.annotate("", xy=(DX, DY), xytext=(X, Y), arrowprops=dict(arrowstyle="-", color='grey', alpha=0.6)) 

        for X,Y,DX,DY in zip(data_GOOD.eiso_MINE, data_GOOD.epeak_MINE, data_GOOD.eiso_AMATI, data_GOOD.epeak_AMATI):
            ax2.annotate("", xy=(DX, DY), xytext=(X, Y), arrowprops=dict(arrowstyle="-", color='grey', alpha=0.6)) 
            
        # ------
        # Full Sample        

        textstr = [' \ \ \ \ \  b  \ \ \ \ \ \ \ m ',
                   '----- \ \ \ \ -----',
                    '   1.958 \ \ 0.549', 
                    '   2.071 \ \ 0.495', 
                    '   2.372 \ \ 0.398', 
                    '   2.069 \ \ 0.511', 
                    '   2.096 \ \ 0.507'] 
                   
        for ii in range(0, len(locs)):
            ax1.text(locs[ii][0], locs[ii][1], 
                     textstr[ii], 
                     transform         = ax1.transAxes, 
                     fontsize          = 10, 
                     color             = clrs[ii], 
                     verticalalignment = 'top', 
                     bbox              = props)

        # ------
        # Good only
        textstr = [' \ \ \ \ \  b  \ \ \ \ \ \ \ m ',
                   '----- \ \ \ \ -----',
                   '   1.958 \ \ 0.549', 
                   '   2.029 \ \ 0.497', 
                   '   2.372 \ \ 0.398', 
                   '   1.913 \ \ 0.533', 
                   '   1.964 \ \ 0.546']
                   
        for ii in range(0, len(locs)):
            ax2.text(locs[ii][0], locs[ii][1], 
                     textstr[ii], 
                     transform         = ax2.transAxes, 
                     fontsize          = 10, 
                     color             = clrs[ii], 
                     verticalalignment = 'top', 
                     bbox              = props)

        for ax in (ax1, ax2):
            ax.legend(loc=2, numpoints=1, 
                      handletextpad=0, 
                      labelspacing=0.1,
                      frameon = False)


        ax2.set_xlabel('$E_{iso}$ ($erg$)', labelpad=-1)
        fig.tight_layout(pad=0.1, w_pad=0.0, h_pad=0.0)
        fig.savefig('X.pdf')
        '''
    return toprint



def print_horizontal_subplots():
    toprint = '''
            Pars1 = {'x':           d1.eiso,
                 'y':           d1.epeak,
                 'xerr':        [d1.eiso_err.tolist(), d1.eiso_err.tolist()],
                 'yerr':        [d1.epeak_err_low.tolist(), d1.epeak_err_up.tolist()],
                 'fmt':          '.', 
                  'ms':           4,         
                 'color':        'darkslategrey',
                 'capsize':      0,
                 'mew':          0,               
                 'mec':          'darkslategrey',    
                 'lw':           1,                
                 'alpha':        0.65,
                 'label':        'Amati (38)',
                    }
                    
        Pars2 = {'x':           d2.eiso,
                 'y':           d2.epeak,
                 'xerr':        [d2.eiso_err.tolist(), d2.eiso_err.tolist()],
                 'yerr':        [d2.epeak_err_low.tolist(), d2.epeak_err_up.tolist()],
                 'fmt':          '.',          
                 'ms':           8,            
                 'color':        'blue',
                 'capsize':      0,
                 'mew':          0.25,                  
                 'mec':          'k',    
                 'lw':           1,                 
                 'alpha':        1,
                 'label':        'Ours, G+L (8)',
                    }
                    


        trendArgs   = dict(color='grey', alpha=0.7, lw=2)
        fillArgs    = dict(color='grey', alpha=0.08)

        # FOR COEFFICIENTS
        clrs   = ['k','k','dimgrey','darkslategrey','blue',
                  'firebrick','mediumturquoise']          
        locs   = [[0.7, 0.50], 
                  [0.7, 0.48],
                  [0.7, 0.43],
                  [0.7, 0.38],
                  [0.7, 0.33],
                  [0.7, 0.28],
                  [0.7, 0.23]]           
        props  = dict(boxstyle='round', facecolor='white', alpha=0.0, lw=0)
        
        
        plt.clf()
        plt.close()
        
        yLims       = [20, .6E4]
        xLims       = [3.0E51, .5E55]
        
        # BEGIN FIGURE
        fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True,
                                       figsize=(5.8, 2.8)) 
        for ax in (ax1, ax2):
            #ax.grid(which='both', axis='both', alpha=0.1)
            PLT_Amati_Relation_2013_Fill(ax = ax,
                                    trendArgs = trendArgs,
                                    fillArgs = fillArgs)  
            ax.minorticks_on()    
            ax.loglog()
            ax.set_xlim(xLims[0], xLims[1])
            ax.set_ylim(yLims[0], yLims[1])


        # line(x, m, b)
        Xdata = np.linspace(48, 55, 100)

        Ydata = line(Xdata-52.0, 0.495, 2.071) 
        ax1.plot(10**Xdata, 10**Ydata, 
                 color='darkslategrey', alpha=0.7, ls='--', lw=1.5)

        Ydata = line(Xdata-52.0, 0.398, 2.372) 
        ax1.plot(10**Xdata, 10**Ydata, 
                 color='blue', alpha=0.8, ls='--', lw=1.5)
                 
        Ydata = line(Xdata-52.0, 0.511, 2.069) 
        ax1.plot(10**Xdata, 10**Ydata, 
                 color='firebrick', alpha=0.8, ls='--', lw=1.5)

        Ydata = line(Xdata-52.0, 0.507, 2.096) 
        ax1.plot(10**Xdata, 10**Ydata, 
                 color='cyan', alpha=0.8, ls='--', lw=1.5)   
                 

        Ydata = line(Xdata-52.0, 0.497, 2.029) 
        ax2.plot(10**Xdata, 10**Ydata, 
                 color='darkslategrey', alpha=0.7, ls='--', lw=1.5)

        Ydata = line(Xdata-52.0, 0.398, 2.372) 
        ax2.plot(10**Xdata, 10**Ydata, 
                 color='blue', alpha=0.8, ls='--', lw=1.5)
                 
        Ydata = line(Xdata-52.0, 0.533, 1.913) 
        ax2.plot(10**Xdata, 10**Ydata, 
                 color='firebrick', alpha=0.8, ls='--', lw=1.5)

        Ydata = line(Xdata-52.0, 0.546, 1.964) 
        ax2.plot(10**Xdata, 10**Ydata, 
                 color='cyan', alpha=0.8, ls='--', lw=1.5)   
                 

        ax1.text(0.35, 0.10, 'Full overlap', transform=ax1.transAxes, fontsize=12, fontstyle='italic')
        ax1.errorbar(**Pars1)
        ax1.errorbar(**Pars2)
        ax1.errorbar(**Pars3)

        ax2.text(0.35, 0.10, 'Constrained only', transform=ax2.transAxes, fontsize=12, fontstyle='italic')
        ax2.errorbar(**Pars4)
        ax2.errorbar(**Pars5)
        ax2.errorbar(**Pars6)

        for X,Y,DX,DY in zip(data_FULL.eiso_MINE, data_FULL.epeak_MINE, data_FULL.eiso_AMATI, data_FULL.epeak_AMATI):
            ax1.annotate("", xy=(DX, DY), xytext=(X, Y), arrowprops=dict(arrowstyle="-", color='grey', alpha=0.6)) 

        for X,Y,DX,DY in zip(data_GOOD.eiso_MINE, data_GOOD.epeak_MINE, data_GOOD.eiso_AMATI, data_GOOD.epeak_AMATI):
            ax2.annotate("", xy=(DX, DY), xytext=(X, Y), arrowprops=dict(arrowstyle="-", color='grey', alpha=0.6)) 
            
        # ------
        # Full Sample        

        textstr = [' \ \ \ \ \  b  \ \ \ \ \ \ \ m ',
                   '----- \ \ \ \ -----',
                    '   1.958 \ \ 0.549', 
                    '   2.071 \ \ 0.495', 
                    '   2.372 \ \ 0.398', 
                    '   2.069 \ \ 0.511', 
                    '   2.096 \ \ 0.507'] 
                   
        for ii in range(0, len(locs)):
            ax1.text(locs[ii][0], locs[ii][1], 
                     textstr[ii], 
                     transform         = ax1.transAxes, 
                     fontsize          = 10, 
                     color             = clrs[ii], 
                     verticalalignment = 'top', 
                     bbox              = props)

        # ------
        # Good only
        textstr = [' \ \ \ \ \  b  \ \ \ \ \ \ \ m ',
                   '----- \ \ \ \ -----',
                   '   1.958 \ \ 0.549', 
                   '   2.029 \ \ 0.497', 
                   '   2.372 \ \ 0.398', 
                   '   1.913 \ \ 0.533', 
                   '   1.964 \ \ 0.546']
                   
        for ii in range(0, len(locs)):
            ax2.text(locs[ii][0], locs[ii][1], 
                     textstr[ii], 
                     transform         = ax2.transAxes, 
                     fontsize          = 10, 
                     color             = clrs[ii], 
                     verticalalignment = 'top', 
                     bbox              = props)

        for ax in (ax1, ax2):
            ax.legend(loc=2, numpoints=1, 
                      handletextpad=0, 
                      labelspacing=0.1,
                      frameon = False)
            ax.set_xlabel('$E_{iso}$ ($erg$)', labelpad=-1)


        ax2.set_yticks([]) # turn off yaxis tick labels
        ax1.set_ylabel('$E^{*}_{pk}$ ($keV$)', labelpad=-2)
        fig.tight_layout(pad=0.1, w_pad=0, h_pad=0)
        fig.subplots_adjust(hspace=0, wspace=0)
        fig.savefig('X.pdf')

        '''
    return toprint
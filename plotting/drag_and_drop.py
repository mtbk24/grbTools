'''
http://scipy-cookbook.readthedocs.io/items/Matplotlib_Drag_n_Drop_Text_Example.html

Event Handling
https://matplotlib.org/users/event_handling.html

'''

from __future__ import division
import numpy as np
import pandas as pd
from Amati.tools import *
from matplotlib import pylab as p
from matplotlib.text import Text


class DragHandler(object):
    """ 
    A simple class to handle Drag n Drop.
    This is a simple example, which works for Text objects only.
    Good for moving labels.  When using the ipython gui, make 
    tight layouts within that instead of the code here.
    Save figure as .pdf or .eps for lossless resolution.
    .png .jpg and .tiff image files suffer resolution loss and
    no settings of dpi or figsize can help.


    Muse use imports: 
    from matplotlib import pylab as p
    from matplotlib.text import Text

    EXAMPLE: 
    --------------------
	p.clf()
	p.close()
	yLims       = [2E1, 6E3]  
	xLims       = [1.5E51, 6E54]
	fig, (ax1, ax2) = p.subplots(1, 2, sharey=True, dpi=100, figsize=(10,5)) 
	for ax in [ax1, ax2]:
	    ax.grid(which='both', axis='both', alpha=0.2)
	    PLT_Amati_Relation_2013_2(color='black', alpha=0.5, ax=ax)
	    ax.hlines([8, 900], xLims[0], xLims[1], color='purple', lw=2, alpha=0.1)                # NaI Boundaries
	    ax.hlines([200, 38000], xLims[0], xLims[1], color='darkgoldenrod', lw=2, alpha=0.1)     # BGO Boundaries
	    ax.fill_between([xLims[0], xLims[1]],8, 900, color='purple', alpha=0.1)                 # NAI
	    ax.fill_between([xLims[0], xLims[1]], 200, 38000, color='darkgoldenrod', alpha=0.1)     # BGO
	    ax.set_xlim(xLims[0], xLims[1])
	    ax.set_ylim(yLims[0], yLims[1]) 
	    ax.loglog()
	    ax.set_xlabel('$E_{iso}$ $(erg)$', fontsize=16, labelpad=-3)

	    x           = df_LAT_X.eiso.tolist()
	    xerr        = [df_LAT_X.eiso_err_low.tolist(), 
	                   df_LAT_X.eiso_err_up.tolist()]
	    y           = df_LAT_X.epeak.tolist()
	    yerr        = [df_LAT_X.epeak_err_low.tolist(), 
	                   df_LAT_X.epeak_err_up.tolist()] 
	    labels      = df_LAT_R.number.tolist()
	    
	    ax.errorbar(x, y, xerr=xerr, yerr=yerr,     
	                   fmt          = '.',
	                   ms           = 11,
	                   color        = 'blue',
	                   mec          = 'k',
	                   capsize      = 0,
	                   alpha        = 0.7,
	                   label        = 'GBM+LAT') 
	                   
	    # add labels and set their picker attribute to True
	    for a,b,l in zip(x, y, labels):
	        ax.text(a, b, l, picker=True)

	    # Create the event hendler 
	    dragh = DragHandler()    
	               
	ax1.set_ylabel(r'$\nu F_{\nu}$  $peak$ $energy$ $(keV)$', fontsize=18, labelpad=-5)   
	ax2.legend(loc=4, numpoints=1, fontsize=16, handletextpad=0, labelspacing=0)
	p.show()


    
    """
    def __init__(self, figure=None) :
        """ Create a new drag handler and connect it to the figure's event system.
        If the figure handler is not given, the current figure is used instead
        """

        if figure is None : figure = p.gcf()
        # simple attibute to store the dragged text object
        self.dragged = None

        # Connect events and callbacks
        figure.canvas.mpl_connect("pick_event", self.on_pick_event)
        figure.canvas.mpl_connect("button_release_event", self.on_release_event)

    def on_pick_event(self, event):
        " Store which text object was picked and were the pick event occurs."

        if isinstance(event.artist, Text):
            self.dragged = event.artist
            self.pick_pos = (event.mouseevent.xdata, event.mouseevent.ydata)
        return True

    def on_release_event(self, event):
        " Update text position and redraw"

        if self.dragged is not None :
            old_pos = self.dragged.get_position()
            new_pos = (old_pos[0] + event.xdata - self.pick_pos[0],
                       old_pos[1] + event.ydata - self.pick_pos[1])
            self.dragged.set_position(new_pos)
            self.dragged = None
            p.draw()
        return True



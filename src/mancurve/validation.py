import os
import pandas as pd
import numpy as np
import core as mc
import datetime as dt
import h5py
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter, DayLocator, HourLocator
import time    

# Plot 
pd.plotting.register_matplotlib_converters()

def import_data(file):
    '''HDF5 Datei importieren und Plots erzeugen. Es werden alle Initialisierungs-
    Zeiten als eine Kurve dargestellt. Daneben wird die Beobachtung aufgetragen.
    
    :param file: Dateiname
    :type file: String
    '''
    f = h5py.File(file, 'r')
    init  = list(f.keys())
    dif   = {}
    n = len(init)
    for idx,init in enumerate(init):
        data  = pd.read_hdf(file, init)
        cs = np.unique([x[:8] for x in data.columns])        
        obj = mc.combine_curves()
    
        if idx == 0:
            for p in cs:
                dif[p] = np.empty(0,dtype=int)

        for pegel in cs:     
            obj.stations = [pegel]
            tp = data[pegel+'_MAN'].dropna().index[-1]
            obj.get_obs(download=False,zeitpunkt=tp+dt.timedelta(minutes=30))
  
            obj.po[pegel].data = obj.po[pegel].data.loc[data.index[0]:]#[::5]
            data[pegel+'_MC'] = data[pegel+'_MC'].loc[:tp]#[::5]
            
            d   = (data[pegel+'_MC'] -
                   obj.po[pegel].data['waterlevel']).dropna()*100
            dif[pegel] = np.append(dif[pegel],d.to_numpy())
                        
            fig = plt.figure(pegel)    
            lw = 2
            
            plt.plot(obj.po[pegel].data.index[::30],
                         obj.po[pegel].data['waterlevel'][::30],
                         color='k',ms = 3,#marker = 'o',
                         lw=lw, label=None,zorder=100)
            plt.plot(data[pegel+'_MC'].index,data[pegel+'_MC'],
                     color=[0.8,1-(idx/n),idx/n])
                          
            line,caps,bars=plt.errorbar(data[pegel+'_MAN'].dropna().index,     # X
                            data[pegel+'_MAN'].dropna(),    # Y
                            yerr= data[pegel+'_MANUNC'].dropna(),# Y-errors
                            color=[.8,1-(idx/n),idx/n],   # format line like for plot()
                            marker = 's',
                            linewidth=0,   # width of plot line
                            elinewidth=2,  # width of error bar line
                            ecolor='grey',    # color of error bar
                            capsize=2,     # cap length for error bar
                            capthick=1,   # cap thickness for error bar
                            label='Manuell '+ init,
                            ms = 4,
                            zorder = 120                        
                            )
            ax = plt.gca()
            
            # Format Time axis -> Hours on top, Days at bottom of plot          
            ax.xaxis.set_major_locator(DayLocator())   
            ax.xaxis.set_major_formatter(DateFormatter('%a %d.%m.%y'))                        
            ax.xaxis.set_minor_locator(HourLocator(byhour=[0,6,12,18]))    
            ax.xaxis.set_minor_formatter(DateFormatter('%H Uhr (UTC)'))    
            ax.grid(which='major', axis='x',)
        
            # rotates and right aligns the x labels, and moves the bottom of the
            # axes up to make room for them
            fig.autofmt_xdate()
            
            #plt.xlabel('Minutes')
            plt.ylabel('Waterlevel')
            plt.title('Manuelle Kurve: ' + pegel)
            plt.legend(loc="upper right", scatterpoints=1, prop={'size': 6})
            #plt.savefig(os.path.join(self.dirname, 'figs/'+key+'.png'),dpi=96)
    plt.show()    
    plt.close('all')        
    return dif


        
def main(data):
    ''' Hauptfunktion. Hier werden alle Intialisierungszeiten in einem Array
    gesammelt, um sie statistisch auszuwerten. Dies könnte auf alle Daten 
    erweitert werden.
    
    :param data: file
    :type data:  string
    '''
    all_data = None
    if type(data) != str:
        for d in data:
            dic = import_data(d)
            if not all_data:
                all_data = dic
            else:
                for k in all_data:
                    all_data[k] = np.append(all_data[k],dic[k])
    else:
        all_data = import_data(data)
        
    for k in all_data:
        x = pd.DataFrame(all_data[k])
        
        if k == 'pgl_506P':
            print(k,'\n',x.describe(percentiles=[.01,.05,.95,.99]))    
        #plt.figure()
        #plt.hist(all_data[k], bins = np.arange(-100,110,5), 
        #         alpha = .8, label = 'Manuelle Kurve',
        #         weights=np.ones(len(all_data[k])) / len(all_data[k]),
        #         color = 'r',edgecolor='grey', zorder = 1) 
    #plt.show()
    
if __name__ == '__main__':
    path = os.path.join(os.path.abspath('../../data/data/'),'')
    file = '20171029.h5'
    file = '20200210.h5'
    data = path+file   
    #data = [path+'20200303.h5',path+'20200304.h5',path+'20200305.h5',
    #        path+'20200306.h5',path+'20200307.h5',]
    stime1 = time.time()    
    main(data)    
    print("Total time elapsed: %.1f s" % (time.time() - stime1))
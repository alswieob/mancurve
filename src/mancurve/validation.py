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

def import_data(file,tp):
    f = h5py.File(file, 'r')
    init  = list(f.keys())
    dif   = {}
    figs  = {} 
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
            
            #dif2 = data[pegel+'_MAN'] - obj.po[pegel].data['waterlevel']
                        
            fig = plt.figure(pegel)    
            lw = 2
            
            plt.plot(obj.po[pegel].data.index[::30],
                         obj.po[pegel].data['waterlevel'][::30],
                         color='k',ms = 3,#marker = 'o',
                         lw=lw-1, label=None,zorder=100)
            plt.plot(data[pegel+'_MC'].index,data[pegel+'_MC'],
                     color=[0.8,1-(idx/n),idx/n])
            #plt.fill_between(data[pegel+'_GPR'].index,
            #                 data[pegel+'_GPR'] - data[pegel+'_GPRSTD'],
            #                 data[pegel+'_GPR'] + data[pegel+'_GPRSTD'],
            #                 color='darkorange',alpha=0.2, zorder=5)
    
            #plt.plot(data[pegel+'_MAN'].dropna().index,
            #         data[pegel+'_MAN'].dropna(),
            #         'x')
                          
            line,caps,bars=plt.errorbar(data[pegel+'_MAN'].dropna().index,     # X
                            data[pegel+'_MAN'].dropna(),    # Y
                            yerr= data[pegel+'_MANUNC'].dropna(),# Y-errors
                            color=[.8,1-(idx/n),idx/n],#"#1f77b4",    # format line like for plot()
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
    #plt.show()    
    plt.close('all')        
    return dif


        
def main(data,zeitpunkt):
    all_data = None
    if type(data) != str:
        for d in data:
            dic = import_data(d,tp = zeitpunkt)
            if not all_data:
                all_data = dic
            else:
                for k in all_data:
                    all_data[k] = np.append(all_data[k],dic[k])
    else:
        all_data = import_data(data,tp = zeitpunkt)
        
    for k in all_data:
        x = pd.DataFrame(all_data[k])
        print(k,'\n',x.describe(percentiles=[.01,.05,.95,.99]))    
        #plt.figure()
        #plt.hist(all_data[k], bins = np.arange(-100,110,5), 
        #         alpha = .8, label = 'Manuelle Kurve',
        #         weights=np.ones(len(all_data[k])) / len(all_data[k]),
        #         color = 'r',edgecolor='grey', zorder = 1) 
    #plt.show()
    
if __name__ == '__main__':
    path = os.path.join(os.path.abspath('../../data/data/'),'')
    file = '20200320.h5'
    data = path+file   
    #data = [path+'20200303.h5',path+'20200304.h5',path+'20200305.h5',
    #        path+'20200306.h5',path+'20200307.h5',]
    zeitpunkt = dt.datetime(2020,3,13,1,00)
    stime1 = time.time()    
    main(data, zeitpunkt)    
    print("Total time elapsed: %.1f s" % (time.time() - stime1))
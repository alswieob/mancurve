# -*- coding: utf-8 -*-
"""
This script runs the operational manual waterlevel forecast curve construction.
It can be executed from the command line.
"""
import logging
import subprocess
import os
import numpy as np
import datetime as dt
from pytz import timezone
import re
import matplotlib.pyplot as plt
import pandas as pd
from itertools import compress
import shutil, glob
import time      
from matplotlib.dates import DateFormatter, DayLocator, HourLocator
import pegelonline_getter as obs
import data_feeder as df
from scipy.signal import savgol_filter

__author__ = "Luis Becker"
__copyright__ = "Luis Becker"
__license__ = "mit"

_logger = logging.getLogger(__name__)

class Object(object):
    pass

class combine_curves():
    def __init__(self):
        self.FNULL = open(os.devnull, 'w')
        dirname = os.path.dirname(__file__)
        self.aim_dir = os.path.abspath('../../data/temp/')
        
    def clear_old_data(self):
        path  = os.path.abspath('../../data/temp/')
        path  = os.path.join(path, '')
        
        files = os.listdir(path)
        for f in files:
            if f.endswith('.dat') or f.endswith('.txt'):
                x = path + f
                os.remove(x) 
                
    def summarize_figures(self):
        import PIL
        path    = os.path.abspath('../../data/figs/')
        list_im = os.listdir(path) 
        list_im = [os.path.abspath('../../data/figs/' + x) for x in list_im]        
        imgs    = [ PIL.Image.open(i) for i in list_im ]
        n = len(list_im)
        # pick the image which is the smallest, 
        # and resize the others to match it (can be arbitrary image shape here)
        min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
        
        vert = []
        for i in np.arange(0,n,2):
            ims = imgs[i:i+2]
            imgs_comb = np.hstack( [np.asarray( i.resize(min_shape) ) for i in ims ] )
        
            # save that beautiful picture
            imgs_comb = PIL.Image.fromarray( imgs_comb)
            outfile   = os.path.abspath('../../data/output/combined'+str(i)+'.png')
            vert.append(outfile)
            imgs_comb.save(outfile)    
        
        # for a vertical stacking it is simple: use vstack
        imgs    = [ PIL.Image.open(i) for i in vert]
        imgs_comb = np.vstack([np.asarray(i) for i in imgs])
        imgs_comb = PIL.Image.fromarray( imgs_comb)
        imgs_comb.save(os.path.abspath('../../data/output/combined.png'))
        
        for x in vert:
            os.remove(x)
          
        # Rename file for log
        f1 = os.path.abspath('../../data/output/combined.png')
        f2 = os.path.abspath('../../data/output/'
                          +dt.datetime.now().strftime('%Y%m%d%H%MZ')+'.png')
        os.rename(f1,f2)
    
    def summarize_data(self):
        x = glob.glob(os.path.abspath('../../data/output/*.h5'))
        today = dt.datetime.now().strftime('%Y%m%d')
        new_file = './data/'+today+'.h5'
        new_f = pd.DataFrame()
        
        for i in x:
            new_f = pd.concat([new_f,
                               pd.read_hdf(i,'CombinedForecastUTC')], axis=1)
        
        # Store Data in file and create frame with init information
        new_f.to_hdf(new_file,mode='a',complevel = 1, 
                     key=dt.datetime.now().strftime('T%Y%m%d%H%M'))
           

    def get_current_NWHW(self, download = False,
                         zeitpunkt = dt.datetime(2020,2,1,12,16)):       
        if download:
            # NWHW
            cmd = 'ssh bm11mos@linwvd10 ls /FEWS/HWNW_Vorhersage/HWNW* | tail -1'    
            p = subprocess.Popen(cmd,shell=True,stdout = subprocess.PIPE
                                 ,stderr = subprocess.PIPE)
            p.wait()
            filename = p.stdout.read().decode()
            p.stdout.close()
            p.stderr.close()
            
            cmd = 'bm11mos@linwvd10:'+filename
            p = subprocess.Popen(["scp","-T",cmd,self.aim_dir],
                                 stdout = self.FNULL,
                                 stderr = subprocess.STDOUT)
            p.wait()
            aim_file =  os.path.join(self.aim_dir, filename.rsplit('/', 1)[-1].strip())
            
        else:            
            x = df.nwhw(zeitpunkt)
            aim_file = x.valid_file
            
        self.filename = aim_file
   
    def get_eformat(self, download = False):
        self.dirname = os.path.dirname(__file__)
        self.aim_dir = os.path.abspath('../../data/temp/')
        # Add trailing slash
        self.aim_dir = os.path.join(self.aim_dir, '') 
        
        self.MHW = {}
        self.MNW = {}
        for pegel in self.stations:
            # Immer dieses und nächstes Jahr berücksichtigen
            year = dt.datetime.now().strftime('%Y')
            
            # EFORMAT: Mittlere Hoch-Niedrigwasser
            if download:
                cmd = ("ssh bm14gast@lingezeiten11 ls /data/vb_hwnw/deu"+year+"/"
                        "*"+pegel[-4:]+"*")  
                p = subprocess.Popen(cmd,shell=True,stdout = subprocess.PIPE
                                    ,stderr = subprocess.PIPE)
                filename = p.stdout.read().decode()
                p.stdout.close()
                p.stderr.close()    
                p.wait()
                
                cmd = 'bm14gast@lingezeiten11:'+filename
                p = subprocess.Popen(["scp","-T",cmd,self.aim_dir],
                                     stdout = self.FNULL,
                                     stderr = subprocess.STDOUT)
                p.wait()
                filename = self.aim_dir + filename.rsplit('/', 1)[-1].strip()
                self.read_NWHW_file(filename,pegel)  
            else:
                break
            
        if not download:
            # Get Special Event / year / archive
            files = os.listdir(os.path.abspath('../../data/arch/EFORMAT'))
            pnames = ['pgl' + x[-13:-8] for x in files]
            files  = [os.path.join(os.path.abspath('../../data/arch/EFORMAT'),x) for x in files]           
            for idx,source_file in enumerate(files):            
                shutil.copy2(source_file,self.aim_dir)
                self.read_NWHW_file(source_file, pnames[idx])
        

    def read_NWHW_file(self,filename,pegel):
        with open(filename, mode = 'rb') as f:
            for idx,line in enumerate(f):
                if idx == 21:
                    self.MHW[pegel] = \
                            float(re.findall(r'\d.\d\d',str(line))[0])
                elif idx == 22:
                    self.MNW[pegel] = \
                            float(re.findall(r'\d.\d\d',str(line))[0])
                elif idx > 22:
                    f.close()
                    break         

    def get_mos(self, download = False,zeitpunkt = dt.datetime(2020,2,1,12,16)):
        self.mos_data = {}
        #self.ast_data = {}        
        self.mos_errstate = {}
        
        if download:
            for pegel in self.stations:
            
                # MOS
                cmd = 'ssh bm11mos@linwvd10 ls /data_MOS/bm11mos/1704/output/*' +\
                       pegel+'* | tail -1'
                p = subprocess.Popen(cmd,shell=True,stdout = subprocess.PIPE
                                     ,stderr = subprocess.PIPE)
                filename = p.stdout.read().decode()
                p.stdout.close()
                p.stderr.close()
                p.wait()
                
                cmd = 'bm11mos@linwvd10:'+filename
                p = subprocess.Popen(["scp","-T",cmd,self.aim_dir],
                                     stdout = self.FNULL,
                                     stderr = subprocess.STDOUT)
                p.wait()

                filename = self.aim_dir+filename.split('/')[-1].strip()
                
                self.mos_filename = filename
                self.pegel_name = pegel
                self.load_mos()  
                   
        else:
            x = df.mos(zeitpunkt = zeitpunkt, freq=15)
            x.main()
            source_file = ("/home/bm1483/Development/manuelle_kurve"
                           "/scenarios/mos_data.h5")          
            data = pd.read_hdf(source_file,'table')
            for pegel in self.stations:
                self.mos_data[pegel] = data[pegel]    
            self.mos_errstate = x.mos_errstate

    def read_NWHW(self):
        self.NWHW = {}
        #x = np.loadtxt(self.filename, dtype=str)       
        x = np.genfromtxt(self.filename,dtype=str, delimiter = [12,14,14,16])
        
        self.stations = ['pgl_'+y[1:5] for y in x[:,0]]
        data = {}
        timearray = {}
        for idx,name in enumerate(self.stations):  
            if name not in data:
                data[name]      = []
                timearray[name] = []
                       
            er    = x[idx,0][-3:-1] 
            dt_obj_Timezoned = timezone('Europe/Berlin').localize(
                    dt.datetime.strptime(x[idx,1],"%Y%m%d%H%M%S"))
            utc_datetime = dt_obj_Timezoned.astimezone(timezone('UTC'))  
            utc_datetime = utc_datetime.replace(tzinfo=None)# Remove timezone Information
            timearray[name].append(utc_datetime)

            values = [float(x.replace(',','.')) for x in 
                      re.findall(r"[+-]?\d+(?:\,\d+)?",x[idx,3])]
            
            if len(values) >1:
                if values[0] > values[1]:
                    value_low,value_high = values[1],values[0]
                else:
                    value_low,value_high = values[0],values[1]
            else:
                value_low,value_high = values[0],values[0]
                
            data[name].append((er,value_low,value_high))
            
        self.stations = np.unique(self.stations)
        for name in self.stations:
            data[name] = np.array(data[name])
            data[name] = pd.DataFrame(data[name],
                     index = timearray[name],
                     columns = ['Ereignis','Von', 'Bis'])      
            data[name][["Von", "Bis"]] = \
                        data[name][["Von", "Bis"]].apply(pd.to_numeric)
            data[name]['avg']   = data[name][['Von', 'Bis']].mean(axis=1)
            
            # Range festlegen aus Textbaustein Vorhersageintervall
            data[name].loc[abs(data[name]['avg']) <= .5, 'range'] = .1
            data[name].loc[(abs(data[name]['avg']) > .5) & (data[name]['avg'] <= 1.5), 'range'] = .125
            data[name].loc[abs(data[name]['avg']) > 1.5, 'range'] = .25
            
            self.NWHW[name] = data[name]             
     
    def load_mos(self):
        """ Loads MOS. Attaches
        the data to attributes of class instance. Calculates trend 
        from data.
        """           
        # MOS and Astronomic Prediction
        fh = open(self.mos_filename,"r")
        init_time = []
        datestamp = []
        hours = []
        minutes = []
        pegel = []
        r0 = []
        r1 = []
        r2 = []
        r3 = []
        r4 = []
        r5 = []
        results = {}
        header1 = fh.readline() 
        for line in fh:
            # Speichern der einzelnen Variablen (Formspezifisch)
            line = line.strip() # entfernen von \t und \n 
            datestamp.append(line[7:15])
            init_time.append(line[15:19])
            hours.append(float(line[23:26]))
            minutes.append(float(line[26:28]))
            pegel.append(float(line[40:43])*0.01)
            r0.append(float(line[43:48])*0.01)
            r1.append(float(line[48:53])*0.01)
            r2.append(float(line[53:58])*0.01)
            r3.append(float(line[58:63])*0.01)
            r4.append(float(line[63:68])*0.01)          
            r5.append(float(line[68:73])*0.01)          
        fh.close()                      
        for keyword in ['pegel','init_time','hours','minutes','r0',
                        'r1','r2','r3','r4','r5', 'datestamp']:
                results[keyword]=np.array(eval(keyword))                   

        self.check_mos(results)
                    
    def plausability(self, mos_dict, new_r):
        """ Checks MOS on plausability.
        
        :param mos_dict: Dictionary with MOS-Data
        :param new_r: Valid "Rückfall-Positionen" (Ohne Fehlerwerte)
        :type mos_dict: Dictionary
        :type new_r: List
        
        :return: If MOS-Data is plausible
        :rtype: Bool
        """
        
        valid_r = np.zeros((len(new_r),len(mos_dict[new_r[0]])))
        for idx,r in enumerate(new_r):
                valid_r[idx] = mos_dict[r]            
        # Wenn das Mittel der ersten 10 Werte von zwei (gültigen)
        # Rückfallpositionen um mehr als 1 Meter abweicht    
        if (len(new_r) > 2) and \
        (abs(np.nanmean(valid_r[0,:10])-np.nanmean(valid_r[1,:10])) > 1) :
                return True 
            
        # Wenn die verwendete Rückfallposition mehr als 80 Zentimeter 
        # vom Staumodell abweicht
        if ('r5' in new_r) and (len(new_r) > 2) and \
        (np.nanmax(abs(valid_r[0]-valid_r[new_r.index('r5')])) > .8):
            return True    
        
        stdds = np.zeros((len(mos_dict[new_r[0]])))
        dr_dt = np.zeros((len(mos_dict[new_r[0]])-1))
        for i, _ in enumerate(mos_dict[new_r[0]]):   
            stdds[i] = np.std(valid_r[:,i])               
            if i < len(mos_dict[new_r[0]])-1: 
                dr_dt[i]    = mos_dict[new_r[0]][i+1]-mos_dict[new_r[0]][i]  
                
        # Wenn die Standardabweichung zw. den Rückfallpositionen > 1m ist 
        # oder die Änderung im Stau >80cm/15min        
        if  np.nanmax(stdds) > 1 or np.nanmax(dr_dt) > 0.8:
                return True       
        return False    
                                   
    def check_mos(self,mos_file_data):
        """ Checks MOS Input (calling plausability method) and creates
        array of datetime-strings which is attached as attribute to 
        the class instance.
        
        :param mos_file_data: Dictionary of MOS-Data 
        :type mos_file_data: Dictionary   
        
        :return: If Mos data is plausible
        :rtype: Bool
           
        """
        all_r  = ['r0','r1','r2','r3','r4','r5']
        boolean_arr = []
        for rueckfall in all_r:
            # Falls Fehlerwert in Zeitreihe, anderer Rückfall
            boolean_arr.append(np.all(abs(mos_file_data[rueckfall]) < abs(50)))    
        new_r  = list(compress(all_r,boolean_arr))
        
        if not new_r or self.plausability(mos_file_data, new_r):    
            mos_WL   = mos_file_data['pegel'] #+ self.PNP     
            print(self.pegel_name + ': MOS X') 
            print(self.pegel_name + ': -> Nur Astronomie')                                                                  
            self.mos_errstate[self.pegel_name] = 1
        else:                
            # Astronomischer Vorhersage
            mos_WL  = (mos_file_data['pegel']  +
                       mos_file_data[new_r[0]] )                  
            self.mos_errstate[self.pegel_name] = 0
            
        # MOS / Astronomischer Zeit Array               
        time_array =  [(dt.datetime.strptime(mos_file_data['datestamp'][x]+
                        mos_file_data['init_time'][x],'%Y%m%d%H%M') + 
                        dt.timedelta(minutes=mos_file_data['minutes'][x]) +
                        dt.timedelta(hours=mos_file_data['hours'][x])) 
                        for x in range(len(mos_file_data['minutes']))]     

        self.mos_data[self.pegel_name] =  pd.DataFrame(mos_WL,
                     index = time_array,
                     columns = ['waterlevel'])   

    
    def get_obs(self,download = False,zeitpunkt = dt.datetime(2020,2,1,12,16)):
        self.po = {}
        self.obs_errstate = {}
        if download:
            for key in self.stations:            
                try:
                    store_dir = os.path.join(os.getcwd(), self.aim_dir)
                    self.po[key] = obs.pegelonline_timeseries(key, store_dir)
                    self.obs_errstate[key] = 0                    
                except:
                    print(key + ': No observation data')
                    self.obs_errstate[key] = 1
        else:            
            p = df.pegelonline(zeitpunkt, freq=1)
            p.store_data_part()
            
        
            inp = pd.read_hdf("/home/bm1483/Development/manuelle_kurve"
                        "/scenarios/pon_data.h5", 'table')*0.01
                              
            for key in self.stations:
                try:    
                    po_o = Object()
                    po_o.data = pd.DataFrame(inp[key].to_numpy(),
                                             index=inp[key].index,
                                             columns = ['waterlevel'])
                    self.po[key] = po_o
                    self.obs_errstate[key] = 0
                except:
                    print(key + ': No observation data')
                    self.obs_errstate[key] = 1
                     
    def prepare_data(self):
        # Dont get lost of data origin
        self.tracedata = {}
        # Resample data:
        #  neuer Index: 
        for key in self.mos_data:           
            # Reorganize MOS _data
            # Resample frequency to 1 Minute
            # Interpolate linearly in time dimension
            self.mos_data[key] = self.mos_data[key].resample('60S').asfreq()   
            self.mos_data[key] = self.mos_data[key].interpolate(method='linear')  
            
            # Add Observation and create one time series
            # Create complete index from start of observation to end of
            # mos / ast forecast
            if self.obs_errstate[key] < 1:
                c_df = pd.concat([self.po[key].data,
                              self.mos_data[key][~self.mos_data[key].index.isin(self.po[key].data.index)]])             
            else:
                c_df = self.mos_data[key]
                
            # Verschiebung der MOS-Kurve beim Übergang von Beobachtung zu
            # Vorhersage
            first_delta = 0
            if self.obs_errstate[key] < 1:
                i = 0
                x = np.nan
                
                # Liegen nans zwischen letzter Beobachtung und erster Vorhersage?
                while np.isnan(x):
                    try:
                        x = self.mos_data[key]['waterlevel'].loc[
                            self.po[key].data.index[-1]
                            + dt.timedelta(minutes=i)]
                    except:
                        i+=1
                        continue
                    i+=1
                    
                # Erste Verschiebung berechnen
                first_delta = (x
                                - self.po[key].data.at[
                                        self.po[key].data.index[-1],
                                        'waterlevel'])
                            
                if np.isnan(first_delta):
                    first_delta = 0
                
                # Wenn Verschiebung sehr groß, dann die NANs auffüllen 
                # (linear interpoliert), und erst dann Kurven verbinden:
                if abs(first_delta) > .5:
                    # Extrapolate
                    odt = (self.po[key].data.at[
                            self.po[key].data.index[-1],
                            'waterlevel']-
                          self.po[key].data.at[
                                  self.po[key].data.index[-2],
                                  'waterlevel'])
                    
                    for m in range(i+1):
                        c_df.loc[self.po[key].data.index[-1] 
                        + dt.timedelta(minutes=i-1)] \
                        = (c_df.loc[self.po[key].data.index[-1]] + m*odt)
           
            # Kurven zusammenführen und Nans auffüllen
            mos = pd.DataFrame(c_df, index = c_df.index.copy())
            mos = mos.interpolate(method='linear') 
            self.mos_data[key] = pd.DataFrame(c_df, 
                                              index = c_df.index.copy(),
                                              columns = c_df.columns.copy())
            
            # Wenn Verschiebung sehr groß, dann nach Interpolationsvorgang
            # die Daten wieder "zurückschieben", da im späteren Schritt
            # die eigentliche Verschiebung stattfindet
            if abs(first_delta) > .5:
                self.mos_data[key]['waterlevel'].loc[
                        (self.po[key].data.index[-1] +
                         dt.timedelta(minutes=1)):
                        (self.po[key].data.index[-1] +
                         dt.timedelta(minutes=i-1))] += first_delta - m*odt 
                mos['waterlevel'].loc[
                        (self.po[key].data.index[-1] +
                         dt.timedelta(minutes=1)):
                        (self.po[key].data.index[-1] +
                         dt.timedelta(minutes=i-1))] += first_delta - m*odt
                
            # Merge MHW, MNW to Manual forecast dataframe
            for i in ['Von','Bis', 'avg']:
                for j in ['NW','HW']:
                    if j == 'NW':
                        self.NWHW[key].loc[self.NWHW[key]['Ereignis'] == j,
                                  i] += self.MNW[key]
                        
                    else:
                        self.NWHW[key].loc[self.NWHW[key]['Ereignis'] == j,
                                  i] += self.MHW[key]  
             
            # Manual forecast on MOS Index
            mf = pd.DataFrame(self.NWHW[key].copy(), index=mos.index.copy(),
                              columns = self.NWHW[key].columns.copy()) 
            
            # Indexcolumn for valid indizes on mos index
            mf['IDX']      = np.arange(len(mf.index))
            mf             = mf.dropna()
            self.NWHW[key] = mf.copy()
            
            # Manuelle Vorhersage auf MOS HW/NW Zeitpunkt verschieben
            new_index = []
            #if self.mos_errstate[key] < 1:
            for idx,i in enumerate(mf['IDX']):   
                if mf.iloc[idx,mf.columns.get_loc('Ereignis')] == 'HW': 
                    # Find local peaks
                    x = mos.iloc[np.max([0,i-150]):i+150,
                                 mos.columns.get_loc('waterlevel')].idxmax()
                else:
                    x = mos.iloc[np.max([0,i-150]):i+150,
                                 mos.columns.get_loc('waterlevel')].idxmin()   
                new_index.append(x)
                
            #Reindex manual forecast
            mf.index = new_index
            mf = pd.DataFrame(mf.copy(), index=mos.index.copy(),
                  columns = mf.columns.copy()) 
            mf['IDX'] = np.arange(len(mf.index))
            
            # Punktauswahl für GPR
            mf             = mf.dropna()
            self.NWHW[key] = mf.copy()
                
            # Ignore 1. Manual forecast, if it is to close to observation
            # or if it is in the past
            if self.obs_errstate[key] < 1:
                while (mf.index[0] - self.po[key].data.index[-1]) < dt.timedelta(hours=2):
                    mf = mf.iloc[1:]
                else:
                    pass
            else:
                if self.mos_errstate[key] < 1:
                    if (mf.index[0] - mos.index[0]) < dt.timedelta(hours=2):
                        mf = mf.iloc[1:]
                else:
                    pass
                
            # Calculate difference between curves at manual forecast 
            # times and difference between obsveration and first value of mos
            # New Dataframe von Beobachtung bis 8h hinter 
            # letztem manuellen Ereignis
            # Observation has highest accuracy apart from manual forecast
            if self.obs_errstate[key] < 1:
                o_idx = len(self.po[key].data.index)-1
                #mos.iloc[:o_idx,mos.columns.get_loc('DY')] = 0.0001
            else:
                #mos.iloc[:60,mos.columns.get_loc('DY')] = 0.0001
                o_idx = 0    
                
            if self.obs_errstate[key] < 1:
                delta = pd.DataFrame(mos.iloc[o_idx+1,mos.columns.get_loc('waterlevel')]
                                   -self.po[key].data['waterlevel'][-1],index=[mos.index[o_idx+1]],
                                   columns = ['waterlevel'])
            else: 
                delta = pd.DataFrame(0,index=[mos.index[o_idx+1]],
                                     columns = ['waterlevel'])
                
            for lm,l in enumerate(mf['IDX']):
                delta = delta.append(pd.DataFrame(mos.iloc[l,mos.columns.get_loc('waterlevel')]
                              -mf.iloc[lm,mf.columns.get_loc('avg')],index=[mos.index[l]],
                              columns = ['waterlevel']))
      
            # 8 hours after last HW/NW no delta
            delta = delta.append(pd.DataFrame(0,
                                 index=[mos.index[l]+dt.timedelta(hours=8)],
                                 columns = ['waterlevel']))                
            delta = delta.dropna()                
                
            # Interpolate differences linearly
            delta = delta.resample('60S').asfreq()   
            delta = delta.interpolate(method='linear')                             
            # Substract difference from mos
            self.old_mos = mos.copy()         
            
            mos['waterlevel'] = mos['waterlevel'].sub(
                                delta['waterlevel'], fill_value = 0)       
            mos = mos.resample('60S').asfreq()   
            mos = mos.interpolate(method='linear')                         
            self.complete_mos = mos.copy()
            
            # Start GPR
            self.GPR_MANU(key, mf, mos) 
              
        
    def GPR_MANU(self, key, mf, mos):
        # Glättung
        smooth = savgol_filter(self.complete_mos['waterlevel'].to_numpy(),
                               135,2)  
        #smooth = self.complete_mos['waterlevel'].to_numpy()
        # Plot results
        pd.plotting.register_matplotlib_converters()
        
        fig = plt.figure(figsize=(10, 5))
        lw = 3
        colors = {'NWHW': "#1f77b4",'lin':'#2ca02c','Pchip':'#9467bd',
                  'mos':'navy'}
        for idx,k in enumerate(self.tracedata):
            plt.scatter(mos.index[np.array(self.tracedata[k][0],
                                                         dtype=int)],
                        self.tracedata[k][1],
                        label=k, s = 6, color = colors[k], zorder = 100-idx)
                
        if self.obs_errstate[key] < 1:
            plt.plot(self.po[key].data.index,self.po[key].data['waterlevel'],
                 color='r',ms = 3,lw=lw, label='OBS',zorder=100)
        
        plt.plot(self.old_mos.index, self.old_mos,
                 color='navy',lw=lw-1.5, label='MOS',zorder=11)
        
        plt.plot(mos.index, smooth,
                 color='c',lw=lw, label='Manuelle Kurve',zorder=10)    
        
        line,caps,bars=plt.errorbar(mos.index[
                        self.NWHW[key]['IDX']],     # X
                        self.NWHW[key]['avg'],    # Y
                        yerr= self.NWHW[key]['range'],# Y-errors
                        color="darkorange",    # format line like for plot()
                        marker = 's',
                        linewidth=0,   # width of plot line
                        elinewidth=2,  # width of error bar line
                        ecolor='grey',    # color of error bar
                        capsize=2,     # cap length for error bar
                        capthick=1,   # cap thickness for error bar
                        label='Manuell',
                        ms = 6,
                        zorder = 111                        
                        )
        ax = plt.gca()
        ax.axhline(self.MNW[key],LineStyle='--',Color='k',
                   Linewidth=0.5,zorder=0, alpha=.6)
        ax.axhline(self.MHW[key],LineStyle='--',Color='k',
                   Linewidth=0.5,zorder=0, alpha=.6)
        
        
        ax.set_xlim([self.NWHW[key].index[0] - dt.timedelta(days = .75),
                     self.NWHW[key].index[0] + dt.timedelta(days = 1.5)])
        
        # Format Time axis -> Hours on top, Days at bottom of plot          
        ax.xaxis.set_major_locator(DayLocator())   
        ax.xaxis.set_major_formatter(DateFormatter('%a %d.%m.%y'))                        
        ax.xaxis.set_minor_locator(HourLocator(byhour=[0,6,12,18]))    
        ax.xaxis.set_minor_formatter(DateFormatter('%H (UTC)'))    
        ax.grid(which='major', axis='x',)

        # rotates and right aligns the x labels, and moves the bottom of the
        # axes up to make room for them
        fig.autofmt_xdate()
        
        #plt.xlabel('Minutes')
        plt.ylabel('Waterlevel')
        plt.title('GPR: ' +key)
        plt.legend(loc="upper right", scatterpoints=1, prop={'size': 6})
        plt.savefig(os.path.abspath('../../data/figs/'+key+'.png'),dpi=96)
        #plt.show()     

        if self.obs_errstate[key] < 1:
            start = self.po[key].data.index[-1] + dt.timedelta(minutes=1)
        else:
            start = mos.index[0]
 
        ypred_df = pd.DataFrame(smooth,
                                index = mos.index,
                                columns = [key+'_MC'])         
        
        store = pd.concat([ypred_df[key+'_MC'],
                           self.NWHW[key]['avg'],
                           self.NWHW[key]['range']], axis=1,
                           keys=[key+'_MC',
                                 key+'_MAN',
                                 key+'_MANUNC'])
       
        
        store.index.name = 'Timestamp'
        
        store = store.loc[start:
                self.NWHW[key].index[0] + dt.timedelta(days = 1.5)]
        store.to_hdf(os.path.abspath('../../data/output/'+key+'.h5'),
                             'CombinedForecastUTC')

        plt.close('all')

if __name__ == '__main__':
    d = True
    #zeitpunkt = dt.datetime(2017,10,29,7,00)
    zeitpunkt = dt.datetime(2020,3,1,9,00)
    
    x = combine_curves()
    
    print('# GET DATA #')
    stime1 = time.time()
    
    x.get_current_NWHW(download = d,zeitpunkt = zeitpunkt)
    x.read_NWHW()
    print("Manuelle Daten \u2713 : %.1f s" % (time.time() - stime1))
    
    stime = time.time()
    x.get_eformat(download = d)
    print("EFORMAT        \u2713 : %.1f s" % (time.time() - stime))
    
    stime = time.time()
    x.get_mos(download = d,zeitpunkt = zeitpunkt)
    print("MOS            \u2713 : %.1f s" % (time.time() - stime))
    
    stime = time.time()
    x.get_obs(download = d,zeitpunkt = zeitpunkt)
    print("Pegelonline    \u2713 : %.1f s" % (time.time() - stime))
    
    
    print('# Kurven zusammenfügen #')
    x.prepare_data()
    x.clear_old_data()
    x.summarize_figures()
    
    #x.summarize_data()
    
    print("Total time elapsed: %.1f s" % (time.time() - stime1))
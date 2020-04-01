# -*- coding: utf-8 -*-
"""
This script runs the operational manual waterlevel forecast curve construction.
It can be executed from the command line.
"""

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

class Object(object):
    '''
    Hilfsobjekt definieren, um später Daten als Attribut ansprechen zu können.
    '''
    pass

class combine_curves():
    '''
    Hauptklasse mit allen Attributen und Methoden, die notwendig sind, um
    eine manuelle Kurve zu erzeugen.
    
    '''
    def __init__(self):
        '''
        Directories initialisieren. Nicht benötigten Standard-Unix output
        in self.FNULL umleiten.
        
        '''
        self.FNULL = open(os.devnull, 'w')
        self.aim_dir = os.path.abspath('../../data/temp/')
        self.aim_dir = os.path.join(self.aim_dir, '')
        self.trace_data = {} # Für BSH-Homepage alle Daten sammeln
        
    def clear_old_data(self):
        '''
        Alte temporäre Daten bereinigen.
        
        '''
        path  = os.path.abspath('../../data/temp/')
        path  = os.path.join(path, '')
        
        files = os.listdir(path)
        for f in files:
            if f.endswith('.dat') or f.endswith('.txt'):
                x = path + f
                os.remove(x) 
    
    def summarize_figures(self):
        '''
        Optional function to combine all figures into 1 pdf: requires PIL
        
        '''
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
            imgs_comb = np.hstack( [np.asarray( i.resize(min_shape) ) 
                                    for i in ims ] )
        
            # save that beautiful picture
            imgs_comb = PIL.Image.fromarray( imgs_comb)
            outfile   = os.path.abspath("../../data/output/combined"+
                                        str(i)+".png")
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
        '''
        Methode, um alle Daten zusammenzufassen. Sie werden in das Hdf5-File
        des aktuellen Tages, mit einem Zeitstempel versehen, geschrieben.
        
        '''
        x = glob.glob(os.path.abspath('../../data/output/*.h5'))
        today = dt.datetime.now().strftime('%Y%m%d')
        new_file = '../../data/data/'+today+'.h5'
        new_f = pd.DataFrame()
        
        for i in x:
            new_f = pd.concat([new_f,
                               pd.read_hdf(i,'CombinedForecastUTC')], axis=1)
        
        # Store Data in file and create frame with init information
        new_f.to_hdf(new_file,mode='a',complevel = 1, 
                     key=dt.datetime.now().strftime('T%Y%m%d%H%M'))
           

    def get_current_NWHW(self, download = True,
                         zeitpunkt = dt.datetime(2020,2,1,12,16)):       
        ''' Neueste manuelle Vorhersage oder archivierte Vorhersage beziehen.
        
        :param download: Schalter, ob download der aktuellsten Daten oder Zeitpunkt
        :param zeitpunkt: Falls kein download, gewünschter Zeitpunkt
        :type download: Bool
        :type zeitpunkt: datetime object
        
        '''
        self.download = download
        if download:
            # NWHW
            cmd = ("ssh bm11mos@linwvd10 ls "
                   "/FEWS/HWNW_Vorhersage/HWNW* | tail -1")    
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
            aim_file =  os.path.join(self.aim_dir,
                                     filename.rsplit('/', 1)[-1].strip())
            
        else:            
            x = df.nwhw(zeitpunkt)
            aim_file = x.valid_file
            
        self.filename = aim_file
   
    def get_eformat(self, download = False):
        ''' E-Format des Jahres beziehen. Für mittleres HW|NW.
        
        :param download: Schalter, ob download der aktuellsten Daten oder Zeitpunkt
        :type download: Bool
        
        '''
        self.dirname = os.path.dirname(__file__)
        self.aim_dir = os.path.abspath('../../data/temp/')
        # Add trailing slash
        self.aim_dir = os.path.join(self.aim_dir, '') 
        
        self.MHW = {}
        self.MNW = {}
        self.PNP = {}
        self.SKN = {}
        for pegel in self.stations:
            # Immer dieses und nächstes Jahr berücksichtigen
            year = dt.datetime.now().strftime('%Y')
            
            # EFORMAT: Mittlere Hoch-Niedrigwasser
            if download:
                cmd = ("ssh bm14gast@lingezeiten11 ls /data/vb_hwnw/deu"
                       +year+"/*"+pegel[-4:]+"*")  
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
                self.read_EFORMAT_file(filename,pegel)  
            else:
                break
            
        if not download:
            # Get Special Event / year / archive
            files = os.listdir(os.path.abspath('../../data/arch/EFORMAT'))
            pnames = ['pgl' + x[-13:-8] for x in files]
            files  = [os.path.join(os.path.abspath(
                      '../../data/arch/EFORMAT'),x) for x in files]           
            for idx,source_file in enumerate(files):            
                shutil.copy2(source_file,self.aim_dir)
                self.read_EFORMAT_file(source_file, pnames[idx])
        

    def read_EFORMAT_file(self,filename,pegel):
        ''' Liest MHW|MNW aus EFORMAT File.
        
        :param filename: Dateiname
        :param pegel: Falls kein download, gewünschter Zeitpunkt
        :type filename: Bool
        :type pegel: datetime object
        
        '''
        with open(filename, mode = 'rb') as f:
            for idx,line in enumerate(f):
                if idx == 16:
                    self.PNP[pegel] = \
                        float(re.findall(r'\d.\d\d',str(line))[0])
                if idx == 17:
                    self.SKN[pegel] = \
                        float(re.findall(r'\d.\d\d',str(line))[0])
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
        ''' Aktuellste MOS-Daten beziehen oder aus Archiv laden.
        
        :param download: Schalter, ob download der aktuellsten Daten oder Zeitpunkt
        :param zeitpunkt: Falls kein download, gewünschter Zeitpunkt
        :type download: Bool
        :type zeitpunkt: datetime object
        
        '''
        self.mos_data = {}
        self.old_mos = {}
        #self.ast_data = {}        
        self.mos_errstate = {}
        
        if download:
            for pegel in self.stations:            
                # MOS
                cmd = ("ssh bm11mos@linwvd10 ls "
                       "/data_MOS/bm11mos/1704/output/*" +
                       pegel+"* | tail -1")
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
                self.old_mos[pegel] = self.mos_data[pegel].copy()  
                        
        else:
            x = df.mos(zeitpunkt = zeitpunkt, freq=15)
            x.main()
            pname = os.path.join(os.path.abspath("../../data/scenarios/"),"")
            source_file = (pname + "mos_data.h5")          
            data = pd.read_hdf(source_file,'table')
            for pegel in self.stations:
                self.mos_data[pegel] = data[pegel]  
                self.old_mos[pegel] = self.mos_data[pegel].copy()  
            self.mos_errstate = x.mos_errstate
            

               
    def read_NWHW(self):
        '''Manuelle Vorhersage-Datei lesen
        
        '''
        self.NWHW = {}     
        x = np.genfromtxt(self.filename,dtype=str, delimiter = [12,14,14,16])
        
        self.stations = ['pgl_'+y[1:5] for y in x[:,0]]
        data = {}
        timearray = {}
        for idx,name in enumerate(self.stations):  
            if name not in data:
                data[name]      = []
                timearray[name] = []
                   
            if idx== 0:
                dt_obj_Timezoned = timezone('Europe/Berlin').localize(
                    dt.datetime.strptime(x[idx,2]," %Y%m%d%H%M "))
                utc_datetime = dt_obj_Timezoned.astimezone(timezone('UTC'))  
                utc_datetime = utc_datetime.replace(tzinfo=None)     
                self.nwhw_init = utc_datetime
                
            er    = x[idx,0][-3:-1] 
            dt_obj_Timezoned = timezone('Europe/Berlin').localize(
                    dt.datetime.strptime(x[idx,1],"%Y%m%d%H%M%S"))
            utc_datetime = dt_obj_Timezoned.astimezone(timezone('UTC'))  
            utc_datetime = utc_datetime.replace(tzinfo=None)                   # Remove timezone Information
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
            data[name].loc[(abs(data[name]['avg']) > .5) &
                           (data[name]['avg'] <= 1.5), 'range'] = .125
            data[name].loc[abs(data[name]['avg']) > 1.5, 'range'] = .25
            
            self.NWHW[name] = data[name]             
     
    def load_mos(self):
        """ Loads MOS. Starts checking routines (plausability).
        Attaches the data to attributes of class instance.
        
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
        
        # Raw data speichern
        self.trace_data['mos'+self.pegel_name] =  \
                                    self.mos_data[self.pegel_name].copy()*100     
        self.trace_data['ast'+self.pegel_name] = \
                                    pd.DataFrame(mos_file_data['pegel']*100,
                                              index = time_array,
                                              columns = ['waterlevel'])     
        
    def get_obs(self,download = False,zeitpunkt = dt.datetime(2020,2,1,12,16)):
        ''' Beobachtungsdaten etnweder direkt von Pegelonline laden oder
        aus dem Archiv beziehen. Für den Download-Fall das Hilfspaket laden
        und dort werden die Plausabilitätsprüfungen durchgeführt.
        
        :param download: Schalter, ob download der aktuellsten Daten oder Zeitpunkt
        :param zeitpunkt: Falls kein download, gewünschter Zeitpunkt
        :type download: Bool
        :type zeitpunkt: datetime object
        
        '''
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
                if np.isnan(self.po[key].data.to_numpy()).all():
                    print(key + ': No observation data')
                    self.obs_errstate[key] = 1   
        else:        
            p = df.pegelonline(zeitpunkt, freq=1)
            p.store_data_part()
            
            pname = os.path.join(os.path.abspath("../../data/scenarios/"),"")
            inp = pd.read_hdf(pname + "pon_data.h5", 'table')*0.01
                              
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
                    
                if np.isnan(self.po[key].data.to_numpy()).all():
                    print(key + ': No observation data')
                    self.obs_errstate[key] = 1              
                     
    def construct_curve(self):
        ''' Hauptmethode, die aus den Eingangsdaten eine manuelle Kurve
        erzeugt.
        
        '''
        # Dont get lost of data 
        self.inp_mos = {}
        for key in self.mos_data:           
            # Reorganize MOS _data
            # Resample frequency to 1 Minute
            # Interpolate linearly in time dimension
            self.mos_data[key] = self.mos_data[key].resample('60S').asfreq()   
            self.mos_data[key] = self.mos_data[key].interpolate(method='linear')
            
            # Eingangsdaten glätten
            self.mos_data[key] = self.mos_data[key].dropna()
            self.mos_data[key]['waterlevel'] = savgol_filter(
                            self.mos_data[key]['waterlevel'].to_numpy(),141,2)  
            
            #Interpolierte Daten behalten
            self.inp_mos[key] = self.mos_data[key]
            
            # Add Observation and create one time series
            # Create complete index from start of observation to end of
            # mos / ast forecast
            if self.obs_errstate[key] < 1:
                c_df = pd.concat([self.po[key].data,
                              self.mos_data[key][
                                      ~self.mos_data[key].index.isin(
                                              self.po[key].data.index)]])             
            else:
                c_df = self.mos_data[key]
                
            # Datenlücken zwischen Beobachtung und MOS/Astronomie behandeln
            if self.obs_errstate[key] < 1:
                i = 0
                x = np.nan
                
                # Liegen nans zwischen letzter Beobachtung und erster Vorhersage?
                while np.isnan(x):
                    try:
                        x = self.mos_data[key]['waterlevel'].loc[
                            self.po[key].data.last_valid_index()
                            + dt.timedelta(minutes=i)]
                    except:
                        i+=1
                        continue
                    i+=1
                    
                # NANs auffüllen und erst dann Kurven verbinden:
                if i < 60:
                    # Extrapolate Observation
                    odt = (self.po[key].data.at[
                            self.po[key].data.index[-1],
                            'waterlevel']-
                          self.po[key].data.at[
                                  self.po[key].data.index[-2],
                                  'waterlevel']) 

                    ts_1 = self.po[key].data.index[-1] 
                    last_val = self.po[key].data['waterlevel'].iloc[-1]
                    
                    for m in range(1,i+1):
                        if m < (i):
                            c_df.loc[ts_1 + dt.timedelta(minutes=m)] \
                            = (last_val+ m*odt)
                        
                        self.po[key].data.loc[ts_1
                        + dt.timedelta(minutes=m)] \
                        = (last_val + m*odt)
         
            # Kurven zusammenführen und Nans auffüllen
            mos = pd.DataFrame(c_df, index = c_df.index.copy())
            mos = mos.resample('60S').asfreq()
            self.mos_data[key] = pd.DataFrame(c_df, 
                                              index = c_df.index.copy(),
                                              columns = c_df.columns.copy())
                
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
            # Hierzu die lokalen Maxima/Minima um den manuellen Zeitpunkt finden
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
            
            # Nans rauswerfen
            mf             = mf.dropna()
            self.NWHW[key] = mf.copy()
                
            # Ignore 1. Manual forecast, if it is to close to observation
            # or if it is in the past
            if self.obs_errstate[key] < 1:
                while ((mf.index[0] - self.po[key].data.index[-1])
                        < dt.timedelta(hours=2)):
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
            
            # Last Obersavtion Index 
            if self.obs_errstate[key] < 1:
                o_idx = self.po[key].data.index[-1]
            else:
                o_idx = mos.index[0]  
            
            # Erster Differenzwert
            if self.obs_errstate[key] < 1:    
                try:                   
                    delta = pd.DataFrame(
                            self.inp_mos[key].at[
                            o_idx,'waterlevel']-      
                            self.po[key].data.at[o_idx,'waterlevel'],
                            index=[o_idx],
                            columns = ['waterlevel'])
                except:
                    delta = pd.DataFrame(0,index=[o_idx],
                                     columns = ['waterlevel'])
            else: 
                delta = pd.DataFrame(0,index=[o_idx],
                                     columns = ['waterlevel'])
                
            # Für die Ereignisse Differenzen berechnen
            for lm,l in enumerate(mf['IDX']):
                delta = delta.append(pd.DataFrame(mos.iloc[
                              l,mos.columns.get_loc('waterlevel')]
                              -mf.iloc[lm,mf.columns.get_loc('avg')],
                              index=[mos.index[l]],
                              columns = ['waterlevel']))
          
            # 8 hours after last HW/NW no delta
            delta = delta.append(pd.DataFrame(0,
                                 index=[mos.index[l]+dt.timedelta(hours=8)],
                                 columns = ['waterlevel']))       
            delta = delta.fillna(0)
                
            # Interpolate differences linearly
            delta = delta.resample('60S').asfreq()   
            delta = delta.interpolate(method='linear')
            
            # Restructure and fill Nans with 0
            delta = pd.DataFrame(delta, index = mos.index,
                                 columns = ['waterlevel'])
            delta = delta.fillna(0)
            
            # Substract difference from mos            
            mos['waterlevel'] = mos['waterlevel'].sub(
                                delta['waterlevel'])       
            
            mos = mos.resample('60S').asfreq()
            self.complete_mos = mos.copy()
            
            # Start GPR
            self.plotting_and_output(key, mf, mos) 
              
        
    def plotting_and_output(self, key, mf, mos):
        '''Daten darstellen und Output schreiben.
        
        :param key: pegelname
        :param mf:  Ursprüngliche manuelle Vorhersage
        :param mos: Gesamte manuelle Daten
        
        :type key: strin
        :type mf: pandas DataFrame
        :type mos: pandas DataFrame
        
        '''
        # Als numpy array für Schreibroutine
        smooth = self.complete_mos['waterlevel'].to_numpy()       

        # Plot results
        pd.plotting.register_matplotlib_converters()
        
        fig = plt.figure(figsize=(10, 5))
        lw = 3
                
        if self.obs_errstate[key] < 1:
            plt.plot(self.po[key].data.index,self.po[key].data['waterlevel'],
                 color='r',ms = 3,lw=lw, label='OBS',zorder=100)
        
        plt.plot(self.old_mos[key].index, self.old_mos[key],
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
        
        
        ax.set_xlim([self.NWHW[key].index[0] - dt.timedelta(days = .25),
                     self.NWHW[key].index[-1] + dt.timedelta(days = .25)])
        
        '''
        ax.set_xlim([self.po[key].data.index[-1] - dt.timedelta(days = .025),
                     self.po[key].data.index[-1] + dt.timedelta(days = .025)]) 
        ax.set_ylim([np.min(self.complete_mos['waterlevel'].loc[
                     self.po[key].data.index[-1] - dt.timedelta(days = .025):
                     self.po[key].data.index[-1] + dt.timedelta(days = .025)])-.25,
                     np.max(self.complete_mos['waterlevel'].loc[
                     self.po[key].data.index[-1] - dt.timedelta(days = .025):
                     self.po[key].data.index[-1] + dt.timedelta(days = .025)])+.25])
        '''
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
        
        # Für Internet
        # Spezielle CSV Datei für Tableau Output: UTC + Zentimeter
        #bshnr	 datum	 messung	 astro	 stauPN	 mosdatum	 pos	 NHN	 SKN
        
        if self.download:
            # Neuer Zeitindex
            beob_start = self.po[key].data.index[-1] - dt.timedelta(hours=14)
            manc_end   = self.NWHW[key].index[-1] + dt.timedelta(days = .5)
            min_freq   = pd.date_range(beob_start,manc_end, freq='60S')
            comb_freq  = min_freq.union(self.trace_data['mos'+key].index[
                        ~self.trace_data['mos'+key].index.isin(min_freq)])
                    
            #1. Create dataframe with all data
            store = pd.DataFrame({'bshnr': key[4:],
                   'datum' : comb_freq,
                   'messung': (
                       self.po[key].data['waterlevel']*100).astype(int),
                   'astro':(
                       self.trace_data['ast'+key]['waterlevel']).astype(int),
                   'stauPN':(
                       self.trace_data['mos'+key]['waterlevel']).astype(int),
                   'mosdatum': dt.datetime.now().replace(microsecond=0), 
                   'pos':'Normal',
                   'NHN': (self.trace_data['mos'+key]['waterlevel']-
                           self.PNP[key]*100).astype(int) ,
                   'SKN': (self.trace_data['mos'+key]['waterlevel']-
                           self.SKN[key]*100).astype(int) ,
                   'manuell': (self.complete_mos['waterlevel'].loc[
                           self.po[key].data.index[-1]:]*100).astype(int)},
                    index=comb_freq)
            
            store.index.name = 'datum'
            
            if os.path.exists('../../data/output/'+'tableau_Nordsee31.csv'):
                header = False
            else:
                header = True
            
            store.to_csv(
                    os.path.abspath('../../data/output/'+
                                    'tableau_Nordsee31.csv'),
                    float_format ="%.0f", mode = 'a', index = False,
                    header = header)
            
            # Create manual Tableau file
            #1. Create dataframe with all data
            #bshnr	 Datum	 Wert	 Abweichung	 Vorhersage	 MOS_UTC	 MOS_GZ	 Warnung	Notiz
            store = pd.DataFrame({'bshnr': key[4:],
                   'Datum' : self.NWHW[key].index,
                   'Wert': (self.NWHW[key]['avg']*100).astype(int),
                   'Abweichung':(self.NWHW[key]['range']*200).astype(int),
                   'Vorhersage': self.nwhw_init,
                   'MOS_UTC': dt.datetime.strptime(
                        self.mos_filename.split('.')[0][-17:-5],'%Y%m%d%H%M'), 
                   'MOS_GZ':  (dt.datetime.strptime(
                        self.mos_filename.split('.')[0][-17:-5],'%Y%m%d%H%M')+
                        dt.timedelta(hours=1)), 
                   'Warnung': 'Normal',
                   'Notiz': ''},
                    index=self.NWHW[key].index)
            
            store.to_csv(
                    os.path.abspath('../../data/output/'+
                                    'tableauHWNW_Vorhersage31.csv'),
                    float_format ="%.0f", mode = 'a', index = False,
                    header = header)

if __name__ == '__main__':
    d = True
        
    try:
        os.remove(os.path.abspath('../../data/output/'+
                                  'tableau_Nordsee31.csv'))
        os.remove(os.path.abspath('../../data/output/'+
                                  'tableauHWNW_Vorhersage31.csv'))
    except:
        pass
    #zeitpunkt = dt.datetime(2017,10,29,7,00)
    zeitpunkt = dt.datetime(2020,2,10,9,00)
    stime1 = time.time()
    x = combine_curves()
    
    print('# GET DATA #')   
    
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
    x.construct_curve()
    x.clear_old_data()
    #x.summarize_figures()   
    #x.summarize_data()
    
    print("Total time elapsed: %.1f s" % (time.time() - stime1))
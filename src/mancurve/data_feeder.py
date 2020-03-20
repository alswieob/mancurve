"""
Create test-scenario for specific date:
# Gets observation data until specific time
# Loads MOS-Forecast for fitting time
-> Store scenario data in directory to be loaded by manu_curve_os_independent
"""

__author__ = 'Luis Becker'
__version__ = '02/2020'

import os
import shutil
import subprocess
import glob
import numpy as np
import pandas as pd
import datetime as dt
from itertools import compress

class pegelonline():
    def __init__(self,endzeitpunkt, freq=1):
        self.end  = endzeitpunkt
        self.freq = str(freq)+'min'
        # Um zwei komplette Tage einzulesen
        self.start = self.end - dt.timedelta(days=1)        
        self.timearray =  [self.start.strftime('%Y_%m_%d'),
                           self.end.strftime('%Y_%m_%d')] 
        
        # Metadaten [pegel_ID, (LAT, LON), PNP]
        self.metadata = np.load("/home/bm1483/Pakete/"
                                "validation_pkg/imonav_validation/"
                                "/stationen/stationen.npy",
                                allow_pickle = True).item()               
        self.PNP = {self.metadata[name][0]: self.metadata[name][2]*100
            for name in self.metadata}  

    
    def store_data_part(self):
        '''
        Erzeugt ein HDF5 file, dass die Beobachtungsdaten für den
        Zeitpunkt -48h enthält. Dafür werden die 24h-Dateien
        eingelesen und ein neuer Array erzeugt.
        '''
        for idx,datum in enumerate(self.timearray):             
            path =  '/work/bm1483/pegelonline/'+datum+'/'        
            if idx == 0:
                pegelonline = self.load_po_data(path,datum) #+ self.PNP   
            else:
                x = self.load_po_data(path,datum) #+ self.PNP   
                pegelonline = pegelonline.append(x)

        # New Index
        self.new_idx = pd.date_range(start=pegelonline.index[0],
                                     end=self.end, freq=self.freq)
        
        # Convert numpy array to dataframe
        po_t = pd.DataFrame(pegelonline,
                            columns = self.pegel,
                            index = self.new_idx)
        
        po_t.to_hdf("/home/bm1483/Development/manuelle_kurve"
                    "/scenarios/pon_data.h5", 'table')
        
    def load_po_data(self,path,datum):
        timestamp = datum.replace('_','') 
        data = pd.read_hdf("{0}{1}".format(path,'/Pegelonline_'+
                           timestamp+'.h5'),'OBS')
        self.pegel = data.columns
        return data
    
class mos():
    def __init__(self,zeitpunkt, freq=1):       
        self.time = zeitpunkt
        self.freq = str(freq)+'min'
        self.mos_data     = {}
        self.mos_errstate = {}
        
        # Um zwei komplette Tage einzulesen
        self.datum = self.time.strftime('%Y_%m_%d')        
        self.init  = (self.time.strftime('%H') 
                      + str((self.time.minute // 15)*15).zfill(2))
        
        
        # Metadaten [pegel_ID, (LAT, LON), PNP]
        self.metadata = np.load("/home/bm1483/Pakete/"
                                "validation_pkg/imonav_validation/"
                                "/stationen/stationen.npy",
                                allow_pickle = True).item()               
        self.PNP = {self.metadata[name][0]: self.metadata[name][2]*100
            for name in self.metadata}
      
    def prepare_mosdata(self):
        mos_path = "{0}{1}{2}".format('/work/bm1483/mos_archiv/',
                                      self.datum,'/')
        
        came_from = os.getcwd()
        os.chdir(mos_path)
        
        tar_archives = os.listdir(mos_path)
        tar_archives = [x for x in tar_archives if x.endswith("tar.gz")]
        # Entpacke Archives
        for archive in tar_archives:
            file_name = "{0}{1}".format(mos_path , archive)
            p = subprocess.Popen(['tar', 'xfz', file_name])
            p.wait()
            
        # Alle Dateien zu ganzer Stunde auslesen
        data = os.listdir(mos_path)
        data = [x for x in data if x.endswith(".dat") and "00_" in x]
        init_times = [x[-13:-9] for x in data]
        init_times = np.unique(init_times)
        
        # Zeiten aussortieren, für die nicht alle Pegel vorhanden sind
        valid_times = []
        for idx,tim in enumerate(init_times):
            if len([x for x in data if "{0}{1}".format(tim,'_') in x]) < 30:
                continue
            valid_times.append(tim)
     
        os.chdir(came_from)
        
        # Init auswählen, das am nächsten am Zeitpunkt liegt
        x,i = False,0
        if self.init in valid_times:
            self.valid_init = self.init
        else:
            while x == False:
                if self.init[:2] in valid_times[i]:
                    self.valid_init = valid_times[i]
                    x = True
                else:
                    i+=1                   
        
        # Gültige Files
        self.mos_files = glob.glob(mos_path+'*'+self.valid_init+'_*')
        
    def mos_get_h5_part(self):     
        for open_file in self.mos_files:     
            # Read in Mos-Data
            self.pegel_name = open_file.split('/')[-1][:8]           
            self.get_mos_dictionary(open_file)  
            
                
        self.mos_data = pd.concat(self.mos_data, axis = 1)
        #self.mos_data.columns = self.mos_data.columns.droplevel(1)
        self.mos_data = self.mos_data.resample('60S').asfreq()   
        self.mos_data = self.mos_data.interpolate(method='linear')       
           
        self.mos_data = self.mos_data.resample(self.freq).asfreq()   
        self.mos_data.to_hdf("/home/bm1483/Development/manuelle_kurve"
                             "/scenarios/mos_data.h5", 'table')
        
    def get_mos_dictionary(self,filename):
        """ Liest mos.txt File ein und gibt ein Dictionary mit Daten raus
        
        Args:
            file:       Filename eines Mos files
            
        Returns:
            results     Dictionary mit den MOS-Daten
        """   
        # MOS and Astronomic Prediction
        fh = open(filename,"r")
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
                 
    def plausa_mos(self, mos_dict, new_r):
        """ Checks MOS on plausability. Treat Zollenspieker & Geesthacht with
        higher tresholds.
        
        :param mos_dict: Dictionary with MOS-Data
        :param new_r: Valid "Rückfall-Positionen" (Ohne Fehlerwerte)
        :type mos_dict: Dictionary
        :type new_r: List
        
        :return: If MOS-Data is plausible
        :rtype: Bool
        """
        if self.pegel_name == 'pgl_732D' or self.pegel_name == 'pgl_731P':
            thhold = 2
        else:
            thhold = 1
            
        valid_r = np.zeros((len(new_r),len(mos_dict[new_r[0]])))
        for idx,r in enumerate(new_r):
                valid_r[idx] = mos_dict[r]            
        # Wenn das Mittel der ersten 10 Werte von zwei (gültigen)
        # Rückfallpositionen um mehr als 1 Meter abweicht    
        if (len(new_r) > 2) and \
        (abs(np.nanmean(valid_r[0,:10])-np.nanmean(valid_r[1,:10])) > thhold):
                return True 
            
        # Wenn die verwendete Rückfallposition mehr als 80 Zentimeter 
        # vom Staumodell abweicht
        if ('r5' in new_r) and (len(new_r) > 2) and \
        (np.nanmax(abs(valid_r[0]-valid_r[new_r.index('r5')])) > thhold):
            return True    
        
        stdds = np.zeros((len(mos_dict[new_r[0]])))
        dr_dt = np.zeros((len(mos_dict[new_r[0]])-1))
        for i, _ in enumerate(mos_dict[new_r[0]]):   
            stdds[i] = np.std(valid_r[:,i])               
            if i < len(mos_dict[new_r[0]])-1: 
                dr_dt[i]    = mos_dict[new_r[0]][i+1]-mos_dict[new_r[0]][i]  
                
        # Wenn die Standardabweichung zw. den Rückfallpositionen > 1m ist 
        # oder die Änderung im Stau >80cm/15min        
        if  np.nanmax(stdds) > thhold or np.nanmax(dr_dt) > 0.8:
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
        
        if not new_r or self.plausa_mos(mos_file_data, new_r):    # True for astronomic simulation
            # Astronomische Vorhersage
            mos_WL   = mos_file_data['pegel'] #+ self.PNP    
            print(self.pegel_name + ': -> Nur Astronomie')                                                                  
            self.mos_errstate[self.pegel_name] = 1
        else:                
            mos_WL  = (mos_file_data['pegel']  +
                       mos_file_data[new_r[0]] )#+
                       #self.PNP[self.pegel_name])                   
           
               
            self.mos_errstate[self.pegel_name] = 0
            
        # MOS / Astronomischer Zeit Array               
        time_array =  [(dt.datetime.strptime(mos_file_data['datestamp'][x]+
                        mos_file_data['init_time'][x],'%Y%m%d%H%M') + 
                        dt.timedelta(minutes=mos_file_data['minutes'][x]) +
                        dt.timedelta(hours=mos_file_data['hours'][x])) 
                        for x in range(len(mos_file_data['minutes']))]     # Potenziell langsam

        self.mos_data[self.pegel_name] =  pd.DataFrame(mos_WL,
                     index = time_array,
                     columns = ['waterlevel'])  
            
    def pack_mosdata(self):
        # Entpackte einzeldateien wieder löschen
        mos_path = "{0}{1}{2}".format('/work/bm1483/mos_archiv/',
                   self.datum,'/')
        
        came_from = os.getcwd()
        os.chdir(mos_path)
        
        data = os.listdir(mos_path)
        data = [x for x in data if x.endswith(".dat")]
        	
        for f in data:
            os.remove(f)   
        
        os.chdir(came_from)
            
    def main(self):
        self.prepare_mosdata()
        self.mos_get_h5_part()
        self.pack_mosdata()

class nwhw():
    def __init__(self,zeitpunkt):          
        self.time = zeitpunkt
        self.folder = zeitpunkt
        self.FNULL = open(os.devnull, 'w')
        vf = self.get_file()
        
        if vf:
            pass
        else:
            self.folder = zeitpunkt - dt.timedelta(days=1)
            vf = self.get_file()
            
        self.valid_file = vf
        self.load_file()
                    
    def get_file(self):
        path = ("/store_linfile11/Vorhersagen/Staumatrix/HWNW_vorhersage"
                "/"+str(self.folder.year)+"/")
        # NWHW
        cmd = ('ssh bm11mos@linwvd10 ls '+ path + '*'
               + self.folder.strftime('_%Y%m%d')  + '*') 

        p = subprocess.Popen(cmd,shell=True,stdout = subprocess.PIPE
                             ,stderr = subprocess.PIPE)
        p.wait()
        filename = p.stdout.read().decode().splitlines()
        p.stdout.close()
        p.stderr.close()
     
        old_dt = np.inf
        valid_file = False
        for x in list(filename):
            delta = ((self.time - 
                 dt.datetime.strptime(x[-12:],'%Y%m%d%H%M')).total_seconds())
            
            if delta > 0 and delta < old_dt:
                valid_file = x
                old_dt = delta
        
        return valid_file
            
    def load_file(self):               
        cmd = 'bm11mos@linwvd10:' + self.valid_file
        p = subprocess.Popen(["scp","-T",cmd,"/home/bm1483/Development"
                              "/manuelle_kurve/scenarios/"],
                             stdout = self.FNULL,
                             stderr = subprocess.STDOUT)
        p.wait()

        shutil.copy2("/home/bm1483/Development/manuelle_kurve/scenarios/"
                     +self.valid_file.split('/')[-1],
                     "/home/bm1483/Development/manuelle_kurve/scenarios"
                     "/HWNW_Vorhersage.txt")
        os.remove("/home/bm1483/Development/manuelle_kurve/scenarios/"
                     +self.valid_file.split('/')[-1])
        
        self.valid_file = ("/home/bm1483/Development/manuelle_kurve/scenarios"
                          "/HWNW_Vorhersage.txt")
        
if __name__ == '__main__':
    x = pegelonline(dt.datetime(2020,2,1,12,55), freq=15)
    y = mos(dt.datetime(2020,2,1,12,16), freq=15)
    z = nwhw(dt.datetime(2020,2,1,12,16))
    
    x.store_data_part()
    y.main()
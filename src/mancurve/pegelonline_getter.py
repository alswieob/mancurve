import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
import warnings
import datetime as dt
import requests

class pegelonline_timeseries():
    def __init__(self,pegel_id, storedir):
        """ Starts Class-Object and directly loads newest data for one location.
        
        :param pegel_id: BSH-Pegelnummer
        :param storedir: Speicherort
        :type pegel_id: string
        :type storedir: string
        """   
        self.pid = pegel_id
        self.a_dir = storedir
        
        self.get_newest_data()
        self.load_data()
                
    def load_data(self):
            """ Loads Pegelonline data. Attaches
            the data to attributes of class instance. 
            """                              
            # Pegelonline                        
            pon_wl,time_array = self.read_pegelonline()           
              
            # Interpolate missing data and smooth it
            # Filter till first NAN by searching in reverse order
            # for first non nan value and cut empty part at the end     
            if len(pon_wl) > 2:                  
                idx = len(pon_wl) - np.isnan(pon_wl[::-1]).argmin()
                time_array = time_array[:idx]  
                pon_wl = pon_wl[:idx]
                
                if np.count_nonzero(np.isnan(pon_wl)) == 0:
                    pon_wl = savgol_filter(pon_wl,151,4)  
                # Wenn weniger als 5% der Daten fehlen, zuerst interpolieren
                elif np.count_nonzero(np.isnan(pon_wl))/len(pon_wl) < 0.5: 
                    pon_wl = pd.Series(pon_wl)
                    # interpolate missing data
                    pon_wl = pon_wl.interpolate(method="linear")
                    # Smooth data
                    pon_wl = savgol_filter(pon_wl,151,4)  
                     
            self.data = pd.DataFrame(pon_wl*0.01, index = time_array,
                                             columns = ['waterlevel'])    
         
    def read_pegelonline(self):
            """ Reads and checks Pegelonline.txt File.
            """          
            with open(self.pon_filename, 'rb') as inp:
                decoded_inp   = [x.decode("ISO-8859-1") for x in inp]
                filtered_inp  = filter(lambda x: '#' in x, decoded_inp)
    
            self.pon_first_date = decoded_inp[0].strip()
            
            if not self.pon_first_date: # Falls keine Gültige Pegelonlinedatei
                return [],[]                    
    
            fil = [i for i in filtered_inp]
            tim_arr = [i.split('#')[0] for i in fil]
            pon_data = [i.split('#')[1].split('\r')[0] for i in fil]
    
            pon_data = [(float(x)) if x.isdigit() else np.nan 
                        for x in pon_data]
    
            # Überprüfung, ob die Daten plausibel sind 
            # Änderunge pro Minute nicht mehr als 50cm
            check_pegel = np.array(pon_data)
            
            # Nur bis zur ersten NAN behalten 
            mask = np.isnan(check_pegel)
            last_valid_idx = np.argmax(~mask[::-1])
            check_pegel=check_pegel[:-last_valid_idx]
            mask = mask[:-last_valid_idx]
            tim_arr = tim_arr[:-last_valid_idx]
            
            # Fill in NaN's with error value...       
            if sum(mask)>0:
                check_pegel[mask] = 10000 #np.interp(np.flatnonzero(mask),
                                   #np.flatnonzero(~mask), check_pegel[~mask])            
            for i in range(2):
                diff = abs(check_pegel[:-1] - check_pegel[1:])
                with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=RuntimeWarning)   
                        # Sprünge ersetzen (+/-1 IDX)
                        where_to_replace  = np.where(diff > 50)
                        where_to_replace2 = np.squeeze(where_to_replace) + 1 
                        where_to_replace3 = np.squeeze(where_to_replace) - 1 
                        if any(where_to_replace3 < 0):                            
                            where_to_replace3 = np.clip(where_to_replace3,
                                                        0,100000)
                            check_pegel[0] = np.nan
                            
                        # Hängengebliebnen pegel ersetzen
                        where_to_replace4 = np.where(diff == 0)                   
                        check_pegel[1:][where_to_replace] = np.nan
                        check_pegel[1:][where_to_replace2] = np.nan
                        check_pegel[1:][where_to_replace3] = np.nan
                        check_pegel[1:][where_to_replace4] = np.nan
            
            pon_data = list(check_pegel)       
            tim_arr    =  np.array(tim_arr) 
            # Umrechnung in UTC
            dates_times = [(dt.datetime.strptime(self.pon_first_date,'%d.%m.%Y') + 
                            dt.timedelta(minutes=x)) + dt.timedelta(hours=-1)
                            for x in range(1,len(tim_arr)+1)]    
            
            inp.close()
            
            return pon_data, dates_times
                    
    def get_newest_data(self):
            """ Download newest numerical oberservation 
            data (pegelonline). Data will be stored in the  "self.a_dir" directory.
            
            Pegelonline data is downloaded from www.pegelonline.wsv.de
            """         
       
            # Pegelonline
            date_array = [(dt.date.today() - 
                           dt.timedelta(days=x)).strftime('%d.%m.%Y')
                           for x in range(1,-1,-1)]
            
            # Names on Pegelonline
            po_elbe = ['HAMBURG+ST.+PAULI','BRUNSB%DCTTEL+MOLE+1',
                       'CUXHAVEN+STEUBENH%D6FT','OSTERIFF+MPM','BROKDORF',
                       'ST%D6R-SPERRWERK+AP','GL%DCCKSTADT','KRAUTSAND','KOLLMAR',
                       'KR%DCCKAU-SPERRWERK+AP','GRAUERORT','PINNAU-SPERRWERK+AP',
                       'STADERSAND','HETLINGEN','L%DCHORT','SCHULAU','CRANZ',
                       'ZOLLENSPIEKER','GEESTHACHT','WEHR+GEESTHACHT+UP',
                       'OTTERNDORF+MPM','BLANKENESE+UF', 'SEEMANNSH%D6FT',
                       'SCH%D6PFSTELLE', 'HAMBURG-HARBURG', 'BUNTHAUS',
                       'OVER', 'ALTENGAMME']
    
            po_nordsee = ['BORKUM+FISCHERBALJE','EIDER-SPERRWERK+AP','BAKE+A',
                          'B%DCSUM','DAGEB%DCLL','HELGOLAND+BINNENHAFEN','HUSUM',
                          'H%D6RNUM','LIST+AUF+SYLT','NORDERNEY+RIFFGAT',
                          'WITTD%DCN','BAKE+C+-+SCHARH%D6RN', 'MITTELGRUND',
                          'BAKE+Z','BORKUM+S%DCDSTRAND','HELGOLAND+S%DCDHAFEN',
                          'LANGEOOG', 'PELLWORM+ANLEGER', 'SPIEKEROOG',
                          'ZEHNERLOCH']
    
            po_ems   = ['LEERORT', 'PAPENBURG','EMDEN+NEUE+SEESCHLEUSE',
                        'EMSH%D6RN', 'DUKEGAT', 'KNOCK','POGUM', 'TERBORG',
                        'WEENER']
    
            po_weser = ['BRAKE','LEUCHTTURM+ALTE+WESER', 'RECHTENFLETH',
                        'BHV+ALTER+LEUCHTTURM', 'DWARSGAT','ROBBENS%DCDSTEERT',
                        'NORDENHAM', 'ELSFLETH', 'FARGE', 'VEGESACK',
                        'OSLEBSHAUSEN']
        
            po_oste  = ['BELUM']
    
            po_jade  = ['WHV+ALTER+VORHAFEN','WANGEROOGE+NORD',
                        'WHV+NEUER+VORHAFEN', 'WANGEROOGE+WEST+','WANGEROOGE+OST',
                        'SCHILLIG', 'MELLUMPLATE', 'HOOKSIELPLATE']
    
            # BSH-Pegel ID-Names
            pegl_ids_elbe    = ['pgl_508P', 'pgl_504P', 'pgl_506P', 'pgl_682P',
                                'pgl_688P', 'pgl_690R','pgl_695P', 'pgl_697P',
                                'pgl_698P', 'pgl_700R', 'pgl_703P', 'pgl_704R',
                                'pgl_709P', 'pgl_711P', 'pgl_712P', 'pgl_714P',
                                'pgl_717P', 'pgl_731P','pgl_732A', 'pgl_732D',
                                'otterndo', 'blankene', 'seemanns', 'schoeste',
                                'hhharbur', 'bunthaus', 'over'    , 'altengam']
                
            pegl_ids_nordsee = ['pgl_101P', 'pgl_664P', 'pgl_677C', 'pgl_505P',
                                'pgl_635P', 'pgl_509A', 'pgl_510P', 'pgl_624P',
                                'pgl_617P', 'pgl_111P', 'pgl_631P', 'bake_c'  , 
                                'mittelgr', 'bake_z',   'borkum_s', 'helgol_s',
                                'langeoog', 'pellworm', 'spiekero', 'zehnerlo']
                
            pegl_ids_ems     = ['pgl_806P', 'pgl_814P', 'pgl_507P', 'emshoern',
                                'dukegat' , 'knock'   , 'pogum'   , 'terborg' ,
                                'weener']
    
            pegl_ids_weser   = ['pgl_743P', 'pgl_734P', 'pgl_741B', 'pgl_103P',
                                'dwarsgat', 'robbenst','nordenha', 'elsfleth',
                                'farge'   , 'vegesack','pgl_502P']
    
            pegl_ids_oste    = ['pgl_683P']
    
            pegl_ids_jade    = ['pgl_512P', 'pgl_754P', 'whv_neuV', 'wange_we',
                                'wange_os', 'schillig', 'mellumpl', 'hooksiel']
            
            def get_pon(url,output,modus='wb'): 
                #if not os.path.isfile(output):
                req = requests.get(url)
                file = open(output, modus)
                for chunk in req.iter_content(100000):
                    file.write(chunk)
                file.close()  
           
            #print('Pegelonlinedownload: '+ self.pid)
            for idx,dates in enumerate(date_array): # Loop over dates     
                # If first date, write modus, else append modus
                if idx == 0:
                    modus = 'wb'
                else:
                    modus = 'ab'   
                               
                # Download elbe data
                if self.pid in pegl_ids_elbe:
                    i = pegl_ids_elbe.index(self.pid)
                    url = ('https://www.pegelonline.wsv.de/webservices/files/'+
                           'Wasserstand+Rohdaten/ELBE/'+po_elbe[i]+'/'+dates+
                           '/down.txt')
                    output = self.a_dir + pegl_ids_elbe[i]+'.txt'
                    get_pon(url,output,modus=modus)
    
                # Download Nordsee data
                elif self.pid in pegl_ids_nordsee:
                    i = pegl_ids_nordsee.index(self.pid)
                    url = ('https://www.pegelonline.wsv.de/webservices/files/'+
                           'Wasserstand+Rohdaten/NORDSEE/'+po_nordsee[i]+
                           '/'+dates+'/down.txt')
                    output = self.a_dir+pegl_ids_nordsee[i]+'.txt'
                    get_pon(url,output,modus=modus)
          
                # Download Ems data
                elif self.pid in pegl_ids_ems:
                    i = pegl_ids_ems.index(self.pid)
                    url = ('https://www.pegelonline.wsv.de/webservices/files/'+
                           'Wasserstand+Rohdaten/EMS/'+po_ems[i]+'/'+dates+'/down.txt')
                    output = self.a_dir+pegl_ids_ems[i]+'.txt'
                    get_pon(url,output,modus=modus)
        
                # Download Weser data
                elif self.pid in pegl_ids_weser:
                    i = pegl_ids_weser.index(self.pid)
                    url = ('https://www.pegelonline.wsv.de/webservices/files/'+
                           'Wasserstand+Rohdaten/WESER/'+po_weser[i]+'/'+
                           dates+'/down.txt')
                    output = self.a_dir+pegl_ids_weser[i]+'.txt'
                    get_pon(url,output,modus=modus)      
        
                # Download Oste data
                elif self.pid in pegl_ids_oste:
                    i = pegl_ids_oste.index(self.pid)
                    url = ('https://www.pegelonline.wsv.de/webservices/files/'+
                           'Wasserstand+Rohdaten/OSTE/'+po_oste[i]+'/'+dates+
                           '/down.txt')
                    output = self.a_dir+pegl_ids_oste[i]+'.txt'                        
                    get_pon(url,output,modus=modus)
                           
                # Download Jade data
                elif self.pid in pegl_ids_jade:
                    i = pegl_ids_jade.index(self.pid)
                    url = ('https://www.pegelonline.wsv.de/webservices/files/'+
                           'Wasserstand+Rohdaten/JADE/'+po_jade[i]+'/'
                           +dates+'/down.txt')
                    output = self.a_dir+pegl_ids_jade[i]+'.txt'
                    get_pon(url,output,modus=modus)  
                                  
            self.pon_filename = output                 
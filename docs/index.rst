.. mancur documentation master file, created by
   sphinx-quickstart on Mon Mar 30 09:17:27 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

    
.. toctree::
   :maxdepth: 3
   

Mancurve
====================

Dokumentation für die Kombination aus
 Pegelonline / MOS | Astronomie und manueller Vorhersage.


Hauptklasse mit zentralen Methoden
==================================

.. autoclass:: mancurve.core.combine_curves
    :members:
    
Klasse für Pegelonline Datenbezug
========================================

.. automethod:: mancurve.pegelonline_getter.pegelonline_timeseries.__init__

.. autoclass::  mancurve.pegelonline_getter.pegelonline_timeseries
     :members:
     
Klassen für Szenariorechnung
======================================== 

.. automethod:: mancurve.data_feeder.pegelonline.__init__
    
.. autoclass::  mancurve.data_feeder.pegelonline
     :members:
     
.. automethod:: mancurve.data_feeder.mos.__init__
 
.. autoclass::  mancurve.data_feeder.mos
     :members:     
     
.. automethod:: mancurve.data_feeder.nwhw.__init__
 
.. autoclass::  mancurve.data_feeder.nwhw
     :members:
     
Paket zur Validation (Plots & Statistik)
======================================== 
.. automodule:: mancurve.validation
     :members:
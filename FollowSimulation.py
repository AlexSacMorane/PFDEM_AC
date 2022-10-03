# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

Following live simulation.
A debuging script. Must be run from a GUI (Spyder for example)
"""

#-----------------------------------------------------------------------------
# Librairies
#-----------------------------------------------------------------------------

import pickle
import matplotlib.pyplot as plt

#-----------------------------------------------------------------------------
# Load data
#-----------------------------------------------------------------------------

i_PF = 2
name_file = 'Debug/DEM_ite/PF_'+str(i_PF)+'/save_tempo' #data to follow
file = open(name_file,'rb')
dict = pickle.load(file)
file.close

#-----------------------------------------------------------------------------
# Work
#-----------------------------------------------------------------------------

#data
Ecin_stop = dict['E_cin_stop']
dk0_stop = dict['dk0_stop']
dy_box_max_stop = dict['dy_box_max_stop']
n_window_stop = dict['n_window_stop']
Ecin_simu = dict['E_cin'][1:]
Force_simu = dict['Force'][1:]
k0_min_simu = dict['k0_xmin_tracker']
k0_max_simu = dict['k0_xmax_tracker']
ymax_simu = dict['y_box_max'][:-1]

#iteration on times to verify if stop criteria are reached
Ecin_stop_L=[]
k0_min_stop_L=[]
k0_max_stop_L=[]
ymax_stop_L=[]
for i in range(n_window_stop-1,len(Ecin_simu)):
    Ecin = Ecin_simu[i]
    k0_min_L = k0_min_simu[i+1-n_window_stop:i+1]
    k0_max_L = k0_max_simu[i+1-n_window_stop:i+1]
    ymax_L = ymax_simu[i+1-n_window_stop:i+1]

    if Ecin<Ecin_stop:
        Ecin_stop_L.append(1)
    else:
        Ecin_stop_L.append(0)
    if max(k0_min_L)-min(k0_min_L)<dk0_stop:
        k0_min_stop_L.append(3)
    else:
        k0_min_stop_L.append(2)
    if max(k0_max_L)-min(k0_max_L)<dk0_stop:
        k0_max_stop_L.append(5)
    else:
        k0_max_stop_L.append(4)
    if max(ymax_L)-min(ymax_L)<dy_box_max_stop:
        ymax_stop_L.append(7)
    else:
        ymax_stop_L.append(6)

#-----------------------------------------------------------------------------
# Plot
#-----------------------------------------------------------------------------

plt.close(1)
plt.figure(1,figsize=(16,9))
plt.plot(range(n_window_stop-1,len(Ecin_simu)),ymax_stop_L,label='ymax')
plt.plot(range(n_window_stop-1,len(Ecin_simu)),k0_max_stop_L,label='k0_max')
plt.plot(range(n_window_stop-1,len(Ecin_simu)),k0_min_stop_L,label='k0_min')
plt.plot(range(n_window_stop-1,len(Ecin_simu)),Ecin_stop_L,label='Ecin')
ax = plt.gca()
ax.set_yticks([0,1,2,3,4,5,6,7])
ax.set_yticklabels(['False','True','False','True','False','True','False','True'])
plt.legend()

#-----------------------------------------------------------------------------

plt.close(2)
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(16,9),num=2)

ax1.set_title('Mean kinetic energy (e-12 J)')
ax1.plot(dict['E_cin'])
ax1.plot([0, len(Ecin_simu)-1],[Ecin_stop, Ecin_stop],'r')

ax2.set_title('Mean force applied (µN)')
ax2.plot(Force_simu)

ax3.set_title('k0s (-)')
ax3.plot(k0_min_simu,label='xmin')
ax3.plot(k0_max_simu,label='xmax')
ax3.legend()

ax4.set_title('About the upper plate')
ax4.plot(ymax_simu, color = 'blue')
ax4.set_ylabel('Coordinate y (µm)', color = 'blue')
ax4.tick_params(axis ='y', labelcolor = 'blue')
ax4a = ax4.twinx()
ax4a.plot([100, len(dict['F_on_ymax'])-1],[dict['F_on_ymax'][-1], dict['F_on_ymax'][-1]], color = 'red')
ax4a.plot(range(100,len(dict['F_on_ymax'])),dict['F_on_ymax'][100:], color = 'orange')
ax4a.set_ylabel('Force applied (µN)', color = 'orange')
ax4a.tick_params(axis ='y', labelcolor = 'orange')

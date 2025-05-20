import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.gridspec import GridSpec

matplotlib.rcParams.update({'errorbar.capsize': 2})
plt.close('all')

def array_prep(arr, arr_conf):
    arr = np.asarray(arr , dtype=float)             
    arr_conf = np.asarray(arr_conf, dtype=float)
    arr_conf[:,0] = arr_conf[:,0] - arr
    arr_conf[:,1] = arr - arr_conf[:,1] 
    arr_conf[arr_conf<0] = 0
    arr_conf = arr_conf.transpose()
    return arr, arr_conf

reversal_number = 100
max_gap = 1000
min_rem_yr = 100

max_gap = max_gap/1e6

lost_magnetozones = []
lost_magnetozones_conf = []

lost_magnetozones_thikness = []
lost_magnetozones_thikness_conf = []

pr_lost_magnetozones = []
pr_lost_magnetozones_conf = []

lost_change_zones = []
lost_change_zones_conf = []

pr_lost_change_zones = []
pr_lost_change_zones_conf = []

lost_change_zones_thikness = []
lost_change_zones_thikness_conf = []

gap_percent_list = [10,20,30,40,50,60,70,80,90,95,98.5]

for file in os.listdir("."):
    if file.endswith(".txt"):
        f_split = file.split()
        if(f_split[2] == "reversals_" + str(reversal_number) and f_split[5] == "gap_" + str(max_gap) and f_split[9] ==  str(min_rem_yr) + '.0.txt') :
            
            with open(os.path.join(".", file), 'r') as file:                
                lines = [line.rstrip() for line in file]

                lost_magnetozones.append(lines[10].replace(',', ' ').split()[-1])                
                lost_magnetozones_conf.append([lines[9].replace(',', ' ').split()[4], lines[9].replace(',', ' ').split()[3]])
                
                pr_lost_magnetozones.append(lines[13].replace(',', ' ').split()[-1])    
                pr_lost_magnetozones_conf.append([lines[12].replace(',', ' ').split()[4], lines[12].replace(',', ' ').split()[3]])
                
                lost_magnetozones_thikness.append(lines[16].replace(',', ' ').split()[-1])    
                lost_magnetozones_thikness_conf.append([lines[15].replace(',', ' ').split()[4], lines[15].replace(',', ' ').split()[3]])
                
                lost_change_zones.append(lines[22].replace(',', ' ').split()[-1])   
                lost_change_zones_conf.append([lines[21].replace(',', ' ').split()[4], lines[21].replace(',', ' ').split()[3]])

                lost_change_zones_thikness.append(lines[19].replace(',', ' ').split()[-1])   
                lost_change_zones_thikness_conf.append([lines[18].replace(',', ' ').split()[4], lines[18].replace(',', ' ').split()[3]])

                pr_lost_change_zones.append(lines[25].replace(',', ' ').split()[-1]) 
                pr_lost_change_zones_conf.append([lines[24].replace(',', ' ').split()[4], lines[24].replace(',', ' ').split()[3]])


lost_magnetozones, lost_magnetozones_conf = array_prep(lost_magnetozones, lost_magnetozones_conf)              
pr_lost_magnetozones, pr_lost_magnetozones_conf = array_prep(pr_lost_magnetozones, pr_lost_magnetozones_conf)
lost_change_zones, lost_change_zones_conf = array_prep(lost_change_zones, lost_change_zones_conf)      
lost_change_zones_thikness, lost_change_zones_thikness_conf = array_prep(lost_change_zones_thikness, lost_change_zones_thikness_conf)
pr_lost_change_zones, pr_lost_change_zones_conf = array_prep(pr_lost_change_zones, pr_lost_change_zones_conf)
lost_magnetozones_thikness, lost_magnetozones_thikness_conf = array_prep(lost_magnetozones_thikness, lost_magnetozones_thikness_conf)

#figure, axis = plt.subplots(4, 2, figsize=(14, 12),  constrained_layout=True)

figure = plt.figure(figsize=(14, 16), constrained_layout=True)
gs = GridSpec(4, 2, figure=figure)

ax1 = figure.add_subplot(gs[0, 0]) 
ax2 = figure.add_subplot(gs[0, 1])  
ax3 = figure.add_subplot(gs[1, 0])  
ax4 = figure.add_subplot(gs[1, 1]) 
ax5 = figure.add_subplot(gs[2, 0])  
ax6 = figure.add_subplot(gs[2, 1]) 
ax7 = figure.add_subplot(gs[3, :]) 

axis = np.array([[ax1, ax2], [ax3, ax4], [ax5, ax6], [ax7, None]])

if (reversal_number == 20): 
    figure.suptitle(str(int(reversal_number/5))+' инверсии за миллион лет. Длительность диастем до ' + str(int(max_gap*1e6)) + ' лет. Минимальный интервал ' + str(min_rem_yr) + ' лет.', fontsize=16)

if (reversal_number == 100): 
    figure.suptitle(str(int(reversal_number/5))+' инверсий за миллион лет. Длительность диастем до ' + str(int(max_gap*1e6)) + ' лет. Минимальный интервал ' + str(min_rem_yr) + ' лет.', fontsize=16)


#figure.tight_layout()

axis[0,0].errorbar(gap_percent_list, lost_magnetozones, yerr=lost_magnetozones_conf, fmt='o')
axis[0,0].set_ylim(0,100)
axis[0,0].set_xticks([20, 40, 60, 80])
axis[0,0].set_title('Потерянные магнитозоны')
axis[0,0].set_ylabel('Процент потерь')
axis[0,0].grid(True)

axis[0,1].errorbar(gap_percent_list, pr_lost_magnetozones, yerr=pr_lost_magnetozones_conf, fmt='o', color = '#ff7f0e')
axis[0,1].set_ylim(0,100)
axis[0,1].set_xticks([20, 40, 60, 80])
axis[0,1].set_title('Фактически потерянные магнитозоны')
axis[0,1].grid(True)

axis[1,0].errorbar(gap_percent_list, lost_change_zones, yerr=lost_change_zones_conf, fmt='o', color = '#2ca02c')
axis[1,0].set_ylim(0,100)
axis[1,0].set_xticks([20, 40, 60, 80])
axis[1,0].set_title('Потерянные переходные зоны')
axis[1,0].set_ylabel('Процент потерь')
axis[1,0].grid(True)

axis[1,1].errorbar(gap_percent_list, pr_lost_change_zones, yerr=pr_lost_change_zones_conf, fmt='o', color = '#d62728')
axis[1,1].set_ylim(0,100)
axis[1,1].set_xticks([20, 40, 60, 80])
axis[1,1].set_title('Фактически потерянные переходные зоны')
axis[1,1].grid(True)

axis[2,0].errorbar(gap_percent_list, lost_change_zones_thikness, yerr=lost_change_zones_thikness_conf, fmt='o', color = '#9467bd')
axis[2,0].set_ylim(0,100)
axis[2,0].set_xticks([20, 40, 60, 80])
#axis[2,0].set_xlabel('Процент перерывов')
axis[2,0].set_title('Потеря мощности переходных зон')
axis[2,0].set_ylabel('Процент потерь')
axis[2,0].grid(True)

axis[2,1].errorbar(gap_percent_list, lost_magnetozones_thikness, yerr=lost_magnetozones_thikness_conf, fmt='o', color = '#8c564b')
axis[2,1].set_ylim(0,100)
axis[2,1].set_xticks([20, 40, 60, 80])
#axis[2,0].set_xlabel('Процент перерывов')
axis[2,1].set_title('Потеря мощности магнитозон')
#axis[2,1].set_ylabel('Процент потерь')
axis[2,1].grid(True)


axis[3,0].errorbar(gap_percent_list, lost_magnetozones, yerr=lost_magnetozones_conf, fmt='o', label = 'Потерянные магнитозоны')
axis[3,0].errorbar(gap_percent_list, pr_lost_magnetozones, yerr=pr_lost_magnetozones_conf, fmt='o', label  = 'Фактически потерянные магнитозоны')
axis[3,0].errorbar(gap_percent_list, lost_change_zones, yerr=lost_change_zones_conf, fmt='o', label = 'Потерянные переходные зоны')
axis[3,0].errorbar(gap_percent_list, pr_lost_change_zones, yerr=pr_lost_change_zones_conf, fmt='o', label = 'Фактически потерянные переходные зоны')
axis[3,0].errorbar(gap_percent_list, lost_change_zones_thikness, yerr=lost_change_zones_thikness_conf, fmt='o', label = 'Потеря мощности переходных зон')
axis[3,0].errorbar(gap_percent_list, lost_magnetozones_thikness, yerr=lost_magnetozones_thikness_conf, fmt='o', label = 'Потеря мощности магнитозон')
axis[3,0].set_ylim(0,100)
axis[3,0].set_xticks([20, 40, 60, 80])
axis[3,0].set_xlabel('Процент перерывов')
axis[3,0].legend(loc = 'upper left')
axis[3,0].set_title('Сравнение')
axis[3,0].grid(True)

plt.show()      

if (reversal_number == 20): 
    plt.savefig("./img/" + str(int(reversal_number/5))+' инверсии за миллион лет. Длительность диастем до ' + str(int(max_gap*1e6)) + ' лет. Минимальный интервал ' + str(min_rem_yr) + ' лет.png')   

if (reversal_number == 100): 
    plt.savefig("./img/" + str(int(reversal_number/5))+' инверсий за миллион лет. Длительность диастем до ' + str(int(max_gap*1e6)) + ' лет. Минимальный интервал ' + str(min_rem_yr) + ' лет.png')   
    
       




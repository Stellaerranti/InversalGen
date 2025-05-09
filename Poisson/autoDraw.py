import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
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
max_gap = 0.0005
min_rem_yr = 20

#max_gap = max_gap/1e6

lost_magnetozones = []
lost_magnetozones_conf = []

pr_lost_magnetozones = []
pr_lost_magnetozones_conf = []

lost_change_zones = []
lost_change_zones_conf = []

pr_lost_change_zones = []
pr_lost_change_zones_conf = []

lost_change_zones_thikness = []
lost_change_zones_thikness_conf = []

gap_percent_list = [10,20,30,40,50,60,70,80,90]

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
                
                lost_change_zones.append(lines[19].replace(',', ' ').split()[-1])   
                lost_change_zones_conf.append([lines[18].replace(',', ' ').split()[4], lines[18].replace(',', ' ').split()[3]])

                lost_change_zones_thikness.append(lines[16].replace(',', ' ').split()[-1])   
                lost_change_zones_thikness_conf.append([lines[15].replace(',', ' ').split()[4], lines[15].replace(',', ' ').split()[3]])

                pr_lost_change_zones.append(lines[22].replace(',', ' ').split()[-1]) 
                pr_lost_change_zones_conf.append([lines[21].replace(',', ' ').split()[4], lines[21].replace(',', ' ').split()[3]])


lost_magnetozones, lost_magnetozones_conf = array_prep(lost_magnetozones, lost_magnetozones_conf)              
pr_lost_magnetozones, pr_lost_magnetozones_conf = array_prep(pr_lost_magnetozones, pr_lost_magnetozones_conf)
lost_change_zones, lost_change_zones_conf = array_prep(lost_change_zones, lost_change_zones_conf)      
lost_change_zones_thikness, lost_change_zones_thikness_conf = array_prep(lost_change_zones_thikness, lost_change_zones_thikness_conf)
pr_lost_change_zones, pr_lost_change_zones_conf = array_prep(pr_lost_change_zones, pr_lost_change_zones_conf)

figure, axis = plt.subplots(3, 2, figsize=(14, 12),  constrained_layout=True)

figure.suptitle(str(int(reversal_number/5))+' инверсий за миллион лет. Средняя длина диастем ' + str(int(max_gap*1e6)) + ' лет. Минимальный интервал ' + str(min_rem_yr) + ' лет.', fontsize=16)


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
axis[2,0].set_xlabel('Процент перерывов')
axis[2,0].set_title('Потеря мощности переходных зон')
axis[2,0].set_ylabel('Процент потерь')
axis[2,0].grid(True)


axis[2,1].errorbar(gap_percent_list, lost_magnetozones, yerr=lost_magnetozones_conf, fmt='o', label = 'Потерянные магнитозоны')
axis[2,1].errorbar(gap_percent_list, pr_lost_magnetozones, yerr=pr_lost_magnetozones_conf, fmt='o', label  = 'Фактически потерянные магнитозоны')
axis[2,1].errorbar(gap_percent_list, lost_change_zones, yerr=lost_change_zones_conf, fmt='o', label = 'Потерянные переходные зоны')
axis[2,1].errorbar(gap_percent_list, pr_lost_change_zones, yerr=pr_lost_change_zones_conf, fmt='o', label = 'Фактически потерянные переходные зоны')
axis[2,1].errorbar(gap_percent_list, lost_change_zones_thikness, yerr=lost_change_zones_thikness_conf, fmt='o', label = 'Потеря мощности переходных зон')
axis[2,1].set_ylim(0,100)
axis[2,1].set_xticks([20, 40, 60, 80])
axis[2,1].set_xlabel('Процент перерывов')
axis[2,1].legend(loc = 'upper left')
axis[2,1].set_title('Сравнение')
axis[2,1].grid(True)

plt.show()      

plt.savefig("./img/" + str(int(reversal_number/5))+' инверсий за миллион лет. Средняя длина диастем ' + str(int(max_gap*1e6)) + ' лет. Минимальный интервал ' + str(min_rem_yr) + ' лет.png')   





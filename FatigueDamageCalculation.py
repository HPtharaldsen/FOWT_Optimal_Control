# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 11:01:24 2022

@author: seragela
"""

#Importing module
import struct
import numpy as np
import matplotlib.pyplot as plt
import rainflow as rfc

# Base path
base_path = r"c:\Users\thara\SIMA Workspaces\TMR03\IEA15MW_VolturnUS_S\Part3\Part3_"

# Function to unpack binary results
def unpackBinaryResults(fileName, Nt, Nchan):
    with open(fileName, mode='rb') as file:
        fileContent = file.read()
        
    numRecords = int((len(fileContent) / 4))
    unpacked = struct.unpack("f" * numRecords, fileContent)
    if Nt > 1:
        Nchan = len(np.asarray(unpacked)) / Nt
    else:
        Nt = len(np.asarray(unpacked)) / Nchan
    
    A = np.reshape(np.asarray(unpacked), (int(Nchan), int(Nt)), order='F')     
    return A 

def process_part(part_number, NtRiflexForce):
    part_path = f"{base_path}{part_number}\\"
    ForcefileName = part_path + 'sima_elmfor.bin'
    A = unpackBinaryResults(ForcefileName, NtRiflexForce, 0)

    TowerBaseAxial   = A[2]  #Axial force
    TowerBaseBMY = A[4]  # Bending moment in local y
    TowerBaseBMZ = A[6]  # Bending moment in local z

    return TowerBaseAxial, TowerBaseBMY, TowerBaseBMZ

# Define the mapping of parts to wind speeds
wind_speed_mapping = {
    4: ["1" , "2", "5", "6"], # 101, 102, 105, 106
    6: ["3", "4", "7", "8"], # 103, 104, 107, 108
    8: ["9", "10", "11", "12", # 109, 110, 201, 202
         "13", "14", "15", "16"], # 203, 204, 111, 112
    10: ["17", "18", "21", "22"], # 113, 114, 117, 118
    12: ["19", "20", "23", "24"],  # 115, 116, 119, 120
    14: ["25", "26"], # 121, 122
    18: ["27", "28", "29", # 123, 124, 125
         "30", "31", "32"], # 126, 127, 128 
    20: ["33", "34"], # 129, 130
    22: ["35", "36"], # 131, 132
    24: ["37", "38"] # 133, 134
}


# Number of stored time steps in the force file
NtRiflexForce = 5000  # adjust this as needed

# Process all parts and condition by wind speed
bending_moments_by_wind_speed = {ws: {'Axial': [], 'BMY': [], 'BMZ': []} for ws in wind_speed_mapping}
for wind_speed, parts in wind_speed_mapping.items():
    for part in parts:
        print(f"Processing Part: {part} for Wind Speed: {wind_speed} m/s")
        Axial, BMY, BMZ = process_part(part, NtRiflexForce)
        bending_moments_by_wind_speed[wind_speed]['Axial'].append(Axial)
        bending_moments_by_wind_speed[wind_speed]['BMY'].append(BMY)
        bending_moments_by_wind_speed[wind_speed]['BMZ'].append(BMZ)

# Defining FatigueDamage calculation
def FatigueDamage(stress_history,thk,thk_ref,k,K1,beta1,stress_lim=52.639,K2=1/10**15.606,beta2=5):
    cc = rfc.count_cycles(stress_history)
    sk = [c[0] for c in cc] # stress range
    n = [c[1] for c in cc] # cycle count
    
    Ns = np.zeros(len(sk)) #initialize damage
    
    for i,s in enumerate(sk):
        if s>stress_lim:
            beta = beta1; K = K1;
        else:
            beta = beta2; K = K2;
        
        Ns[i] = 1/K*(s*(thk/thk_ref)**0.2)**(-beta)
        
    FD = np.sum(n/Ns)
    
    return FD   

# Set-up fatigue imput
thk     = 82.95
thk_ref = 25
r       = 10/2
A       = np.pi * (r**2 - (r-thk)**2)
I_yy    = np.pi/4 * (r**4 - (r-thk)**4)
I_zz    = np.pi/4 * (r**4 - (r-thk)**4)
theta   = np.radians(270)

k = 0.20
K1 = 1/(10**12.164)
beta1 = 3

# Calculating stress history
stress_history = {ws: {'Axial': None, 'BMY': None, 'BMZ': None} for ws in wind_speed_mapping}

for wind_speed, moments in bending_moments_by_wind_speed.items():
    stress_history[wind_speed] = (Axial/A + BMY/I_yy*r*np.sin(theta) + BMZ/I_zz*r*np.cos(theta))*10e-6
        

fatigue_damage_1h_per_seed = {}
wave_seed_to_number = {"1": 101, "2": 102, "3": 103, "4": 104, "5": 105, "6": 106, 
                       "7": 107, "8": 108, "9": 109, "10": 110, "11": 201, "12": 202, 
                       "13": 203, "14": 204, "15": 111, "16": 112, "17": 113, "18": 114, 
                       "19": 115, "20": 116, "21": 117, "22": 118, "23": 119, "24": 120,
                       "25": 121, "26": 122, "27": 123, "28": 124, "29": 125, "30": 126, 
                       "31": 127, "32": 128, "33": 129, "34": 130, "35": 131, "36": 132, 
                       "37": 133, "38": 134}



for wind_speed, parts in wind_speed_mapping.items():
    for part in parts:
        wave_seed_number = wave_seed_to_number[part]
        Axial, BMY, BMZ = process_part(part, NtRiflexForce)
        stress_history = (Axial/A + BMY/I_yy*r*np.sin(theta) + BMZ/I_zz*r*np.cos(theta))*10e-6
        fatigue_damage_1h_per_seed[wave_seed_number] = FatigueDamage(stress_history,thk,thk_ref,k,K1,beta1,stress_lim=52.639,K2=1/10**15.606,beta2=5)

print(list(fatigue_damage_1h_per_seed.keys()))
# Creating Figure for 1 hour fatigue damage per wave seed
plt.figure()
plt.title('1h Fatigue Damage per Wave Seed')
plt.plot(list(fatigue_damage_1h_per_seed.keys()), list(fatigue_damage_1h_per_seed.values()), marker='o', linestyle='-', color='blue', label='1h Fatigue per Seed')
plt.xlabel('Wave Seed Number')
plt.xticks(range(len(list(fatigue_damage_1h_per_seed.keys()))), ['101','102','103','104','105','106','107','108','109','110',
                                                                 '201', '202', '203', '204',
                                                                 '111','112','113','114','115','116','117','118',
                                                                 '119', '120', '121', '122', '123', '124', '125', '126', '127', '128', '129', '130',
                                                                 '131', '132', '133', '134'])
plt.ylabel('Fatigue Damage (-)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Define the conditions of wave seeds
conditions = {
    "1": [101, 102],
    "2": [103, 104],
    "3": [105, 106],
    "4": [107, 108],
    "5": [109, 110, 201, 202, 203, 204],
    "6": [111, 112],
    "7": [113, 114],
    "8": [115, 116],
    "9": [117, 118],
    "10": [119, 120],
    "11": [121, 122],
    "12": [123, 124],
    "13": [125, 126],
    "14": [127, 128],
    "15": [129, 130],
    "16": [131, 132],
    "17": [133, 134]
}

# Initialize a dictionary to store mean fatigue damage per condition
mean_fatigue_damage_per_condition = {condition: 0 for condition in conditions}

# Iterate over wind_speed_mapping and calculate mean fatigue damage for each condition
for wind_speed, parts in wind_speed_mapping.items():
    for part in parts:
        wave_seed_number = wave_seed_to_number[part]
        Axial, BMY, BMZ = process_part(part, NtRiflexForce)
        stress_history = (Axial/A + BMY/I_yy*r*np.sin(theta) + BMZ/I_zz*r*np.cos(theta))*10e-6
        fatigue_damage = FatigueDamage(stress_history, thk, thk_ref, k, K1, beta1, stress_lim=52.639, K2=1/10**15.606, beta2=5)
        
        # Determine which condition this wave seed belongs to
        condition = None
        for key, value in conditions.items():
            if wave_seed_number in value:
                condition = key
                break
        
        # Update the mean fatigue damage for the condition
        if condition:
            mean_fatigue_damage_per_condition[condition] += fatigue_damage

# Calculate the mean fatigue damage for each condition
for condition, total_damage in mean_fatigue_damage_per_condition.items():
    num_seeds_in_condition = len(conditions[condition])
    mean_fatigue_damage_per_condition[condition] = total_damage / num_seeds_in_condition

# Creating Figure for mean 1 hour fatigue damage per condition
plt.figure()
plt.title('Mean 1h Fatigue Damage per condition of Wave Seeds')
plt.bar(mean_fatigue_damage_per_condition.keys(), mean_fatigue_damage_per_condition.values(), color='blue', label='Mean 1h Fatigue per condition')
plt.xlabel('condition')
plt.ylabel('Mean Fatigue Damage (-)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

mean_fatigue_damage_1h = {}
for wind_speed, moments in bending_moments_by_wind_speed.items():
    total_fatigue_damage = 0
    count = 0
    for Axial, BMY, BMZ in zip(moments['Axial'], moments['BMY'], moments['BMZ']):
        stress_history = (Axial/A + BMY/I_yy*r*np.sin(theta) + BMZ/I_zz*r*np.cos(theta))*10e-6
        total_fatigue_damage += FatigueDamage(stress_history,thk,thk_ref,k,K1,beta1,stress_lim=52.639,K2=1/10**15.606,beta2=5)
        count += 1
    mean_fatigue_damage_1h[wind_speed] = total_fatigue_damage / count

# Creating Figure for Mean 1 hour fatigue damage
plt.figure()
plt.title('Mean 1h Fatigue Damage')
plt.plot(list(mean_fatigue_damage_1h.keys()), list(mean_fatigue_damage_1h.values()), marker='o', linestyle='-', color='blue', label='Mean 1h Fatigue')
plt.xlabel('Wind Speed (m/s)')
plt.ylabel('Mean Fatigue Damage (-)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

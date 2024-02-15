# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 11:01:24 2022

@author: seragela
"""

# Importing modules
import struct
import numpy as np
import matplotlib.pyplot as plt
import rainflow as rfc

# Function to unpack binary results from a single file
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

# Function to process the single file
def process_file(filePath, NtRiflexForce):
    A = unpackBinaryResults(filePath, NtRiflexForce, 0)

    TowerBaseAxial = A[2]  # Axial force
    TowerBaseBMY = A[4]    # Bending moment in local y
    TowerBaseBMZ = A[6]    # Bending moment in local z

    return TowerBaseAxial, TowerBaseBMY, TowerBaseBMZ

# Path to the single input file
inputFilePath = r"c:\Path\To\Your\File\sima_elmfor.bin" # Adjust this path to your single input file

# Number of stored time steps in the force file
NtRiflexForce = 5000  # Adjust this as needed

# Process the single file
print("Processing the input file...")
Axial, BMY, BMZ = process_file(inputFilePath, NtRiflexForce)
print("Processing completed.")

# Defining FatigueDamage calculation
def FatigueDamage(stress_history, thk, thk_ref, k, K1, beta1, stress_lim=52.639, K2=1/10**15.606, beta2=5):
    cc = rfc.count_cycles(stress_history)
    sk = [c[0] for c in cc] # stress range
    n = [c[1] for c in cc] # cycle count
    
    Ns = np.zeros(len(sk)) # initialize damage
    
    for i, s in enumerate(sk):
        if s > stress_lim:
            beta = beta1
            K = K1
        else:
            beta = beta2
            K = K2
        
        Ns[i] = 1/K * (s * (thk / thk_ref) ** 0.2) ** (-beta)
        
    FD = np.sum(n / Ns)
    
    return FD   

# Example calculation for Fatigue Damage using the bending moments from one file for a wind turbine tower
thk = 82.95
thk_ref = 25
r = 10 / 2
A = np.pi * (r**2 - (r - thk)**2)
I_yy = np.pi / 4 * (r**4 - (r - thk)**4)
I_zz = np.pi / 4 * (r**4 - (r - thk)**4)
theta = np.radians(270)

k = 0.20
K1 = 1 / (10**12.164)
beta1 = 3

# Calculating stress history from processed file data
stress_history = (Axial/A + BMY/I_yy*r*np.sin(theta) + BMZ/I_zz*r*np.cos(theta)) * 10e-6

# Calculating Fatigue Damage
fatigue_damage = FatigueDamage(stress_history, thk, thk_ref, k, K1, beta1)
print(f"Fatigue Damage: {fatigue_damage}")

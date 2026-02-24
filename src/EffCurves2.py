import argparse 
import matplotlib.pyplot as plt
import numpy as np
import glob
import os

def eff_curve(S_1, S_2, S_3, D, T, name):


    if name == "PMT01":

        DF = S_2 * S_3 * 1e-11 

        err_DF = DF * np.sqrt(1/S_2 + 1/S_3)

    if name == "PMT02":

        DF = S_1 * S_3 * 1e-11

        err_DF = DF * np.sqrt(1/S_1 + 1/S_3) 
    
    if name == "PMT03":

        DF = S_2 * S_1 * 1e-11

        err_DF = DF * np.sqrt(1/S_2 + 1/S_1)

    eff = T/(D-DF) 

    err_eff = np.sqrt(eff*(1-eff)/(D-DF))

    return eff * 100 , err_eff * 100, DF, err_DF

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "This program accepts a path to a folder, it estimates efficency and associated error, and plots efficency as a function of HV")

    parser.add_argument("folder", type = str, help = "Path.txt")

    args = parser.parse_args()

    folder = args.folder

    files = glob.glob(os.path.join(folder, "PMT*.txt"))

    fig, ax = plt.subplots()
    
    ax.set_xlabel('HV [V]')
    
    ax.set_ylabel('T/D [%]')
    
    ax.grid(True)

    fig_n, ax_n = plt.subplots()
    
    ax_n.set_xlabel('HV [V]')
    
    ax_n.set_ylabel('DF/D [%]')
    
    ax_n.grid(True)

    labels = []

    for f in files:

        title = os.path.splitext(os.path.basename(f))[0]
        
        HV, S_1, S_2, S_3, D, T = np.loadtxt(f, unpack = True)

        eff, eff_err, DF, err_DF = eff_curve(S_1, S_2, S_3, D, T, title)

        r = DF / D

        err_r = r * np.sqrt((err_DF/DF)**2 + 1/D) 
        
        ax.errorbar(HV, eff, eff_err, fmt = '.' , capsize = 3)

        ax_n.errorbar(HV, r * 100, err_r * 100, fmt = '.' , capsize = 3)
    
        labels += [title]
        
    ax.legend(labels)

    ax_n.legend(labels)
    
    folder_name = os.path.basename(os.path.normpath(folder))
    
    ax.set_title(f'Curve efficienza - {folder_name}')

    ax_n.set_title(f'Contaminazione delle doppie - {folder_name}')
    
    plt.show()

    
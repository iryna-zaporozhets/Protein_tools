import numpy as np
import os
import re
import time
import mdtraj as md
import md_nmr2 as nmr2

# Test 1. Average_mode of analysis 
def test_rdc_average():
    print("")
    print("Testing function calculate_rdc_large with mode='average'")
    path_traj = 'trajectory.xtc'
    path_top  = 'topology.pdb'
    RDC_inp_file = "experimental_data.txt" # Experimental data - allignment A28               
    exp_rdc,D_av = nmr2.calculate_rdc_large(path_traj,
                                                   path_top,
                                                   RDC_inp_file,
                                                   minimize_rmsd=False, # Always use the trajectory where RMSD is already minimized
                                                   mode='average'
                                                   )     
    
    reference = np.loadtxt("reference_calculated_average.txt")
    difference  = np.subtract(D_av, reference)
    print("Maximum difference between calculated and average value:  %f" %np.max(difference))
    print("Mean difference between calculated and average value:     %f" %np.mean(difference))
    return()


def test_rdc_full():
    print("")
    print("Testing function calculate_rdc_large with mode='full'")
    path_traj = 'trajectory.xtc'
    path_top  = 'topology.pdb'
    RDC_inp_file = "experimental_data.txt" # Experimental data - allignment A28               
    exp_rdc,D_av = nmr2.calculate_rdc_large(path_traj,
                                                   path_top,
                                                   RDC_inp_file,
                                                   minimize_rmsd=False, # Always use the trajectory where RMSD is already minimized
                                                   mode='full'
                                                   )     
    
    reference = np.loadtxt("reference_calculated_full.txt")
    difference  = np.subtract(D_av, reference)
    print("Maximum difference between calculated and average value:  %f" %np.max(difference))
    print("Mean difference between calculated and average value:     %f" %np.mean(difference))
    return()

def test_rdc_full_amide():
    print(" ")
    print("Testing function calculate_rdc_amide_large with mode='full'")
    path_traj = 'trajectory.xtc'
    path_top  = 'topology.pdb'
    RDC_inp_file = "experimental_data.txt" # Experimental data - allignment A28               
    exp_rdc,D_av = nmr2.calculate_rdc_amide_large(path_traj,
                                                   path_top,
                                                   RDC_inp_file,
                                                   minimize_rmsd=False, # Always use the trajectory where RMSD is already minimized
                                                   mode='full'
                                                   )

    reference = np.loadtxt("reference_calculated_full.txt")
    difference  = np.subtract(D_av, reference)
    print("Maximum difference between calculated and average value:  %f" %np.max(difference))
    print("Mean difference between calculated and average value:     %f" %np.mean(difference))
    return()

def test_rdc_average_amide():
    print(" ")
    print("Testing function calculate_rdc_amide_large with mode='average'")
    path_traj = 'trajectory.xtc'
    path_top  = 'topology.pdb'
    RDC_inp_file = "experimental_data.txt" # Experimental data - allignment A28               
    exp_rdc,D_av = nmr2.calculate_rdc_amide_large(path_traj,
                                                   path_top,
                                                   RDC_inp_file,
                                                   minimize_rmsd=False, # Always use the trajectory where RMSD is already minimized
                                                   mode='average'
                                                   )

    reference = np.loadtxt("reference_calculated_average.txt")
    difference  = np.subtract(D_av, reference)
    print("Maximum difference between calculated and average value:  %f" %np.max(difference))
    print("Mean difference between calculated and average value:     %f" %np.mean(difference))
    return()



if __name__=='__main__':
    test_rdc_average()
    test_rdc_full()
    test_rdc_full_amide()
    test_rdc_average_amide()

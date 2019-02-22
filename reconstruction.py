import os
import time
import mdtraj as md

def reconstruct_pulchra(traj,
                        save,
                        pulchra,
                        mode=" ",
                        add_H=False,
                        clear=False,
                        save_traj=False,
                        out_name="PULCHRA_reconst.xtc",
                        verbose=True):
    """
    Function creates a pulchra-reconstructed trajectory based on all-atom trajectory
    
    The function works only if the whole trajectory can be completely loaded in memory.
    
    This code used the launches PULCHRA executable. Reference:
    Rotkiewicz, P., & Skolnick, J. (2008). Fast procedure for reconstruction
    of full-atom protein models from reduced representations. 
    Journal of computational chemistry, 29(9), 1460-5.
    
    It looks like Pulchra flag -h does not work, so hydrogens are added via gromacs 
    pdb2gmx program. 
    
    Args: 
        
        traj      - md_traj trajectory
        save      - path to save ouptut. Should end with "/"
        pulchra   - path to a pulchra executable (eg. ../Software/pulchra) 
        mode      - (default " ") pulchra reconstruction mode,string of flags
                    according to pulchra docs, without a leading hyphen(eg "cs").
        add_H     - (default False) if true, adds hydrogens using pdb2gmx tool. 
                    Works only if mode does not include -s flag (skip sidechain)
        clear     - (default False) if true, remove all the files from the scratch
                    folder after reconstruction

        save_traj - (default False) if true, the reconstructed trajectory is 
                    saved in pdb format
        out_name  - (default "PULCHRA_reconst.pdb") If save_traj=True, the output trajectory
                     will have be named <out_name>. File extension must be included in 
                    the output_name.
        verbose   -  (default True) Enable more output
        
    Return:
    
       rec_traj - reconstructed trajectory

    
    """
    
    
    
    
    if verbose:
        print("NOTE: The function uses Scratch directory:")
        print("%sScratch" %save)
        print("Make sure, that the same direcory is not used simultaneously")
        print("by different jobs")
        print("Input trajectory      :", traj)
        print("Pulchra executable    :", pulchra)
        print("Reconstruction mode   :", mode)
        print("Add hydrogens(pdb2gmx):", add_H)
        print("Delate scratch files  :", clear)
        print("Save pdb              :", save_traj)
        print("  ")
        
    if add_H:  #Need to temporarly disable gmx backup
        gmx_backup = os.system("echo $GMX_MAXBACKUP")
        print(gmx_backup)
        os.environ["GMX_MAXBACKUP"] = "-1"
        print("WARNING! GMX_MAXBACKUP was set to -1. GMX BACKUP disabled")
        print("starting value of gmx_backup was %s" %gmx_backup)
    
    atoms_to_keep = [a.index for a in traj.topology.atoms if a.name == 'CA']
    traj_alpha = traj.atom_slice(atoms_to_keep)
    
    pulchra_start = time.time()
    
    os.system('rm -r %sScratch' %save)
    os.system('mkdir %sScratch'% save)
    rec_traj_framename=[]
    
    for i in range(0,traj_alpha.n_frames): 
        savename=("%sScratch/input%i.pdb"%(save,i))
        traj_alpha[i].save_pdb(savename)
        os.system("./%s -%s %sScratch/input%i.pdb" %(pulchra,mode,save,i))
        if not add_H:
            rec_traj_framename.append('%sScratch/input%i.rebuilt.pdb' %(save,i))
        else:
            
            gmx_backup = os.system("echo $GMX_MAXBACKUP")   #Need to disable gmx backup temporaly
            os.system("echo 8 |pdb2gmx -water none -f %sScratch/input%i.rebuilt.pdb -o %sScratch/output%i.pdb" %(save,i,save,i))
            rec_traj_framename.append('%sScratch/output%i.pdb' %(save,i))
            
        
    # Creating a reconstructed trajectory
    traj_reconst = md.load(rec_traj_framename)

    if add_H:
        os.environ["GMX_MAXBACKUP"] = str(gmx_backup)
        print("GMX_MAXBACKUP set back to %s" %str(gmx_backup))
            
    # Removing scratch files:
    if clear:
        os.system('rm  %sScratch/*' %save)
    
    # Saving a reconstructed trajectory and writing topology
    if save_traj:
        traj_reconst.save("%s%s" %(save,out_name))
        if add_H:
            os.system('cp %sScratch/output0.pdb %stop_reconst.pdb'  %(save,save))
        else:
            os.system('cp %sScratch/input0.rebuilt.pdb  %stop_reconst.pdb'  %(save,save))
    
        
    print(traj_reconst)   

    print("Time to perform Pulchra reconstruction: ", time.time() - pulchra_start, " seconds")
    return(traj_reconst)

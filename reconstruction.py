import os
import time
import mdtraj as md

def reconstruct_pulchra(traj,
                        save,
                        pulchra,
                        mode=" ",
                        add_H=False,
                        minimize=False,
                        minimization_params=None,
                        clear=False,
                        save_traj=False,
                        out_name="PULCHRA_reconst.xtc",
                        verbose=True,
                        rank=0):
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
        minimize  - (default False) if true, energy minimization is performed, with parameters,
                    described in the file `minimization_params'
        minimizaiton_params - (default None)  name of the file with minimization parameters, used if minimize = True
        clear     - (default False) if true, remove all the files from the scratch
                    folder after reconstruction

        save_traj - (default False) if true, the reconstructed trajectory is
                    saved in pdb format
        out_name  - (default "PULCHRA_reconst.pdb") If save_traj=True, the output trajectory
                     will have be named <out_name>. File extension must be included in
                    the output_name.
        verbose   -  (default True) Enable more output
        rank      -  (default 0) rank of process

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
        print("Minimize              :", minimize)
        print("Minimization parameter:", minimization_params)
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

    os.system('rm -r %sScratch%i' %(save,rank))
    os.system('mkdir %sScratch%i'% (save,rank))
    rec_traj_framename=[]

    for i in range(0,traj_alpha.n_frames):
        savename=("%sScratch%i/input%i.pdb"%(save,rank,i))
        traj_alpha[i].save_pdb(savename)
        if verbose:
            os.system("./%s -%s %sScratch%i/input%i.pdb >> pulchra.log" %(pulchra,mode,save,rank,i))
        else:
            os.system("./%s -%s %sScratch%i/input%i.pdb >> /dev/null" %(pulchra,mode,save,rank,i))
        if not add_H:
            rec_traj_framename.append('%sScratch%i/input%i.rebuilt.pdb' %(save,rank,i))
        else:
            gmx_backup = os.system("echo $GMX_MAXBACKUP")   #Need to disable gmx backup temporaly
            if verbose:
                os.system("printf '8\n1\n' |pdb2gmx -water none -his -f %sScratch%i/input%i.rebuilt.pdb -o %sScratch%i/output%i.pdb" %(save,rank,i,save,rank,i))
            else:
                os.system("printf '8\n1\n' |pdb2gmx -water none -his -f %sScratch%i/input%i.rebuilt.pdb -o %sScratch%i/output%i.pdb > /dev/null " %(save,rank,i,save,rank,i))
            if minimize:
                os.system("printf '8\n1\n1\n' |pdb2gmx -f  %sScratch%i/output%i.pdb -ignh -ter -o %sScratch%i/all.gro -i %sScratch%i/all.itp -water none -p %sScratch%i/all.top" %( save,rank,i,save,rank,save,rank,save,rank ))
                os.system("grompp  %s  -c %sScratch%i/all.gro -p %sScratch%i/all.top -o %sScratch%i/all.tpr" %(minimization_params, save,rank,save,rank,save,rank))
                root = os.getcwd()
                os.system("cd %sScratch%i")
                result_code = os.system("mdrun -nov -compact -deffnm all -c all_post_%i.pdb". %(save  rank, i))
                assert result_code == 0 "Incorect exit code: %i " % result_code
                os.chdir(root)
                rec_traj_framename.append("%sScratch%i/all_post_%i.pdb" %(save, rank, i))
            else:
                rec_traj_framename.append('%sScratch%i/output%i.pdb' %(save, rank, i))


    # Creating a reconstructed trajectory
    traj_reconst = md.load(rec_traj_framename)

    if add_H:
        os.environ["GMX_MAXBACKUP"] = str(gmx_backup)
        print("GMX_MAXBACKUP set back to %s" %str(gmx_backup))


    # Saving a reconstructed trajectory and writing topology
    if save_traj:
        traj_reconst.save("%s%s" %(save,out_name))
        if add_H:
            os.system('cp %sScratch%i/output0.pdb %stop_reconst.pdb'  %(save,rank,save))
        else:
            os.system('cp %sScratch%i/input0.rebuilt.pdb  %stop_reconst.pdb'  %(save,rank,save))

    # Removing scratch files:
    if clear:
        os.system('rm  %sScratch%i/*' %(save,rank))

    print(traj_reconst)

    print("Time to perform Pulchra reconstruction: ", time.time() - pulchra_start, " seconds")
    return(traj_reconst)

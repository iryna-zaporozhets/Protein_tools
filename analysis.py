import mdtraj as md
import numpy  as np

def RMSD1to1(traj1,traj2):
    """
    The function calculate rmsd for each pair of frames
    
    traj1 and traj2 must have the same number of frames.
    traj1 is a trajectory for which rmsd is calculated.
    traj2 is a reference trajectory.
    The function returns 1D numpy array, which contains
    one real rmsd value for each frame
    """
    
    traj1_nframes = traj1.n_frames
    traj2_nframes = traj2.n_frames
    print ("trajectory 1 contains %i  frames."  %traj1_nframes)
    print ("trajectory 2 contains %i  frames." %traj2_nframes)
    if traj1_nframes != traj2_nframes:
        raise IOError("Input trajectories contain different number of frames")
    
    rmsd = np.empty(traj1_nframes)
    for i in range(0,traj1_nframes):
        rmsd[i]=md.rmsd(traj2[i],traj1[i])
    return(rmsd)
        
        
    

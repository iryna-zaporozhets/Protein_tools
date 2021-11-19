import mdtraj as md
import numpy  as np
from itertools import combinations

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



def index2res(structure):
    """
    The function creates list of residues based on files with indexes, used to calculate
    J3 coupling constants. This stuff should have been done during J3 calculations, because
    outputing atom indexes is meaningless, as they depend on topology.For example, amide
    nitrogen of the second residue may have different atomic number in pdb, containing only
    backbone,than in pdb, containing both backbone and sidechain.

    In the input file, each line correspond to 4 indexes of atoms, which form a dihedral
    angle, used to compute J3. For Halpha-HN, this list includs atoms, forming phi
    dihedral angle. These are C(i-1),N(i),Ca(i),C(i). We are interested in residue,
    Nitrogen from which participates in forming a bond. So we are looling for the 2 atom
    (1 by Python numeration)

    Input:

           structure - mdtraj frame

    Output:

           residue_index - list of Python-numerated residue ids.
           residue_name

    """

    index, J3_inp = md.compute_J3_HN_HA(structure)
    assert(index.shape[0]==structure.n_residues-1)
    residue_index = []
    residue_name = []

    for i in range(0,index.shape[0]):
        atom_index = int(index[i,1])
        assert(structure.top.atom(atom_index).name=='N')
        residue_index.append(structure.top.atom(atom_index).residue.index)
        residue_name.append(structure.top.atom(atom_index).residue.name)

    return(residue_index,residue_name)



def rms(a,normalization='n'):
    """
    The function returns rms for a given 1d dataset

    Arg:
        a            : 1D numpy array
        normalization: (Default 'n'). Can be either 'n' or 'n-1'.
                       define denominator of rms


    """
    if normalization == 'n':
        norm = a.shape[0]
    elif normalization == 'n-1':
        norm = a.shape[0]-1
    else:
        print("Normalization is not correct, set to default falue")
        norm=a.shape[0]


    calculated_rms = np.sqrt(np.sum(np.square(a))/norm)
    return(calculated_rms)


def q_factor(D_measured,D_calculated,normalization='n'):
    """
    The function calculate Q-factor, based on np.arrays of experimental and back-
    calculated RDC values

    Args:

        D_measured  : 1D numpy array with experimental RDC
        D_calculated: 1D numpy array with calculated RDC
        normalization: (Default 'n'). Can be either 'n' or 'n-1'.
                       define denominator of rms



    Returns:
        q : float, quality factor.

    Notes: function gives different results from rdc_svd.py. rdc_svd.py uses
    another denominator
    """

    q = rms(np.subtract(D_measured,D_calculated),normalization)/rms(D_measured,normalization)
    return q


def best_hummer_q(traj, native=None, beta_const=50, lambda_const=1.2, native_cutoff=1.2):
    """Compute the fraction of native contacts according the definition from
    Best, Hummer and Eaton [1]
    Cde modified from (https://mdtraj.org/1.9.3/examples/native-contact.html)
    Parameters for different coarse-grained models taken from 
    Habibi M, Rottler J, Plotkin SS. As Simple As Possible, but Not Simpler:
    Exploring the Fidelity of Coarse-Grained Protein Models for 
    Simulated Force Spectroscopy. PLoS Comput Biol. 2016;12(11):e1005211. 
    Published 2016 Nov 29. doi:10.1371/journal.pcbi.1005211

    Parameters
    ----------
    traj : md.Trajectory
        The trajectory to do the computation for
    native : md.Trajectory
        The 'native state'. This can be an entire trajecory, or just a single frame.
        Only the first conformation is used
    beta_const : float, 1/nm
          Beta constant. 50 1/nm for  AA/HA-Go, AWSEM, Calpha-Go

    lambda_constant : float
                   1.8 for AA and HA-Go, 1.2 for AWSEM and Calpha-Go

    native_cutoff : float, nm
                    0.48 for AA and HA-Go, 0.6 for AWSEM, 1.2 for Ca-Go
        
    Returns
    -------
    q : np.array, shape=(len(traj),)
        The fraction of native contacts in each frame of `traj`
        
    References
    ----------
    ..[1] Best, Hummer, and Eaton, "Native contacts determine protein folding
          mechanisms in atomistic simulations" PNAS (2013)
    """
  
    # get the indices of all of the heavy atoms
    heavy = native.topology.select_atom_indices('heavy')
    # get the pairs of heavy atoms which are farther than 3
    # residues apart
    heavy_pairs = np.array(
        [(i,j) for (i,j) in combinations(heavy, 2)
            if abs(native.topology.atom(i).residue.index - \
                   native.topology.atom(j).residue.index) > 3])
    
    # compute the distances between these pairs in the native state
    heavy_pairs_distances = md.compute_distances(native[0], heavy_pairs)[0]
    # and get the pairs s.t. the distance is less than NATIVE_CUTOFF
    native_contacts = heavy_pairs[heavy_pairs_distances < native_cutoff]
    print("Number of native contacts", len(native_contacts)) 
    # now compute these distances for the whole trajectory
    r = md.compute_distances(traj, native_contacts)
    # and recompute them for just the native state
    r0 = md.compute_distances(native[0], native_contacts)
    
    q = np.mean(1.0 / (1 + np.exp(beta_const * (r - lambda_const * r0))), axis=1)
    return np.expand_dims(q, axis=1)  

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

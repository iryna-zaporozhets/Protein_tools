import numpy as np
import os
import re
import mdtraj as md

########################################################################################
#         
#   The module contains classes and functions needed to calculate NMR data from
#   md traj trajectories.
#   Now implemented calculate_rdc function 
#
#
#
############# RDC calculating function and all the dependencies #########################

class Bond:
    
    def __init__(self, resid_i, resname_i, atomname_i,
                       resid_j, resname_j, atomname_j):
        self.resid_i = resid_i
        self.resname_i = resname_i
        self.atomname_i = atomname_i
        self.resid_j = resid_j
        self.resname_j = resname_j
        self.atomname_j = atomname_j
        
    def display(self):
        print("resid_i = ", self.resid_i)
        print("resname_i = ", self.resname_i)
        print("atomname_i = ", self.atomname_i)
        print("resid_j = ", self.resid_j)
        print("resname_j = ", self.resname_j)
        print("atomname_i = ", self.atomname_j)


############################################################################################

def vector(atom_pair, frame):
    """
    Calculate normalized vectors for a specific pair if atoms in a given frame
    """
    
    vec = frame.xyz[0,atom_pair[1],:]-frame.xyz[0,atom_pair[0],:]
    vec = np.divide(vec, np.sqrt(vec[0]**2 
                               + vec[1]**2 
                               + vec[2]**2)
                      )
    return vec

##############################################################################################

def bilin(vector):
    """
    Calculate vector of bilinear components based on coordinate vector
    
    """
    bil = np.empty(5)
    bil[0] =   vector[0]**2-vector[2]**2
    bil[1] =   vector[1]**2-vector[2]**2
    bil[2] = 2*vector[0]*vector[1]
    bil[3] = 2*vector[0]*vector[2]
    bil[4] = 2*vector[1]*vector[2]
    
    return bil

def test_bilin():
    test_vector = np.array([1,2,4])
    result = np.array([-15, -12,   4,   8,  16])
    assert(np.array_equal(bilin(test_vector),result))
    print ("Passed successfully")

##################################################################################################

def bilin_matrix(bond_selections, test_frame):

    """
    Creates a matrix containing bilinear terms for all bonds,
    listed in bond_selection list"
    
    Input:   1) bond_selections
                    List, which contain 1 element,
                    for each bond, for which experimental RDC is known.
                    Each element is a list of lenght 2.
                    The first element of a list - integer,
                    corresponding to index of the first atom in the bond.
                    The second element of a list - integer,
                    corresponding to index of the second atom in the bond.
                    Numeration - according to the numeration in MDTraj topology
                    
             2)  test_frame
                    A single MDtraj frame
                    
                    
    Output:   F - numpy array with shape ((len(bond_selections),5))
              Includes bilinear components x^2-z^2, y^2-z^2, 2xy, 2xz, 2yz
              
             
    """
    F = np.empty((len(bond_selections),5))
    i=0
    for (pair) in bond_selections:
        vec = vector(pair,test_frame)
        bil = bilin(vec)
        F[i,:]=bil
        i+=1
    return F
#####################################################################################################

#####################################################################################################
def calculate_rdc_large(traj,topology,RDC_inp_file, minimize_rmsd=True,mode='average'):
    
    """
    Calculate residual dipolar couplings based on SVD for a long trajectory,
    when the trajectory cannot be loaded in the memory as a whole. 
    
    Args:
    
       traj      : trajectory file in any format, supported by md_traj
       
       topology  : topology file (the same as one for mdtraj)
        
       RDC_inp_file  : file with experimental RDC values. 
                      File format is close to NMRPipe, but NOT EXACTLY the same.
                      https://www.ibbr.umd.edu/nmrpipe/install.html
                    
                      Coulumn descriptions:
                      1:  Residue i id 
                      2:  Residue i name
                      3:  Atom    i name
                      4:  Residue j id
                      5:  Residue j name
                      6:  Atom    j name
                      7:  RDC for the bond i-j
                    
        minimize_RMSD: Determing, whether superimposion with respect to the frame,
                       minimizing RMSD, will be done.
                       
                       if minimize_RMSD=True (default)  a frame, minimizing
                       sum of  C_alpha RMSD with respect to all other frames is found
                       
                       if minimize_RMSD=False  0-th frame is used to as a reference to super
                       impose all other frames
        

        
    Return: two numpy arrays: 
              exp_rdc - experimental values of RDCs
              D_av    -  back-calculated  RDCs
    
    
    Dependencies: 
                Packages  :   re, np, md
                Classes   :   Bond
                Functions :   bilin_matrix, vector, bilin
    """
    structure=md.load(topology)
    RDC_input = open(RDC_inp_file,'r')
    
    RDCs=[]
    bonds=[]
    for line in RDC_input.readlines():
        match = re.search(r'^\s*(?P<resid_i>[0-9]+)'
                            '\s*[A-Z]{3}'
                            '\s*(?P<name_i>[A-Z\#]{1,4})'
                            '\s*(?P<resid_j>[0-9]+)'
                            '\s*[A-Z]{3}'
                            '\s*(?P<name_j>[A-Z\#]{1,4})'
                            '\s*(?P<rdc>-?[0-9\.]+)', line)
        if match is None:
            continue
        RDCs.append(float(line.split()[6]))  # Work only when RDC are in the 7 coulumn!
        bonds.append(Bond( int(line.split()[0]),
                              (line.split()[1]),
                              (line.split()[2]),
                           int(line.split()[3]),
                              (line.split()[4]),
                              (line.split()[5]))
                    )
    
    RDC_input.close()

  
    bond_selections=[]
    for bond in bonds:
        # -1 correspond to transition between PDB numeration and MDTRAJ numeration
        # Names of atoms in input files should be the same as one used by mdtraj
        selection_i = structure.top.select('resid %i and name %s' %(bond.resid_i-1, bond.atomname_i))
        assert(selection_i.size != 0)
        selection_j = structure.top.select('resid %i and name %s' %(bond.resid_j-1, bond.atomname_j))
        assert(selection_j.size != 0)
        bond_selections.append([selection_i[0], selection_j[0]])
    
    # According to the procedure, described in Olsson2017 papper, need to find a frame, 
    # which minimizes sum of  C_alpha RMSD with respect to all other frames
    # Quadratic algorithm  (O(N2))???
    # Should use C-alpha trajectory!
    #  RMSD matrix is symmetric 
    
    ### From this part, need to use iterator. How can we find superimposed configuration?
    ### Solution: Use a trajectory, that has already been superimposed
    print("NOTE: input trajectory should be superimposed")
    if minimize_rmsd: 
        print("WARNING! RMSD minimization is not implemented in current function yet")
        print("Use superimposed trajectory as an input")

    F_av = np.zeros((len(bond_selections),5))
    n_of_frames=0
    print(topology)
    for chunks in md.iterload(traj, chunk=1,top=topology):
        n_of_frames+=1
        F_av = F_av +  bilin_matrix(bond_selections, chunks)
    F_av = np.divide(F_av,n_of_frames)
    A_av, residuals,  rank,s = np.linalg.lstsq(F_av,np.array(RDCs),rcond=-1)


    if mode=='average':
        D_av=np.dot(F_av,A_av)
        
    if mode=='full':
        D_full=[]
        for chunks in md.iterload(traj, chunk=1,top=topology):
            F=bilin_matrix(bond_selections,chunks)
            D=np.dot(F,A_av)
            D_full.append(D)
        D_av=np.array(D_full)
    
    exp_rdc=np.array(RDCs)
    return(exp_rdc,D_av)
####################################################################################################

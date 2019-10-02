"""
The module contains functions for parsing SMOG contacts and converting them
to coarse-grained representation
"""

# Functions  check_equal, check_pair_in_list, add_pair,create_CACB_exclusions,
# make_pairwise_files were originally done by Justing Chen, with further
# modifications

def check_equal(p1, p2):
    if p1[0] == p2[0] and p1[1] == p2[1]:
        return True
    else:
        return False

def check_pair_in_list(pair_list,check):
    matching = False
    for last_pair in pair_list:
        if check_equal(last_pair, check ):
            matching = True

    return matching

def add_pair(pairs, atm1, atm2,pairs_dict={}):
    new_pair = [atm1.index+1, atm2.index+1]
    if len(pairs) == 0:
        pairs.append(new_pair)
        pairs_dict[tuple(new_pair)] = 1
    else:
        matching = False
        for last_pair in pairs:
            if check_equal(new_pair, last_pair):
                matching = True
                pairs_dict[tuple(new_pair)] += 1
        if not matching:
            pairs.append(new_pair)
            pairs_dict[tuple(new_pair)] = 1

def create_CACB_exclusions(all_atom_pdb, cacb_atom_pdb, contacts, cutAA=4, cutBB=2, cutAB=2):
    """Parse an all-atom contact file fro SMOG and covert to CaCb contacts
    Only keep Ca-Ca pairs and Cb-Cb pairs, exclude all others

    Rules for exclusions are based off the following reference:
    Cheung,et.al. 'Exploring the interplay between topology and secondary structural
        formation in the protein folding problem'. J.Chem.Phys. B. 2003.

    Parameters
    ----------
    all_atom_pdb : str
        String corresponding to the all-atom pdb file for loading in mdtraj
    cacb_atom_pdb : str
        String corresponding to cacb-pdb pdb file for loading in mdtraj
    contacts : array Nx2
        Array containing atom contacts, index base 1 from SMOG.

    """

    # Exclude neighbors closer in sequence than:
    #   |i - j| < 4 for CA_i CA_i pairs
    #   |i - j| < 2 for CB_i CB_j pairs
    #   |i - j| < 2 for CA_i CB_j pairs
    cutAA = cutAA
    cutAB = cutAB
    cutBB = cutBB
    traj = md.load(all_atom_pdb)
    top = traj.top
    cacbtraj = md.load(cacb_atom_pdb)
    cacb_top = cacbtraj.top
    contacts_zero = contacts - 1
    assert np.shape(contacts_zero)[1] == 2
    num_pairs = np.shape(contacts_zero)[0]
    # Loop over all pairs, then add to the new_pairs list
    new_pairs = []
    backbone_idx = top.select("backbone")
    print(backbone_idx)
    sidechain_idx = top.select("sidechain")
    last_pair = [-10, -10]

    pairs_dictionary = {} # Dictionary, wich matches CG pairs and number of
                         # all-atom contacts
    for pair_idx in range(num_pairs):
        idx1 = contacts_zero[pair_idx, 0]
        idx2 = contacts_zero[pair_idx, 1]

        atm1 = top.atom(idx1)
        atm2 = top.atom(idx2)
        resid1 = atm1.residue.index
        resid2 = atm2.residue.index

        if (atm1.index in backbone_idx) and (atm2.index in backbone_idx):
            # both backbone atoms, add calpha-interactions
            print("Adding Backbone Interaction")
            if (resid2 - resid1) >= cutAA:
                new_atm1 = cacb_top.residue(resid1).atom(0)
                new_atm2 = cacb_top.residue(resid2).atom(0)
                add_pair(new_pairs, new_atm1, new_atm2, pairs_dict=pairs_dictionary)
        elif (atm1.index in sidechain_idx) and (atm2.index in sidechain_idx):
            # both side chain atoms, add a cbeta-interaction
            print("Adding Sidechain Interaction")
            if (resid2 - resid1) >= cutBB:
                new_atm1 = cacb_top.residue(resid1).atom(1)
                new_atm2 = cacb_top.residue(resid2).atom(1)
                add_pair(new_pairs, new_atm1, new_atm2,pairs_dict=pairs_dictionary)
        else:
            # a sidechain-backbone interaction. Print warning and move on
            print("Not Adding sidechain-backbone interaction. Ignoring")

    new_pairs = np.array(new_pairs)
    print(new_pairs)
    old_size =  np.shape(new_pairs)[0]
    # sort the array
    sorted_indices = np.argsort(new_pairs[:,0])
    new_pairs_list = []
    next_sort_pairs = []
    next_sort_indices = []
    current_idx = new_pairs[sorted_indices[0], 0]
    temporary = []
    make_next_idx = False
    for sort_idx in sorted_indices:
        if new_pairs[sort_idx,0] == current_idx:
            temporary.append(new_pairs[sort_idx,:])
        else:
            new_array = np.array(temporary)
            try:
                assert np.shape(new_array)[1] == 2
            except:
                print(np.shape(new_array))
            next_sort_pairs.append(np.copy(new_array))
            new_sorted = np.argsort(new_array[:,1])
            next_sort_indices.append(new_sorted)
            temporary = [new_pairs[sort_idx,:]]
            current_idx = new_pairs[sort_idx,0]
    new_array = np.array(temporary)
    next_sort_pairs.append(np.copy(new_array))
    new_sorted = np.argsort(new_array[:,1])
    next_sort_indices.append(new_sorted)

    for idx,pair_sets in enumerate(next_sort_pairs):
        assert np.shape(pair_sets)[0] == np.shape(next_sort_indices[idx])[0]
        for jdx in next_sort_indices[idx]:
            new_pairs_list.append(pair_sets[jdx,:])

    new_pairs = np.array(new_pairs_list)
    try:
        assert np.shape(new_pairs)[0] == old_size
    except:
        print(np.shape(new_pairs)[0])
        print(old_size)
        #raise
    number_of_aa_contacts = []
    for pair in new_pairs:
        number_of_aa_contacts.append(pairs_dictionary[tuple(pair)])
    return new_pairs, pairs_dictionary


def make_pairwise_files(cacb_pdb,
                        pairs,
                        pairwise_params='pairwise_params',
                        model_params='model_params',type='native'):
    # Exclude neighbors closer in sequence than:
    #   |i - j| < 4 for CA_i CA_i pairs
    #   |i - j| < 2 for CB_i CB_j pairs
    #   |i - j| < 2 for CA_i CB_j pairs
    cutAA = 4
    cutAB = 2
    cutBB = 2

    traj = md.load(cacb_pdb)
    top = traj.top

    pairs_zero_idx = pairs - 1
    fpp = open(pairwise_params, "w")
    fmp = open(model_params, "w")
    fpp.write("#    pairs         param         potential_type      other_params\n")
    fmp.write("# model params\n")
    backbone =  top.select("backbone")
    sidechain =  top.select("sidechain")
    count = 0
    for i in range(top.n_atoms):
        for j in range(i, top.n_atoms):
            atm1 = top.atom(i)
            atm2 = top.atom(j)
            res_diff = atm2.residue.index - atm1.residue.index
            compute = False
            inter_type = None
            if (atm1.index in backbone) and (atm2.index in backbone):
                inter_type = "backbone"
                if res_diff >= cutAA:
                    compute = True
            elif (atm1.index in sidechain) and (atm2.index in sidechain):
                inter_type = "sidechain"
                if res_diff >= cutBB:
                    compute = True

            if compute:
                if inter_type == "backbone":
                    excluded_volume = 0.266 #(1.9*1.4)
                elif inter_type == "sidechain":
                    v1 = atom_types.residue_cacb_effective_interaction[atm1.residue.name]
                    v2 = atom_types.residue_cacb_effective_interaction[atm2.residue.name]
                    excluded_volume = (v1 + v2) / 2.

                if check_pair_in_list(pairs_zero_idx,[i,j]):
                    # then a native pair
                    dist = md.compute_distances(traj, [[i,j]])[0,0]
                    native_param = 1.0
                    potential = "    LJ12GAUSSIAN"
                else:
                    continue
                    dist = excluded_volume + 0.2
                    native_param = 0.0
                    potential = "LJ12GAUSSIANTANH"
                if dist < (excluded_volume+0.2):
                    print("Warning: Distance for %d %d is less than the excluded volume, %f versus %f" % (i+1, j+1, dist, excluded_volume))
                # now write out the pairwise params files
                fpp.write("%6d  %6d  %12d       %s     %.6f  %.6f  %.6f\n" % (i+1, j+1, count, potential, excluded_volume, dist, 0.05))
                fmp.write("%.6f\n" % native_param)
                count += 1
    fpp.close()
    fmp.close()


def write_pairs_dictionary(pairs_dictionary,filename):
    """
    The function writes pairs_dictionary to the file filename
    with the format resid1 resid2 number_of_contacts
    """
    out_file = open(filename,'wt')
    for pair, num_of_count in pairs_dictionary.items():
        out_file.write("%5d  %5d   %5d \n" %(pair[0], pair[1],num_of_count))
    out_file.close()
    return 0





def read_pairs_dictionary(filename):
    """
    The function reads pairs_dictionary, that discribes number of all-atom contacts, that correspond to
    a particular  mutant. The result is dictionary of format [(atom i, atom j) : number_of_all_atom_contacts]

    parameters
    ---------

    filename : str
     path to the file, that contains dictionary in the format atom i, atom j, number_of_all_atom_contacts


    return
    ------
    pairs_dictionary : dict
    """

    pairs_dictionary = {}
    inp_file = open(filename,'rt')
    for lines in inp_file:
        data = lines.split()
        pair = (int(data[0]),int(data[1]))
        pairs_dictionary[pair] = int(data[2])
    return pairs_dictionary

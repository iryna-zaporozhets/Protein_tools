import pdb_mutator
import SMOG_contact_parser
import md_nmr2 as nmr
import numpy as np


def test_find_atoms_to_delete():
    """
    Set of tests for find_atoms_to_delete function
    """

    assert pdb_mutator.find_atoms_to_delete('ILE','VAL') == ['CD1']
    assert pdb_mutator.find_atoms_to_delete('ALA','GLY') == ['CB']
    assert pdb_mutator.find_atoms_to_delete('LYS','ALA') == ['CG','CD','CE','NZ']
    assert pdb_mutator.find_atoms_to_delete('VAL','ALA') == ['CG1','CG2']
    assert pdb_mutator.find_atoms_to_delete('THR','ALA') == ['OG1','CG2']
    assert pdb_mutator.find_atoms_to_delete('LEU','ALA') == ['CG','CD1','CD2']
    assert pdb_mutator.find_atoms_to_delete('GLN','ALA') == ['CG','CD','OE1','NE2']

def test_decode_mutation():
    """
    decode_mutation
    """
    assert pdb_mutator.decode_mutation('V11498I') == (11498,'VAL','ILE')
    assert pdb_mutator.decode_mutation('Q1R') == (1,'GLN','ARG')

def test_find_mutation_contacts():
    """
    find_mutation_contacts
    """
    "Test uses a single topology file that involves one aminoacid"
    topfile = 'test_find_mutation_contacts/topology_LYS_THR.pdb'
    contacts = np.array([[1,2],
                        [1,4],
                        [1,7],
                        [1,9],
                        [1,13],
                        [1,14],
                        [2,10],
                       [2,14],
                       [2,16],
                       [3,8],
                       [3,15],
                       [3,12],
                       [4,5],
                       [4,7],
                       [4,9],
                       [4,11],
                       [4,15],
                       [5,7],
                       [5,9],
                       [5,15],
                       [8,16],
                       [6,8],
                       [6,11],
                       [6,14],
                       [6,15],
                       [7,9],
                       [7,12],
                       [7,13],
                       [7,14],
                       [8,9],
                       [8,10],
                       [8,14],
                       [8,16],
                       [9,10],
                       [9,13],
                       [9,15],
                       [10,11],
                       [10,13],
                       [10,16],
                       [11,13],
                       [11,15]])

    K1G_contacts = np.array ([[1,2],
                            [1,4],
                            [1,13],
                            [1,14],
                            [2,10],
                           [2,14],
                           [2,16],
                           [3,15],
                           [3,12],
                           [4,11],
                           [4,15],
                           [10,11],
                           [10,13],
                           [10,16],
                           [11,13],
                           [11,15]
                           ])

    K1A_contacts = np.array([[1,2],
                            [1,4],
                            [1,13],
                            [1,14],
                            [2,10],
                           [2,14],
                           [2,16],
                           [3,15],
                           [3,12],
                           [4,5],
                           [4,11],
                           [4,15],
                           [5,15],
                           [10,11],
                           [10,13],
                           [10,16],
                           [11,13],
                           [11,15]])

    T2A_contacts = np.array([[1,2],
                            [1,4],
                            [1,7],
                            [1,9],
                            [1,13],
                            [1,14],
                            [2,10],
                           [2,14],
                           [3,8],
                           [3,12],
                           [4,5],
                           [4,7],
                           [4,9],
                           [4,11],
                           [5,7],
                           [5,9],
                           [6,8],
                           [6,11],
                           [6,14],
                           [7,9],
                           [7,12],
                           [7,13],
                           [7,14],
                           [8,9],
                           [8,10],
                           [8,14],
                           [9,10],
                           [9,13],
                           [10,11],
                           [10,13],
                           [11,13]
                           ])
    T2G_contacts = np.array([[1,2],
                                [1,4],
                                [1,7],
                                [1,9],
                                [1,13],
                                [2,10],
                               [3,8],
                               [3,12],
                               [4,5],
                               [4,7],
                               [4,9],
                               [4,11],
                               [5,7],
                               [5,9],
                               [6,8],
                               [6,11],
                               [7,9],
                               [7,12],
                               [7,13],
                               [8,9],
                               [8,10],
                               [9,10],
                               [9,13],
                               [10,11],
                               [10,13],
                               [11,13]
                               ])

    assert  np.all(np.equal(SMOG_contact_parser.find_mutation_contacts('K1G',contacts,topfile)
                           ,K1G_contacts)
                            )
    assert  np.all(np.equal(SMOG_contact_parser.find_mutation_contacts('K1A',contacts,topfile)
                               ,K1A_contacts)
                                )
    assert  np.all(np.equal(SMOG_contact_parser.find_mutation_contacts('T2A',contacts,topfile)
                               ,T2A_contacts)
                                )
    assert  np.all(np.equal(SMOG_contact_parser.find_mutation_contacts('T2G',contacts,topfile)
                               ,T2G_contacts)
                                )

def test_create_CACB_exclusions():
    '''
    SMOG_contact_parser.create_CACB_exclusions
    '''
    contacts = np.array ([
                        [9,19],[12,18],[11,21],
                        [17,22],
                        [10,23],[17,20],
                        [1,18], [4,20], [4,21],
                        [5,22],[5,23],[5,24],
                        [1,25],
                        [1,26], [2,26], [3,29], [4,27],
                        [5,30], [5,36],
                        [3,30], [5,26],
                        [2,37], [1,37], [3,38],
                        [5,43], [8,43], [8,42], [7,42],
                        [3,41]
                       ])

    # [9,19],[12,18],[11,21],         # CaCa 2-3: [3,5] (3 contacts)
    # [17,22],                        # CbCb 2-3: [4,6] (1 contacts)
    # [10,23],[17,20],                # CaCb 2-3: [3,6],[4,5] (2 contacts)
    # [1,18], [4,20], [4,21],         # CaCa 1-3: [1,5]  (3 contacts)
    # [5,22],[5,23],[5,24],           # CbCb 1-3: [2,6]  (3 contacts)
    # [1,25],                         # CaCb 1-3: [1,6]  (1 contact)
    # [1,26], [2,26], [3,29], [4,27]  # CaCa 1-4: [1,7] (4 contacts)
    # [5,30], [5,36],                 # CbCb 1-4: [2,8] (2 contacts)
    # [3,30], [5,26],                 # CaCb 1-4: [1,8],[2,7] (2 contacts)
    # [2,37], [1,37], [3,38],         # CaCa 1-5: [1,9]   (3 contacts)
    # [5,43], [8,43], [8,42], [7,42], # CbCb 1-5:  [2,10] (4 contacts)
     #[3,41],                         # CaCb 1-5: [1,10] ( 1 contact)
     # ]

     # Trial 1. cutAA, cutBB, cutAB =2
    new_pairs, pairs_dictionary = SMOG_contact_parser.create_CACB_exclusions(
                            'test_create_CACB_exclusions/All_atom.pdb',
                            'test_create_CACB_exclusions/CaCb.pdb',
                            contacts,
                            cutAA=2,
                            cutBB=2,
                            cutAB=2)
    print (new_pairs)
    assert np.all(np.equal(new_pairs, np.array([[1,5],
                                                [1,7],
                                                [1,9],
                                                [2,6],
                                                [2,8],
                                                [2,10]]
                                                )))

    target_pair_dictionary = {
    (1,5): 3,
    (1,7): 4,
    (1,9): 3,
    (2,6): 3,
    (2,8): 2,
    (2,10): 4
    }

    assert len(pairs_dictionary)==len(target_pair_dictionary)
    for pair in pairs_dictionary:
        assert pairs_dictionary[pair] == target_pair_dictionary[pair]

    # Trial 2. cutAA = 4 cutBB, cutAB =1
    new_pairs, pairs_dictionary = SMOG_contact_parser.create_CACB_exclusions(
                            'test_create_CACB_exclusions/All_atom.pdb',
                            'test_create_CACB_exclusions/CaCb.pdb',
                            contacts,
                            cutAA=4,
                            cutBB=3,
                            cutAB=1)
    print (new_pairs)
    assert np.all(np.equal(new_pairs, np.array([[1,9],
                                                [2,8],
                                                [2,10]]
                                                )))

    target_pair_dictionary = {
    (1,9): 3,
    (2,8): 2,
    (2,10): 4
    }

    assert len(pairs_dictionary)==len(target_pair_dictionary)
    for pair in pairs_dictionary:
        assert pairs_dictionary[pair] == target_pair_dictionary[pair]

def test_normalize():
    print(dir(nmr))
    cutoff = 1e-14
    vec1 = [3.0,4.0,5.0]
    normalized = nmr.normalize(vec1)
    vector_norm = np.linalg.norm(vec1)
    for i in range(3):
        assert (normalized[i] - vec1[i]/vector_norm ) < cutoff

    vec2 = [0.00,0.00,1.00]
    normalized = nmr.normalize(vec2)
    for i in range(3):
        assert (normalized[i] - vec2[i]) < cutoff

def test_bilin():
    test_vector = np.array([1,2,4])
    result = np.array([-15, -12,   4,   8,  16])
    assert(np.array_equal(nmr.bilin(test_vector),result))
    print ("Passed successfully")

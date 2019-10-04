import pdb_mutator
import SMOG_contact_parser
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

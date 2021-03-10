import mdtraj as md
import re
from .pdb_dicts import dict_of_heavy_atoms, names_3_to_1, names_1_to_3


def find_atoms_to_delete(res1,res2):
    """
    The function compare list of havy atoms of res1 and res2  and finds,
     which atoms are present in res1 but not in res2.
    The function also checkes, that all the atoms in res2 are present in res2.

    parameters
    ----------
    res1 : str
    3-letter name of the first residue, upper-case
    res 2 : str
    3-letter name of the second residue, upper-case


    return
    ------
    atoms_to_delete_names : list
    list of strings, each entry corresponds to the name of atom, which is
    present in the first residue, but is absend in the second one

    dependency
    ----------
    dict_of_heavy_atoms is required to be present in global name space

    """
    for atom in dict_of_heavy_atoms[res2]:
        if atom not in dict_of_heavy_atoms[res1]:
            raise ValueError("Mutation is not allowed. The second residue has atoms, not present in the first residue")


    atoms_to_delete_names = []
    for atom in dict_of_heavy_atoms[res1]:
        if atom not in dict_of_heavy_atoms[res2]:
            atoms_to_delete_names.append(atom)

    return atoms_to_delete_names


def remove_pdb_entry(pdb_file,entry_list,output_file):
    """
    The function generate a new pdb_file with
    entries, listed in entry_list, removed

    parameters
    ----------
    pdb_file : str
    File, that should be corrected

    entry_list : list
    list of pdb entries, that should be removed

    output_file : name of the output file

    """

    pdb = open(pdb_file,'rt')
    output_pdb = open(output_file,'wt')
    for line in pdb:
        output = True
        for entry in entry_list:
            if re.match(entry,line) != None:
                output = False
        if output:
            output_pdb.write(line)
        else:
            print ("The line removed: %s" %line)

    pdb.close()
    output_pdb.close()

    return 0


def mutate_structure(input_pdb,resid_to_mutate,resname_before,resname_after,wt_structure_name,clean=True):
    """
    The script takes pdb file and write a pdb with mutation.
    Only mutations, that can be obtained by removing one or several heavy atoms are acceptable


    parameters
    ----------

    input_pdb : str
    path to a pdb file, that contains a  wild-type structure. The input structure should contain only heavy atoms (i.e., non-H atoms)
    resid_to_mutate : int
    number of residue, that should be mutated. Only one mutation per run. Number of the residue should be given in PDB indexing
    (starting from 1)

    resname_before : str
    3-letter name of the residue, that should be mutated.
    Should be the same, as a residue name in input_pdb structure. Otherwise, AssertionError will apper.

    resname_after : str
    3-letter name of the target residue.
    One should be able to obtain the target residue by deleting heavy atoms from resname_before.
    Otherwise, ValueError will appear

    wt_structure_name : str
    root name to use for output pdb files.

    clean : bool (True)
    If clean, additional file without entries TER and ENDMDL is generated.
    Cleaned files can further be used with SMOG webserver.

    """
    structure = md.load(input_pdb)
    atoms_to_delete_names = find_atoms_to_delete(resname_before,resname_after)


    #================================Convert everyting to python notation=================
    resid_to_mutate_python = resid_to_mutate - 1

    for residue in structure.top.residues:
        if residue.index == resid_to_mutate_python: # Need to take into account transition from  PDB notation to Python notation
            assert residue.name == resname_before, "Name and index of a residue, that should be mutated, are inconsistent"
            residue.name = resname_after


    atoms_to_save_ndx = []
    for atom in structure.top.atoms:
        if not ((atom.residue.index == resid_to_mutate_python ) and(atom.name in atoms_to_delete_names)):
            atoms_to_save_ndx.append(atom.index)

    mutant_structure = structure.atom_slice(atoms_to_save_ndx)
    structure_name = '%s_%s%d%s.pdb' %(wt_structure_name,names_3_to_1[resname_before],resid_to_mutate, names_3_to_1[resname_after])
    mutant_structure.save(structure_name)

    # To work with SMOG, need to remove two entries from the pdb file:
    if clean:
        revome_pdb_entry(structure_name,['TER','ENDMDL'],'%s_%s%d%s_cleaned.pdb' %(wt_structure_name,names_3_to_1[resname_before],resid_to_mutate, names_3_to_1[resname_after]))
    return 0

def decode_mutation(mutation_code):
    """
    The function decodes mutation code in format X00Y, where X- one-letter code
    of initial aminoacid, Y - one-letter code of aminoacid upon mutation, 00 represents
    number of aminoacid residue (in PDB format, indexed from 1)
    returns a tuple (number of aminoacid residue 00, in PDB format,
     X in 3-letter format, Y in 3-letter format)
    """
    return (int(mutation_code[1:len(mutation_code)-1]),
            names_1_to_3[mutation_code[0]],
            names_1_to_3[mutation_code[-1]])

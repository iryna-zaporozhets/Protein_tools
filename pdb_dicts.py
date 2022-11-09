"""
This module has dictionaries and constants, usefull when working with protein
structures and pdb files.
"""

"""
Dictionary dict_of_heavy_atoms maps 3-letter aminoacid code to the list of
all heavy-atom namesas used in pdb files. S
The names are taken from:
    https://cdn.rcsb.org/wwpdb/docs/documentation/file-format/PDB_format_1992.pdf

"""
backbone_atoms = ['N', 'CA', 'C', 'O']
sidechain_atoms = { "ALA" : ['CB'],
                    "ARG" : ['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
                    "ASN" : ['CB', 'CG', 'OD1', 'ND2'],
                    "ASP" : ['CB', 'CG', 'OD1', 'OD2'],
                    "CYS" : ['CB', 'SG'],
                    "GLU" : ['CB', 'CG', 'CD', 'OE1', 'OE2'],
                    "GLN" : ['CB', 'CG', 'CD', 'OE1', 'NE2'],
                    "GLY" : [],
                    "HIS" : ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'],
                    "ILE" : ['CB', 'CG1', 'CG2', 'CD1'],
                    "LEU" : ['CB', 'CG', 'CD1', 'CD2'],
                    "LYS" : ['CB', 'CG', 'CD', 'CE', 'NZ'],
                    "MET" : ['CB', 'CG', 'SD', 'CE'],
                    "PHE" : ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
                    "PRO" : ['CB', 'CG', 'CD'],
                    "SER" : ["CB", 'OG'],
                    "THR" : ['CB', 'OG1', 'CG2'],
                    "TRP" : ['CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
                    "TYR" : ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
                    "VAL" : ['CB', 'CG1', 'CG2']
                    }
dict_of_heavy_atoms = {item : backbone_atoms + value for item, value in sidechain_atoms.items() }


"""
Dictionary names_3_to_1 maps 3-letter residue codes to 1-letter codes.
Reference: https://cdn.rcsb.org/wwpdb/docs/documentation/file-format/PDB_format_1996.pdf
"""
names_3_to_1 = { "ALA" : "A",
                 "ARG" : "R",
                 "ASN" : "N",
                 "ASP" : "D",
                 "ASX" : "B",
                 "CYS" : "C",
                 "GLN" : "Q",
                 "GLU" : "E",
                 "GLX" : "Z",
                 "GLY" : "G",
                 "HIS" : "H",
                 "ILE" : "I",
                 "LEU" : "L",
                 "LYS" : "K",
                 "MET" : "M",
                 "PHE" : "F",
                 "PRO" : "P",
                 "SER" : "S",
                 "THR" : "T",
                 "TRP" : "W",
                 "TYR" : "Y",
                 "VAL" : "V"
                 }

"""
Dictionary names_1_to_3 maps 1-letter residue codes to 3-letter codes.
Generated based on names_3_to_1 dictionary
Reference: https://cdn.rcsb.org/wwpdb/docs/documentation/file-format/PDB_format_1996.pdf
"""
names_1_to_3 = {value:key for (key,value) in  names_3_to_1.items()}

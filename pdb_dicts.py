"""
This module has dictionaries and constants, usefull when working with protein
structures and pdb files.
"""

"""
Dictionary dict_of_heavy_atoms maps 3-letter aminoacid code to the list of
all heavy-atom namesas used in pdb files. Only ALA, VAL, GLN, GLY, ILE,
LEU, LYS, THR are included at this point
"""
dict_of_heavy_atoms = { "ALA" : ['N','CA','C','O','CB'],
                        "VAL" : ['N','CA','C','O','CB','CG1','CG2'],
                        "GLN" : ['N','CA','C','O','CB','CG','CD','OE1','NE2'],
                        "GLY" : ['N','CA','C','O'],
                        "ILE" : ['N','CA','C','O','CB','CG1','CG2','CD1'],
                        "LEU" : ['N','CA','C','O','CB','CG','CD1','CD2'],
                        "LYS" : ['N','CA','C','O','CB','CG','CD','CE','NZ'],
                        "THR" : ['N','CA','C','O','CB','OG1','CG2']
    }


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
names_3_to_1 = {value:key for (key,value) in  names_1_to_3.items()}

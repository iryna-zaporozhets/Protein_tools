import re
import numpy as np


def convert_top(raw_topology,save_folder='.'):
    """ The function takes gromacs topology file and splits it into two 
        other files. The first one includes topology components, that can
        be directly loaded into openMM. 
        The second file includes definition of pairs. This function to be used when
        fuction type 6 is used in pairs.
        This function type is present in "extended" GROMACS, but is not supported
        by OpenMM so it needs to be loaded as a custom type.
        This function also creates a pairwise params file in the same format as one used previously.

       args
       ----
       raw_topology - GROMACS topology file, prepared by SMOG web-server
                      http://smog-server.org/
       save_folder  - Folder to save outputs (without trailing slash).

       returns
       -------

       topology   - name of the topology file, that can be loaded automatically,
                    using GromacsTopFile class of OpenMM
       pairs      - name of the file, which contains pairs to be loaded manually
   
    """

    top_file = open(raw_topology,"r")
    topology = open("%s/topology.top" %save_folder,"w")
    pairs    = open("%s/pairs.top" %save_folder,"w")
    section = None
    topology.write(" ; Pairs section removed by convert_top function \n")
    pairs.write(" #  Pairs parameters derived from %s \n" %raw_topology)
    for line in top_file:
        if re.match('.*\[.*\]',line):
            section = line.split()[1]
        if (section == 'pairs'):
            if  (re.match('.+\d+\s+\d+\s+\d+.+',line)):
                new_line=line.split()
                if  int(new_line[2]) == 6:
                    pairs.write(line)
                else:
                    raise ValueError("Function type %s is not supported in this script" %line.split()[2])
            else:
                pairs.write('# '+line)
        else:
            topology.write(line)
    top_file.close()
    topology.close()
    pairs.close()
    return()

def parse_pairs(input_file, pairwise_params='pairwise_params', model_params='model_params'):
    """
    The function is designed to parse "pairs" section of GROMACS topology file and
    generate pairwise_params file and model_params files in the format, used by pyODEM
    At this point, only Gaussian potential for native contacts is implemented. To be used with 
    al-atom structure-based models 

    parameters
    ----------
    
    input_file : string
    A string with path and name of file, that includes section of gromacs topology with information
    about pairs
    Structure of input file (can be found here: http://smog-server.org/extension/gauss.html
     Column 1: atom i, indexed from 1 in the CG representation
     Column 2: atom j, indexed from 1 in the CG representation
     Column 3: function type (only 6 is allowed in current version)
     Column 4: interaction strength (modified during ODEM optimization)
            5: r0 
            6: sigma
            7: a (excluded volume)
    
    model_params: string, optional 
    A string with a path and name of the output file with pairwise interaction descriptions. 
    File structure (Description of columns 1-5 taken from https://github.com/ClementiGroup/pyODEM/blob/master/examples/example_2/Protein_Optimization.ipynb). In this description, the original notation, different from one described
    above is used. Each line describes one pair.
   
    Column 1: Atom i, indexed from 1 in the CG representation.
    Column 2: Atom j, indexed from 1 in the CG representation.
    Column 3: A count of the number of pairwise function's index (indexing from 0 as it is used internally).
    Column 4: A string giving the type of potential (e.g. LJ1210, LJ126, LJ12GAUSSIAN)
    Column 5-: All columns starting with the fifth column are used to parameterize the pairwise potential.
    For this program, column 4 will contain LJ12GAUSSIAN only. For this type of function, corresponding description
    of parameters shown below:
    Column 5: Gives the excluded valume r0 used in the (r0/r)^12 part of a function. (Connected with column 7 in
              the input file, r0 in output = a^(1/12 ) in input.
Column 6: Gives the minima of the Gaussian well. (The same as column 5 in input file)
Column 7: Gives the standard-width of the Gaussian well. (The same as column 6 in input file)                              
     
   
    pairwise_params: string, optional
    A string with a path and name of the output file with model parameters, which includes model parameters 
    interactin strenght, one numer per line. Order should be the same as in  model_params. 
 

    """

    input_file = open(input_file,"r")
    parameters = open(model_params,"w")
    pairs      = open(pairwise_params,"w")
    counter = 0 # Counts number of parameters written to the file
    for line in input_file:
        if not line.strip().startswith("#"):
            data = line.split()
            params_line = data[3] +'\n'
            if int(data[2])==6:
                function_type = 'LJ12GAUSSIAN'
            else:
                raise(ValueError,'Function type %d is not implemented' %int(data[2]))
            delimiter = ' '
            excluded_volume_new = float(data[6])**(1.0/12.0)
            pairs_line = delimiter.join([data[0],data[1],str(counter),function_type,'%f'%excluded_volume_new,data[4],data[5],'\n'])
            parameters.write(params_line)
            pairs.write(pairs_line)
            counter+=1
    
    input_file.close()
    parameters.close()
    pairs.close()
    return()


if __name__ == '__main__':
    parse_pairs(input_file='pairs.top', pairwise_params='pairwise_params', model_params='model_params')

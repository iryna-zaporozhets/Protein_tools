import numpy as np

class Interaction():
    """ parent class for all the possible interaction
     classes that will be added in future
    """
    def __init__(self):
        pass


class NonBondedInteraction(Interaction):
    """
    Class defines instances of pairvise interactions
    """

    # Parameter dictionary defines a set of parameters in appropriate order for
    # each function type. {Function (str) : [tuple_of_parameter_names (str) in
    #  correct order]}
    # Current version of names in parameter dictionary is based on
    # Chen et al, J. Chem Theory Comput 2018, 14, 3849-3858
    parameter_dictionary  = {'LJ12GAUSSIAN':('excluded_volume','distance','gaussian_width'),
                            'LJ12GAUSSIANTANH':('excluded_volume','distance','offset_sigma_t')}

    def __init__(self,pair,potential_type,strength_value, parameters,
                 parameter_names):
        """
         Initialize an instance of pairwise interactions.

         parameters
         ----------

          pair : list of integers, len 2
          list of two integers, each corresponds to a number of atom in
          potential_type : str
          potential type code, than will then be used by openmm parser
          strength_value : float
          Represents strength of a particular interacton

          parameters  : dictionary.
          Dicrionary of parameters. key should be the parameter name,
          value should be the correspinding value

        """

        self.pair = [min(pair),max(pair)]
        self.potential_type = potential_type
        self.strength_value = strength_value
        self.parameters = parameters
        self.parameter_names = parameter_names

    def __str__(self):
        string  = "Non-bonded interaction between particles %d and %d\nPotential type : %s\n" %(self.pair[0], self.pair[1],self.potential_type)
        return string

    def set_strength_value(self,strength_value):
        self.strength_value = strength_value

    @classmethod
    def read(cls,input_string,strength=None,parameter_dict=parameter_dictionary):
        """
        Allows to create an interaction based on a string input with the following
        format:
        particle_i particle_j number potential_type other_parametsrs
        Strength of the parameters will be supplied separatly
        """
        input_data = input_string.split()
        pair = [int(input_data[0]),int(input_data[1])]
        potential_type = input_data[3]
        parameters = [float(i) for i in input_data[4:]]
        parameter_names = parameter_dict[potential_type]
        return cls(pair,potential_type,strength, parameters,
                     parameter_names)

    def write(self,number):
        """
        The method returns two strings: one contains a strength, another conatains
        all the pairwise information

        parameters
        ----------

        self - instance of a class
        number - number of the interaction in order

        returns
        --------
        string1 - str
             contains strength of interaction
        string2 -str
             contains pairwise parameters
        """
        if len(self.parameters) != 3:
            print("WARDING: only 3 first parameters will be written")
        string1 = str(self.strength_value) + '\n'
        string2 = "%6d  %6d  %12d       %s     %.6f  %.6f  %.6f\n" %(self.pair[0],
                                                                     self.pair[1],
                                                                     number,
                                                                     self.potential_type,
                                                                     self.parameters[0],
                                                                     self.parameters[1],
                                                                     self.parameters[2])
        return (string1, string2)


class ModelParameters():
    """
    Class holds parameters of the model.
    """
    def __init__(self):
        """
        Initizlaize instance of ModelParameters class. Should be a list
        of objects, that define individual model parameters. So far, only
        class NonBondedInteraction has been implemented.
        For each object in a list, the following methods should be available:
        read(self,string) - creates an instance of an Interaction class based
        on a string.
        """
        self.model_params = []

    def add_parameter(self,parameter):
       """
       parameter : object, instance of the class Interaction
       Object, that define an individual parameter of the model
       """
       pair = parameter.pair
       if self.includes_intraction(pair) == -1:
           self.model_params.append(parameter)
       else:
           raise ValueError("Interaction for particles %d and %d already exists in the model"
                            %(pair[0],pair[1]))

    def read(self,input_model_params,input_pairwise_params):
        """
        Reads parameter from the input file.

        Parameters:
        -----------
        input_model_param : str
                           Path to the file with strength of interactions

        input_pairwise_params : str
                              Path to the file with interaction description
        """
        strength_array = np.loadtxt(input_model_params)
        with open(input_pairwise_params) as file:
            for lines in file:
                if lines.find("#") == -1:
                    self.add_parameter(NonBondedInteraction.read(lines,strength_array[int(lines.split()[2])]))

    def write(self,
              model_params ='model_params',
              pairwise_params = 'pairwise_params'):
        """ Method writes files with model parameters and pairwise parameters
        """
        strength = open(model_params,"w")
        properties = open(pairwise_params,"w")
        properties.write('#    pairs         param         potential_type      other_params \n')

        for number, parameters in enumerate(self.model_params):
            string1, string2 = parameters.write(number)
            strength.write(string1)
            properties.write(string2)

        strength.close()
        properties.close()


    def includes_intraction(self,pair):
        """
        Method checks, whether interaction of a particular pair is included in the model.
        If yes, returns index of this parameter,
        If no, returns -1

        parameters
        ----------
        pair - list of two integers
        """
        for interaction in self.model_params:
            if ((interaction.pair[0] == min(pair)) and (interaction.pair[1] == max(pair))):
                return self.model_params.index(interaction)
        return -1

    def modify_strength(self,pair,new_strength):
        """ The method modifies sthength of already existing interactions in the model
        parameters
        ----------
        pair : list of int
        Contains list of two integers, representing atoms.
        """

        index = self.includes_intraction(pair)
        if index != -1:
            parameter = self.model_params[index]
            parameter.set_strength_value(new_strength)
        else:
            raise ValueError("Interaction between particles %d and %d does not exists in the model" %(pair[0],pair[1]))

        return 0

    def get_interaction_strength(self,pair):
        """
        The methods returns interaction strength (epsilon value) of a particular
        interaction

        """

        index = self.includes_intraction(pair)
        if index != -1:
            return self.model_params[index]
        else:
            raise ValueError("Interaction between particles %d and %d does not exists in the model" %(pair[0],pair[1]))





############# Testing area #####################################################

#interaction = NonBondedInteraction([53,101],'LJ12GAUSSIAN',1.0, [1,2,3], NonBondedInteraction.parameter_dictionary['LJ12GAUSSIAN'])
# print(interaction)
#
#
# interaction = NonBondedInteraction.read('1      28             0       LJ12GAUSSIAN     0.266000  0.820535  0.050000')
# print(interaction)
#
# string1, string2 = interaction.write(1)
#
# print(string1)
# print(string2)
#
#model = ModelParameters()
#model.read('model_params','pairwise_params')
#model.add_parameter(interaction)
#model.modify_strength([53,101],1000)
#model.write('model_params_test','pairwise_params_test')

"""
Module contains custom openMM tools:

generate_system_report(system object) : print a detailed human-readable report about the content of system object  
get_force_by_name(<system object>, <force name string>) : returns reference to a force by its name
"""

def get_force_by_name(system,force_name):
    """
    Return a list of references to forces in the system that have a name specified by the force_name
    The force name must be specified exactly as in OpenMM documentation.

    """
    forces = system.getForces()
    forces_to_return = []
    count = 0
    for force in forces:
        if type(force).__name__ == force_name:
            count += 1
            forces_to_return.append(force)
    print("{} forces with name {} were found.".format(count,force_name))
    return forces_to_return

def get_single_force_index_by_name(system, force_name):
    """
    Returns an index of the first force with a specified name
    """
    forces = system.getForces()
    for ndx, force in enumerate(forces):
        if type(force).__name__ == force_name:
            return  ndx
    return None

def _report_NonbondedForce(force):
    """
    Generate a report for NonbondedForce in the system
    """
    fstr="{:30}   {}"
    not_impl = "Inquiry not implemented. See docs for "
    print(fstr.format("Cutoff Distance", force.getCutoffDistance()))
    print(fstr.format("Ewald Error Tolerance", force.getEwaldErrorTolerance()))
    print(fstr.format("Reaction Field Dielectric", force.getReactionFieldDielectric()))
    print(fstr.format("Reciprocal Space Force Group", force.getReciprocalSpaceForceGroup()))
    print(fstr.format("Switching Distance", force.getSwitchingDistance()))
    print(fstr.format("Use Dispersion Correction", force.getUseDispersionCorrection()))
    print(fstr.format("Use Switching Function", force.getUseSwitchingFunction()))
    print(fstr.format("Force Group",force.getForceGroup()))
    print(fstr.format("LJ PME Parameters",force.getLJPMEParameters()))
    print(fstr.format("LJ PME Parameters In Context", "Context-specific. See docs for getLJPMEParametersInContext"))

    print(fstr.format("PME Parameters", force.getPMEParameters()))
    print(fstr.format("PME Parameters In Context", "Context_specific. See docs for getPMEParametersInContext"))

    print(fstr.format("Nonbonded Method", force.getNonbondedMethod()))

    num_exception_parameter_offsets = force.getNumExceptionParameterOffsets()
    print(fstr.format("Num Exception Parameter Offsets",num_exception_parameter_offsets))
    print("offsetIndex    parameter    exceptionIndex    chargeProdScale    sigmaScale   epsilonScale")
    for offset_ndx in range(num_exception_parameter_offsets):
        print("{:11d}    {:9}    {:14d}    {}   {}    {}".format(*force.getExceptionParameterOffset(offset_ndx)))
    print(" -                -                -                 -                -             - ")

    num_particle_parameter_offsets = force.getNumParticleParameterOffsets()
    print(fstr.format("Num Particle Parameter Offsets", num_particle_parameter_offsets ))
    print("offsetIndex    parameter    particleIndex    chargeScale    sigmaScale    epsilonScale")
    for offset_ndx in range(num_particle_parameter_offsets):
        print("{:11d}   {:9}    {:13d}    {.3e}    {}    {}".format(*force.getParticleParameterOffset(offset_ndx)))
    print("    -             -              -                -              -               - ")

    num_exceptions = force.getNumExceptions()
    print(fstr.format("Num Exceptions", num_exceptions ))
    print("index    particle i    particle j    chargeProd    sigma    epsilon")
    for exception_ndx in range(num_exceptions):
        print("{:5d}    {:10d}    {:10d}    {}    {}    {}".format(exception_ndx,*force.getExceptionParameters(exception_ndx)))
    print("  -          -             -             -          -          -")
    
    num_particles = force.getNumParticles()
    print(fstr.format("Num Particles", num_particles))
    print("index    charge    sigma    epsilon")
    for particle_ndx in range(num_particles):
        print("{:5d}    {}    {}    {}".format(particle_ndx, *force.getParticleParameters(particle_ndx)))
    print(" -          -       -       -")
    
    num_global_parameters = force.getNumGlobalParameters()
    print(fstr.format("Num Global Parameters", num_global_parameters))
    print("parameter index    name    default value")
    for global_parameter_ndx in range(num_global_parameters):
        print("{}    {}    {}".format(global_parameter_ndx,
                                     force.getGlobalParameterName(global_parameter_ndx),
                                     force.getGlobalParameterDefaultValue(global_parameter_ndx)
                                     ))
    print("      -              -              -")

    return



def _report_CustomNonbondedForce(force):
    """
    Generate a report for CustomNonbondedForce in the system
    """
    fstr="{:30}   {}"
    print(fstr.format("Cutoff Distance", force.getCutoffDistance()))
    print(fstr.format("Energy Function", force.getEnergyFunction()))
    print(fstr.format("Force Group", force.getForceGroup()))
    print(fstr.format("Use Long Range Correction",force.getUseLongRangeCorrection()))
    print(fstr.format("Use Switching Functioin", force.getUseSwitchingFunction()))
    print(fstr.format("Switching Distance", force.getSwitchingDistance()))
    print(fstr.format("NonbondedMethod",force.getNonbondedMethod()))
    num_energy_parameter_derivatives = force.getNumEnergyParameterDerivatives()
    print(fstr.format("Num Energy Parameter Derivatives", num_energy_parameter_derivatives))
    print("parameter derivative index,   energy parameter derivative name")
    for energy_parameter_derivative_ndx in range(num_energy_parameter_derivatives):
        print("{:26d}    {:26}".format(energy_parameter_derivative_ndx,
                                        force.getEnergyParameterDerivativeName))
    print("              -                              -                 ")
    
    num_particles = force.getNumParticles()
    print(fstr.format("Num Particles", num_particles))
    num_per_particle_parameters = force.getNumPerParticleParameters()
    print(fstr.format("Num Per Particle Parameters", num_per_particle_parameters))
    header_string = 'particle index,  '
    for per_particle_param_ndx in range(num_per_particle_parameters):
        header_string += "parameter {:10}".format(force.getPerParticleParameterName(per_particle_param_ndx))
    print(header_string)
    for particle_ndx in range(num_particles):
        template_string = (num_per_particle_parameters+1)*"      {:10}"
        print (template_string.format(particle_ndx,*force.getParticleParameters(particle_ndx)))
    print("      -      "*(num_per_particle_parameters+1))
     
    num_global_parameters = force.getNumGlobalParameters()
    print(fstr.format("Num Global Parameters", num_global_parameters))
    print("parameter index    name    default value")
    for global_parameter_ndx in range(num_global_parameters):
        print("{}    {}    {}".format(global_parameter_ndx,
                                     force.getGlobalParameterName(global_parameter_ndx),
                                     force.getGlobalParameterDefaultValue(global_parameter_ndx)
                                     ))
    print("      -              -              -")
    
    
    num_exclusions = force.getNumExclusions()
    print(fstr.format("Num Exclusions", num_exclusions ))
    print("index    particle i    particle j")
    for exclusion_ndx in range(num_exclusions):
        print("{:5d}    {:10d}    {:10d}".format(exclusion_ndx,*force.getExclusionParticles(exclusion_ndx)))
    print("  -          -              -     ")

    num_tabulated_functions = force.getNumTabulatedFunctions()
    print(fstr.format("Num Tabulated Functions",num_tabulated_functions))
    print("Tabulated function index,     name")
    for tabulated_function_ndx in range(num_tabulated_functions):
        print("{:25d}    {}".format(tabulated_function_ndx, force.getTabulatedFunctionName()))
    print("                -              - ")
    print("To get references to tabulated function object, use getTabulatedFunction(index)")

    num_interaction_groups = force.getNumInteractionGroups()
    print(fstr.format("Num Interaction Groups", num_interaction_groups))

    for interaction_group_ndx in range(num_interaction_groups):
        set_1, set_2 = force.getInteractionGroupParameters(interaction_group_ndx)
        print ("index:", interaction_group_ndx, "; set_1: ", set_1, "; set_2: ", set_2)
     
    
    return

def _report_HarmonicBondForce(force):
    """
    Generate a report for HarmonicBondForce in the system
    """
    fstr="{:30}   {}"
    print(fstr.format("Force Group", force.getForceGroup()))
    num_bonds = force.getNumBonds()
    print(fstr.format("Num Bonds", num_bonds))
    print("bond index     particle i    particle j    length             k")
    for bond_ndx in range(num_bonds):
        print("{:10d}    {:10d}    {:10d}    {}    {}".format(bond_ndx, *force.getBondParameters(bond_ndx)))

def _report_HarmonicAngleForce(force):
    """
    Generate a report  for HarmonicAngleForce
    """
    fstr="{:30}   {}"
    print(fstr.format("Force Group", force.getForceGroup()))
    num_angles = force.getNumAngles()
    print(fstr.format("Num Angles", num_angles))
    print("angles index     particle i    particle j    particle k    angle                  k")
    for angle_ndx in range(num_angles):
        print("{:10d}    {:10d}    {:10d}    {:10d}    {}    {}".format(angle_ndx, *force.getAngleParameters(angle_ndx)))

def _report_PeriodicTorsionForce(force):
    """
    Generate a report for PeriodicTorsionForce in the system
    """
    fstr="{:30}   {}"
    print(fstr.format("Force Group", force.getForceGroup()))
    num_periodic_torsions = force.getNumTorsions()
    print(fstr.format("Num Torsions", num_periodic_torsions))
    print("torsion index     particle i    particle j    particle k    particle l    periodicity   phase    k")
    for torsion_ndx in range(num_periodic_torsions):
        print("{:10d}    {:10d}    {:10d}    {:10d}    {:10d}       {:10d}    {}  {}".format(torsion_ndx, 
                                                                          *force.getTorsionParameters(torsion_ndx)))

def _report_CustomBondForce(force):
    """
    Generate a report for CustomBondForce in the system
    """
    fstr="{:30}   {}"
    print(fstr.format("Force Group", force.getForceGroup()))
    print(fstr.format("Energy Function", force.getEnergyFunction()))


    num_global_parameters = force.getNumGlobalParameters()
    print(fstr.format("Num Global Parameters", num_global_parameters))
    print("parameter index    name    default value")
    for global_parameter_ndx in range(num_global_parameters):
        print("{}    {}    {}".format(global_parameter_ndx,
                                     force.getGlobalParameterName(global_parameter_ndx),
                                     force.getGlobalParameterDefaultValue(global_parameter_ndx)
                                     ))
    print("      -              -              -")
    
    
    num_bonds = force.getNumBonds()
    print(fstr.format("Num Bonds", num_bonds))
    num_per_bond_parameters = force.getNumPerBondParameters()
    print(fstr.format("Num Per Bond Parameters", num_per_bond_parameters))
    header_string = 'bond index,   particle i,    particle j,      '
    for per_bond_param_ndx in range(num_per_bond_parameters):
        header_string += "{:10}".format(force.getPerBondParameterName(per_bond_param_ndx))
    print(header_string)
    for bond_ndx in range(num_bonds):
        template_string = (num_per_bond_parameters+3)*"         {}"
        particle_i,particle_j,params = force.getBondParameters(bond_ndx)
        print (template_string.format(bond_ndx,particle_i, particle_j, *params))
    print("      -      "*(num_per_bond_parameters+1))

    
    num_energy_parameter_derivatives = force.getNumEnergyParameterDerivatives()
    print(fstr.format("Num Energy Parameter Derivatives", num_energy_parameter_derivatives))
    print("parameter derivative index,   energy parameter derivative name")
    for energy_parameter_derivative_ndx in range(num_energy_parameter_derivatives):
        print("{:26d}    {:26}".format(energy_parameter_derivative_ndx,
                                        force.getEnergyParameterDerivativeName))
    print("              -                              -                 ")


def _report_other(force):
    """
    Generate a report for forces, for which custom generation is not supported
    """
    for method in dir(force):
        if method[0:3]=='get':
            try:
                print("{:30}   {}".format(method[3:],getattr(force,method)()))
            except:
                print("{:30} : Automatic inquiry is not yet implemented.  See documentation for  {}  method usage".format(method[3:],method))
    return()

def generate_system_report(system):
    """
    Generate detailed report about the openMM system: 

    input: openMM system object
    """
    number_of_particles = system.getNumParticles()
    number_of_constraints = system.getNumConstraints()


    print("==================openMM system properties========================")
    print("  ")
    
    print("PARTICLES")
    print(" ")
    print("Total number of particles: {}".format(number_of_particles))
    print ("Particle index        Particle mass       Is virtual site")
    for particle_ndx in range(number_of_particles):
        print("{:14d}        {:19}       {:15}".format(particle_ndx,   
                                                      str(system.getParticleMass(particle_ndx)),
                                                      str(system.isVirtualSite(particle_ndx))
                                                    ))
    print(66*"-")
    print(" ")

    print("CONSTRAINTS")
    print(" ")
    print("Total number of distance constraints: {}".format(number_of_constraints))
    print("particle i        particle j        distance")        
    for constraint_ndx in range(number_of_constraints):
        print("{}        {}        {}".format(*system.getConstraintParameters()))
    print(66*"-")
    print(" ")

    print("PERIODIC BOUNDARY CONDITIONS")
    print(" ")
    print("Periodic boundary conditions enabled:  {}".format(system.usesPeriodicBoundaryConditions()))
    print("Periodic box vectors:\n a {}; \n b {}; \n c {};".format(*system.getDefaultPeriodicBoxVectors()))
    print(66*("-"))
    print(" ")

    # ---------------------Forces ------------------------------------------------------
    print("FORCES")
    print(" ")
    force_dictionary={"NonbondedForce" : _report_NonbondedForce,
                      "CustomNonbondedForce" : _report_CustomNonbondedForce,
                      "HarmonicBondForce" : _report_HarmonicBondForce,
                      "HarmonicAngleForce": _report_HarmonicAngleForce,
                      "PeriodicTorsionForce" : _report_PeriodicTorsionForce,
                      "CustomBondForce" : _report_CustomBondForce
                      }
    number_of_forces = system.getNumForces()
    print("Number of Force objects in the system: {}".format(number_of_forces))
    print("")
    for force in system.getForces():
        force_name = type(force).__name__
        print("-->  {}".format(force_name))
        print(" ")
        if force_name in force_dictionary:
            report_function = force_dictionary[force_name]
        else:
            report_function = _report_other
        report_function(force)
        print(67*"_")
    print(" END OF THE REPORT ")
    return






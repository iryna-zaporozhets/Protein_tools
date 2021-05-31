"""
A set of utilits for analyzing outcome of odem optimization
"""
import datetime
import matplotlib.pyplot as plt
import matplotlib
import os
import datetime
import numpy as np
import pandas as pd
import mdtraj as md


def parse_info_file(filename):
    """
    Function parces info file with string name filename
    and returns temperature, old and new Q0-value
    Assumes, that infofile has 3 lines in order T, oldQ, newQ.
    """
    file = open(filename,'rt')

    line = file.readline()
    data = line.split(':')
    assert 'temperature' in data[0], 'Wrong file layout: check temperature field'
    T = float(data[1])

    line = file.readline()
    data = line.split(':')
    assert 'Qold' in data[0], 'Wrong file layout: check Qold field'
    oldQ = float(data[1])

    line = file.readline()
    data = line.split(':')
    assert 'Qnew' in data[0], 'Wrong file layout: check Qnew field'
    newQ = float(data[1])
    return T, oldQ, newQ

def generate_metadata(title, author, description=None, enable_warning=False):
    """
    The function generates metadata for saving a file in png format.
    It pastes current working directory as a source of data, writes
    data and time of plot generation,

    Parameters:
    ----------
    title : str
            Short (one line) title or caption for image
    author: str
            Author
    description: str
            Description of image (possibly long)
    comments : str
             Any miscellaneous comments

    Dependencies
    ------------
    os
    datetime

    Returns:
    -------
    metadata  : dict
    """
    metadata = {'Title':  title,
                'Author': author
                }
    if description is not None:
        metadata['Description'] = description
    if enable_warning:
        metadata['Warning'] = 'Caution: science!'

    # Getting time
    now = datetime.datetime.now()
    string_time = now.strftime("%Y-%m-%d %H:%M")
    metadata['Creation Time'] = string_time

    metadata['Comment']= "Path to the source: {}".format(os.getcwd())
    return metadata


def find_max_iteration(folder,
                       subfolder_root='newton',
                       limit=10000,
                       start_iteration=1,
                       check_continuity=True):
    """
    The function finds number of finished iterations, in a particular folder.
    An iteration with index `ndx` is considered finished, if a folder with name
    <subfolder_root>_<iteration ndx> exists.
    Search for iteration will start from `start_iteration` and will be limited by `limit`
    The highest found ndx is returned as number of iterations.
    If check continuity is True, function checks that all the iterations in range
    [start_iteration ..max_iteration] exist. If any iteration is missing, an AssertionError
    is raised
    """
    # Determine max present iteration
    max_iteration = 0
    for i in range(start_iteration, limit):
        if os.path.isdir("{}/{}_{}".format(folder, subfolder_root, i)):
            max_iteration = i

    # If check_continuity, determine if all the iterations
    # in range start_iteration, max_iteration, including max iteration, exist.
    if check_continuity:
        for i in range(start_iteration, max_iteration+1):
            assert os.path.isdir("{}/{}_{}".format(folder, subfolder_root, i)), 'iteration {} does not exists'.format(i)
    return (max_iteration)


def get_temperatures(folder, info_file_name='info.txt',subfolder_root='newton'):
    """
    The function returns a numpy array, where element with index [i]
    stores temperature, at which iteration [i] was performed.
    """
    max_iter = find_max_iteration(folder,subfolder_root=subfolder_root,limit=10000)
    temperatures = []
    Q_olds = []
    Q_news = []
    for iteration_ndx in range(1,max_iter+1):
        T, Q_old,Q_new = parse_info_file('{}/{}_{}/{}'.format(folder,subfolder_root,iteration_ndx,info_file_name))
        temperatures.append(T)
        Q_olds.append(Q_old)
        Q_news.append(Q_new)
    return temperatures, Q_olds, Q_news


def plot_modeling_temperatures(folder_list, ax=None, label_list=None):
    """
    Function plots temperatures, at which each iteration was performed.
    Params:
    folder_list: list of folders to analyze.
    ax : Matplotlib axis obect to plot on
    label_list: list of labels. If none, folders will be used as labels

    """
    if ax is None:
        fig, ax = plt.subplots()
    if label_list is None:
        label_list = folder_list

    # Folding temerature for different runs

    for folder,label in zip(folder_list, label_list):
        temperatures, _, _ = get_temperatures(folder)
        ax.plot(temperatures,label=label)
        ax.set_xlabel("iteration")
        ax.set_ylabel("T, model units")
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
              ncol=3, mode="expand", borderaxespad=0.)
    print("Graph of modeling temperature as a function of iteration number.")
    return ax


def plot_heat_capacity(folder, ax=None, labeling = 'cb', cmap_name ='viridis', cmap=None):
    """
    Functoin loads and plots heat capacity for an odem cycle stored in a particular folder
    labeling : {'cb', 'legend'}
               Defines, how curves will be mapped to iterations. If labeling='cb', colorbar
               will be used with colormap `cmap`.
               If labeling='legend', a legend will be used with labeled with numbers.
    """
    max_iter = find_max_iteration(folder)
    if cmap == None:
        cmap = matplotlib.cm.get_cmap(cmap_name)
    if ax is None:
        if labeling == 'cb':
            with plt.rc_context({'axes.prop_cycle':  matplotlib.cycler(color=cmap(np.linspace(0, 1, max_iter)))}):
                fig, ax = plt.subplots()
        else:
            fig, ax = plt.subplots()


    if labeling == 'cb':
        bounds = list(range(1, max_iter+1))
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
        cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
        cb.set_label("Iteration")

    for iter_ndx in range(1,max_iter+1):
        cv =  list(np.loadtxt("{}/Cv_iteration_{}.txt".format(folder,iter_ndx)))
        cv.sort(key=lambda x: x[0])
        cv = np.array(cv)
        ax.scatter(cv[:,0],cv[:,1],label=' {}'.format(iter_ndx))
        ax.plot(cv[:,0],cv[:,1])

    if labeling == 'label':
        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
          ncol=4, mode="expand", borderaxespad=0.)

    ax.set_xlabel('Temperature, model units')
    ax.set_ylabel("Cv, kJ/(mol*K)")
    return ax


def plot_Q_function(folder_list, ax=None, label_list=None, Q_type='old'):
    """
    Function plots Q function, that corresponds to a particular iteration
    Params:
    folder_list: list of folders to analyze.
    ax : Matplotlib axis obect to plot on
    label_list: list of labels. If none, folders will be used as labels
    Q_type : str, {'old', 'new'}
             Q_type='old' Q-value calculated directly from model parameters
                    and the trajectory, calculated with these model parameters.
                    Gives a real representation of how good this model is
             Q_type='new' Q-value, obtained with the same trajectory, as in case
                     of Q_type='old', but with model parameters modified with odem
                     procedure. Cannot be used to assess model, but may give insights
                     into optimization algorithm behavior

    """
    if ax is None:
        fig, ax = plt.subplots()
    if label_list is None:
        label_list = folder_list

    # Folding temerature for different runs

    for folder,label in zip(folder_list, label_list):
        _, Q_old,Q_new = get_temperatures(folder)
        if Q_type == 'old':
            ax.plot(Q_old, label=label)
        elif Q_type == 'new':
            ax.plot(Q_new, label=label)
        ax.set_xlabel("iteration")
        ax.set_ylabel("-ln Q")
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
              ncol=3, mode="expand", borderaxespad=0.)
    print("Graph of -ln Q ({}) as a function of number of interations.".format(Q_type))
    return ax


def rmse(target,prediction):
    """
    Return rmse between target and prediction
    """
    return np.sqrt(np.mean((prediction-target)**2))


def get_run_rmse(folder, experiment, prediction_root_name='ddG_calculated'):
    """
    For each iterations in folder, calculates RMSE with respect to experimental data `experiment`
    Parameters:
    -----------
    folder : str
             folder, were all the prediction files are located
    experiment : numpy 1D array
                array with experimental data
    prediction_root_name : str
                           It is assumed, that prediction files have name
                           <predictoin_root_name>_<iteration number>.txt

    Returns:
    --------
    rmse_list : list of rmse, 1 value per iteration

    """
    rmse_list = []
    max_iter = find_max_iteration(folder)
    for iter_ndx in range(1,max_iter+1):
        prediction = np.loadtxt("{}/{}_{}.txt".format(folder, prediction_root_name, iter_ndx))
        rmse_list.append(rmse(experiment,prediction))
    return rmse_list


def plot_rmse(rmse_list, units, ax=None, label_list=None):
    """
    Function plots rmse as a function of iteration index for each run folder in folder_list.
    Params:

    rmse_list : list of length N, where N - number of runs
                Each element of rmse_list is another list. The latter contains one element
    units      : str
                Units of data. Is displayed on a plot
    ax : Matplotlib axis obect to plot on
    label_list: list of labels. If none, folders will be used as labels

    """
    if ax is None:
        fig, ax = plt.subplots()
    if label_list is None:
        label_list = ['' for run_rmse in rmse_list]

    # Folding temerature for different runs
    for run_rmse, label in zip(rmse_list, label_list):
        ax.plot(run_rmse,label=label)
        ax.set_xlabel("iteration")
        ax.set_ylabel("RMSE, {}".format(units))
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
              ncol=3, mode="expand", borderaxespad=0.)
    print("Graph of modeling temperature as a function of iteration number.")
    return ax


def get_rmse_summary(rmse_list, labels=None):
    """
    Creates a pandas dataframe with rmse statistics for each run.
    """
    # Run summary:
    num_of_iters = []
    initial_rmse = []
    final_rmse = []
    min_rmse = []
    iter_min_rmse = []
    for run_rmse in rmse_list:
        num_of_iters.append(len(run_rmse))
        initial_rmse.append(run_rmse[0])
        final_rmse.append(run_rmse[-1])
        min_rmse.append(np.min(run_rmse))
        iter_min_rmse.append(np.argmin(run_rmse)+1)

    run_summary = pd.DataFrame(data={"Number_of_iterations":num_of_iters,
                                     "Initial RMSE": initial_rmse,
                                     "Final RMSE": final_rmse,
                                     "Minimum RMSE" : min_rmse,
                                     "Minimum RMSE at iteration" : iter_min_rmse,
                                }
                          )
    run_summary['delta RMSE'] = run_summary['Minimum RMSE'] - run_summary['Initial RMSE']
    if labels is not None:
        run_summary['labels'] = labels
    return run_summary


def get_observable_values(folder_list, target_iterations, prediction_root_name='ddG_calculated'):
    """
    Get observable values for each each folder in a folder_list. The iteration after which observable
    values were calculated, is specified with target_iterations.

    Parameters
    ----------
    folder_list : list of str
                  list of run folders
    target_iterations : list of ints
                  list of iteration numbers to get observable values for.
                  should be in the same order as folder_list.
  prediction_root_name : str
                    File, that contains calculated values.
                   It is assumed, that prediction files have name
                   <predictoin_root_name>_<iteration number>.txt

    Returns:
    --------

    obs_list : list of 1d numpy arrays
               List contain 1 numpy array per folder in folder_list
               Each element corresponds to calculated observable values obtained at a particular iteration,
               specified in target_iterations

    """
    obs_list = []
    assert len(folder_list) == len(target_iterations), "Number of run folders does not match number of target iterations"
    for folder, iteration in zip(folder_list, target_iterations):
        prediction = np.loadtxt("{}/{}_{}.txt".format(folder, prediction_root_name, iteration))
        obs_list.append(prediction)
    return obs_list


def plot_correlation(obs_list,
                     reference,
                     ax = None,
                     mode='cumulative',
                     plot_std='True',
                     annotate_points=False,
                     point_annotation = None,
                     label_list=None,
                     xerr=None,
                     yerr=None,
                     obs_name= r'$\frac{\Delta \Delta G}{RT}$',
                     errorbar_kwargs={'fmt':'o'},
                     annotation_kwkargs={},
                     arrowprops=dict(arrowstyle="-", color='r', lw=0.5)):

    """
    Create a correlation plot based on data in obs_list

    Parameters
    ----------
    obs_list : list of 1D numpy arrays
               list of arrays, Each array contains calculated data
    reference : 1D numpy array
                Reference data to correlate observables in obs_list with. Should be
                in the same order as in arrays in obs_list
    mode: str, {'all','cumulative'}
          Calculation mode.
          'all' : Each list in obs_list is plotted separatly, with labels specified by labels.
                  Errorbars are specified by xerrs or yerrs. Annotation option is unavailable in
                  this mode
          'cumulative' : Average values are plotted, if plot_std is True, y errorbars will be
                         set to sandard deviation, and any specified yerr will be overrided.
    """
    if ax is None:
        fig, ax = plt.subplots()

    if mode == 'all':
        if label_list is not None:
            assert len(obs_list) == len(label_list), "Length of obs_list does not match length of label_list"
        for obs, label in zip(obs_list,label_list):
            ax.errorbar(reference, obs, yerr=yerr, xerr=xerr, label=label, **errorbar_kwargs)
    elif mode == 'cumulative':
            obs_array = np.array(obs_list)
            mean = np.mean(obs_array, axis=0)
            if  plot_std:
                yerr = np.std(obs_array, axis=0)
            if label_list is None:
                label_list = ['']
            ax.errorbar(reference, mean, yerr=yerr, xerr=xerr, label=label_list[0], **errorbar_kwargs)
            if annotate_points:
                assert len(point_annotation)
                from adjustText import adjust_text
                point_annotations = []
                for exp, mean, mutant in zip(experiment, means_initial, point_annotation):
                    point_annotations.append(plt.text(exp, mean, mutant,**annotation_kwkargs))
                adjust_text(point_annotation,arrowprops)
    xlabel = obs_name + ' experimental'
    ylabel = obs_name + ' calculated'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()
    return ax


def  get_model_params(folder_list,
                     target_iterations,
                     param_folder_root_name='newton',
                     param_file_name='model_params',
                     weighted=False,
                     info_file_name='info.txt'):
    """
    Get model parameters, used for performing a particular iteration, specified by `target iteration`
    for each folder in a folder_list (One set for each folder).

    Parameters:
    -----------
    folder_list : list of str
                  list of run folders
    target_iterations : list of ints
                  list of iteration numbers to get model parameters values for.
                  should be in the same order as folder_list. For folder folder_list[i],
                  model parameters used to perform modeling at iteration target_iterations[i]
                  will be retrieved.

    param_folder_root_name : str
                It is assumed, that model params files are located at
               <param_folder_root_name>_<iteration number -1 >/<param_file_name>
               Notice, that for iteration <iteration number> corresponding parameters are retrieved
               from <iteration number -1 > folder.
    param_file_name : str
                Name of the file, containing model parameters. See param_folder_root_name description
                for more details
    weighted : bool
               If True, model parameter are devided by modeling temperature.
    info_file_name : str
               File that contains information about modeling. It is assumed to be located at
               <param_folder_root_name>_<iteration number>/<info_file_name>

    Returns:
    --------
    param_list : list of 1D numpy arrays
                List contain 1 numpy array per folder in folder_list
                Each element corresponds to model parameters used to perform modeling
                at a particular iteration, specified in target_iterations

    """
    param_list = []
    assert len(folder_list) == len(target_iterations), "Number of run folders does not match number of target iterations"
    for folder, iteration in zip(folder_list, target_iterations):
        params = np.loadtxt("{}/{}_{}/{}".format(folder, param_folder_root_name, iteration-1, param_file_name))
        if weighted:
            T, _,_ = parse_info_file('{}/{}_{}/{}'.format(folder,param_folder_root_name,iteration,info_file_name))
            params_weighted = params/T
            param_list.append(params_weighted)

        else:
            param_list.append(params)
    return param_list


def get_native_contact_ndx(reference_file):
    """
    The function returns indexes of native contacts in the array of
    all contacts.  If a value in a reference file is 1, the contact is condidered
    to be native, and nonnative otherwise
    ASSUMPTION: ORDER OF CONTACTS IS THE SAME FOR ALL THE MODELINGS!
    """

    params = np.loadtxt(reference_file)
    native_contact_ndx = np.argwhere(params > 0.999999).T[0]
    return(native_contact_ndx)

def get_contact_probability(traj,frames,pairs,cut_off_distance_nm):
    """
    Calculate contact probatility, that is defined as a fraction
    of frames, where distance between atoms defined in pairs is less
    that cut_off_value

    Parameters:
    traj : mdtraj traj object
           Trajectory to extract frames from

    frames : list of int
            Indexes of frames to be used in analysis.
            0-based indexing should be used

    pairs : list of iterables
          List of pairs of atoms. Numbering should match trajectory.
          (generally 0-based) Each element of the list should have 2
          elements.

    cut_off_distance_nm : float
                          Value in nm. If  distance between two atoms, specified in `pairs`, is
                          less than cut_off_distance_nm, contact considered formed.
    """
    target_ensemble = traj[list(frames)]
    distances = md.compute_distances(target_ensemble, pairs)
    contacts = np.where(distances <cut_off_distance_nm, 1,0)
    probability = np.mean(contacts,axis=0)
    return(probability)

def get_frame_ndx_by_committor(dtrajs, committor, limits):
    """
    Returns frame indexes, that have committor values in the range specified by limits.

    Parameters:
    -----------

    dtrajs : str or 1D numpy array
             Path to a file or 1D numpy array that
             contains a microstate index (0-based) for each frame

    committor : str or 1D numpy array
                Path to a file or 1D numpy array that
                contains committor values for each microstate. Microstates should be
                0-based and correspond to microstates in dtrajs.

    limits : iterable with two floats
             Lower and upper limit of committor to include. Both lower and upper limits are
             included

    Returns:
    --------

    target_ndx : 1D numpy array
                 Indexes of frames which have committor  value <= limits[1] and >= limits[2]
    """

    if isinstance(dtrajs , str):
        dtrajs_arr = np.loadtxt(dtrajs,dtype=int)
    else:
        dtrajs_arr = dtrajs

    if isinstance(committor,str):
        committor_arr = np.loadtxt(committor)
    else:
        committor_arr = committor

    target_microstates = np.where((committor_arr >= limits[0]) & (committor_arr <= limits[1]))
    target_ndx = np.where(np.isin(dtrajs_arr, target_microstates[0]))[0]
    return target_ndx

def float_to_path(number, dlm='.'):
    """
    Convert a floating-point number to a string,
    that can be used in file name or path.
    Trailing 0th are removed.
    """
    number_str = str(number)
    number_parts = number_str.split(dlm)
    assert len(number_parts) <= 2, "More than 1 delimiter found"
    # check decimal part for nonzero characters:
    has_decimal = False
    for char in number_parts[1]:
        if char != '0':
            has_decimal = True
    if has_decimal:
        return number_parts[0] + '_' + number_parts[1]
    else:
        return number_parts[0]

def get_contact_probability_folder(folder,
                                   iteration,
                                   committor_limits,
                                   pairs,
                                   cut_off_distance_nm=0.8,
                                   traj_extension='xtc'
                                   ):
    """
    Find probability that atom pairs listed in `pairs` are in contact in ensemble
    defined by `committor_limits`. Locations of trajectory file and clusterization
    files are determined based on "standart" layout of odem.
    TO DO: make layout changable

    folder : str
             odem run folder
    iteration : int
               iteration number to get results for.
    committor_limits : iterable with two floats
    pairs : iterable
            Iterables contains pairs of atoms for which do calculations
    cut_off_distance_nm : float
            Two atoms are considered in contact, if distance between them is less than cut_off_distance_nm.
            Units should be nm.
    traj_extension : str, default xtc
            trajectory extension


    ODEM LAYOUT SPECIFICATIONS:
    file                :   path
    info file           :  folder/newton_<iteration ndx>/info.txt
    trajectory          :  folder/iteration_<iteration ndx>/traj.<traj_extension>
    topology            :  folder/ref.pdb
    committor file      :  folder/discretization_<iteration ndx>/commitor.txt
    discrite trajectory :  folder/discretization_<iteration ndx>/dtrajs.txt



    """

    temperature, _, _ =  parse_info_file('{}/newton_{}/info.txt'.format(folder,iteration))
    traj = md.load('{}/iteration_{}/{}/traj.{}'.format(folder,iteration, float_to_path(temperature), traj_extension),
     top='{}/ref.pdb'.format(folder))
    committor_file = '{}/discretization_{}/{}/commitor.txt'.format(folder, iteration, float_to_path(temperature))
    dtrajs_file = '{}/discretization_{}/{}/dtrajs.txt'.format(folder, iteration,  float_to_path(temperature))
    ensemble_frames = get_frame_ndx_by_committor(dtrajs_file, committor_file, committor_limits)
    num_ensemble_frames = len(ensemble_frames)
    probability = get_contact_probability(traj,ensemble_frames, pairs,cut_off_distance_nm)

    return probability, num_ensemble_frames


def get_contact_probabilities_folder_list(folder_list, iteration_list, committor_limits, pairs, cut_off_distance_nm, traj_extension='dcd'):
    """
    Get contact probabilities and number of frames in ensembles for all the trajectories matching folder_list
    and iteration_list
    """
    probabilities = []
    num_ensemble_frames = []
    for folder, iteration in zip(folder_list, iteration_list):
        probability, frames = get_contact_probability_folder(folder,
                                                                   iteration=iteration,
                                                                   committor_limits=committor_limits,
                                                                   pairs=pairs,
                                                                   cut_off_distance_nm=cut_off_distance_nm,
                                                                   traj_extension=traj_extension
                                                                   )
        probabilities.append(probability)
        num_ensemble_frames.append(frames)
    return probabilities, num_ensemble_frames


def get_overall_probability(probabilities, num_ensemble_frames):
    overall_probability = np.zeros(len(probabilities[0]))
    total_frame_count = 0
    for sample, frame_number in zip(probabilities, num_ensemble_frames):
        total_frame_count += frame_number
        overall_probability += sample*frame_number
    overall_probability /= total_frame_count
    return(overall_probability)

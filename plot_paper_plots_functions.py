import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import scipy.stats as stats
#from scipy import ndimage
import os as os
from plot_paper_object import *
from pyFrustration.util import determine_close_residues_from_file
import argparse
"""
Execute All Plotting Macros in __main__
"""

def determine_zero_frustration_result(native_energy_loaded, decoy_avg_loaded, decoy_sd_loaded):
    avg_results = stats.mode(decoy_avg_loaded[np.where(decoy_avg_loaded != 0)], axis=None)[0][0]
    sd_results = stats.mode(decoy_sd_loaded[np.where(decoy_avg_loaded != 0)], axis=None)[0][0]

    zscore_shift = avg_results / sd_results

    return zscore_shift


def calculate_z_distribution(all_z_values):
    # all_z_value is a list of non-zero z-scores
    n_positive = 0
    n_mid_positive = 0
    n_negative = 0
    n_mid_negative = 0

    for zval in all_z_values:
        if zval > 1:
            n_positive += 1
        elif zval > 0:
            n_mid_positive += 1
        elif zval >= -1:
            n_mid_negative += 1
        else:
            n_negative += 1

    return n_positive, n_mid_positive, n_negative, n_mid_negative


def assert_all_zero(mat):
    biggest = np.shape(mat)[0]
    for i in range(biggest):
        for j in range(i+1, biggest):
            assert mat[j,i] == 0

def get_fip35_awsem_strings():
    main_dir = "awsem_frustratometer/FrustrationData"
    prepend_string = "%s/201821420102721686.pdb_" % main_dir

    return main_dir, prepend_string

def get_alpha3d_awsem_strings():
    main_dir = "awsem_frustration_alpha3D/FrustrationData"
    prepend_string = "%s/201832913444214843.pdb_" % main_dir

    return main_dir, prepend_string

def load_frustration_results(main_dir, ndim=2, rescale_native=False, min_native=0.0000):
    decoy_avg = np.loadtxt("%s/decoy_avg.dat" % main_dir)
    decoy_sd = np.loadtxt("%s/decoy_sd.dat" % main_dir)
    if ndim == 2:
        native_energy = np.loadtxt("%s/native_pairwise.dat" % main_dir)
    elif ndim == 1:
        native_energy = np.loadtxt("%s/native_single_residue.dat" % main_dir)

    if rescale_native:
        native_scale = np.loadtxt("%s/native_close_contacts.dat" % main_dir).astype(float)
        assert native_energy.ndim == native_scale.ndim

        if ndim == 2:
            for idx in range(np.shape(native_scale)[0]):
                for jdx in range(np.shape(native_scale)[1]):
                    if (native_scale[idx,jdx] != 0) and (native_scale[idx,jdx] >= min_native):
                        native_energy[idx,jdx] /= native_scale[idx,jdx]
        else:
            pass

    return_dict = {"decoy_avg":decoy_avg, "decoy_sd":decoy_sd, "native_energy":native_energy}

    return return_dict

def compute_frustration(native_energy_loaded, decoy_avg_loaded, decoy_sd_loaded, min_value=0.1, high_bound=None, low_bound=None, zero_bounds=None, nresidues=35, apply_shift=False, shift_value=None, min_res_spacing=2):
    native_energy = np.copy(native_energy_loaded)
    decoy_avg = np.copy(decoy_avg_loaded)
    decoy_sd = np.copy(decoy_sd_loaded)
    #min_value =  np.min(decoy_sd[np.where(decoy_sd > 0)])
    decoy_sd[np.where(decoy_sd < min_value)] = min_value
    #assert_all_zero(decoy_avg)
    #assert_all_zero(native_energy)
    #assert_all_zero(decoy_sd)

    zscore = (decoy_avg - native_energy) / decoy_sd
    #assert_all_zero(zscore)
    for i_zs in range(nresidues):
        for j_zs in range(nresidues):
            if np.abs(i_zs - j_zs) <= 2:
                zscore[i_zs, j_zs] = 0

    if apply_shift:
        mode_results = stats.mode(zscore, axis=None)
        if shift_value is None:
            mode_value = mode_results[0][0]
            mode_count = mode_results[1][0]
            shift_value = mode_value
        shift_indices = np.where(zscore != 0)
        zscore[shift_indices] = zscore[shift_indices] - shift_value

    # set all just-off diagonals to zero
    if zscore.ndim == 2:
        for idx in range(nresidues):
            for jdx in range(nresidues):
                if np.abs(jdx - idx) <= min_res_spacing:
                    zscore[idx, jdx] = 0

    # low_bound and high_bound set if not None:
    if high_bound is not None:
        zscore[np.where(zscore > high_bound)] = high_bound
    if low_bound is not None:
        zscore[np.where(zscore < low_bound)] = low_bound
    if zero_bounds is not None:
        zscore[np.where(np.logical_and(zscore<zero_bounds[1], zscore>zero_bounds[0]))] = 0

    max_native_idxs = np.where(np.abs(native_energy) == np.max(np.abs(native_energy)))
    max_native = native_energy[max_native_idxs]

    max_decoy_idxs = np.where(np.abs(decoy_avg) == np.max(np.abs(decoy_avg)))
    max_decoy = decoy_avg[max_decoy_idxs]

    print "max native energy: %s" % str(max_native)
    print "max decoy energy: %s" % str(max_decoy)

    print "max sd: %f    min sd: %f" % (np.max(decoy_sd), np.min(decoy_sd))

    print "max zscore: %f    min zscore: %f" % (np.max(zscore), np.min(zscore))
    print "max zscore and min zscore locations:"
    print np.where(zscore == np.max(zscore))
    print np.where(zscore == np.min(zscore))

    return zscore, shift_value

def compute_cg_frustration(main_dir):
    fold_e = np.loadtxt("%s/fold_avg.dat" % main_dir)
    other_e = np.loadtxt("%s/other_avg.dat" % main_dir)
    other_std = np.loadtxt("%s/other_std.dat" % main_dir)
    min_value =  np.min(other_std[np.where(other_std >0)])
    other_std[np.where(other_std==0)] = min_value
    zscore = (other_e - fold_e) / other_std

    print "E terms:"
    print np.max(zscore)
    print np.min(zscore)

    return zscore

def load_awsem_frustration(main_dir, prepend_string, nresidues=35):

    configurational = np.loadtxt("%sconfigurational" % prepend_string, comments="#", usecols=(0, 1, 8, 9 ,10))
    mutational = np.loadtxt("%smutational" % prepend_string, comments="#", usecols=(0, 1, 8, 9 ,10))

    config_matrix = calculate_awsem_matrix(configurational, nresidues)
    mutat_matrix = calculate_awsem_matrix(mutational, nresidues)

    return config_matrix, mutat_matrix

def load_awsem_single_residue(main_dir, prepend_string, nresidues=35):

    singleres = np.loadtxt("%ssingleresidue" % prepend_string, comments="#", usecols=(0, 7))

    residues_oneindexed = singleres[:,0]
    frustrations = singleres[:,1]
    return residues_oneindexed, frustrations

def calculate_awsem_matrix(data, nresidues=35):
    zscore = np.zeros((nresidues,nresidues))
    for row in range(np.shape(data)[0]):
        idx = data[row,0] - 1
        jdx = data[row,1] - 1
        zscore[idx, jdx] = (-data[row,2] + data[row,3]) / data[row,4]
        zscore[jdx, idx] = zscore[idx, jdx]

    print "For AWSEM"
    print np.max(np.abs(zscore))
    return zscore


def run_plot_awsem_frustration(plot_obj, native_idxs, dca_idxs, nresidues=35):
    config_matrix, mutat_matrix = load_awsem_frustration(nresidues=nresidues)
    edges = np.arange(0.5, 36, 1.)
    maxplot = 2.5
    plot_obj.plot_epsilon_map([config_matrix], [edges], native_idxs, [""], "configurational_pair_frustration", maxval=maxplot, dca_pairs=dca_idxs, ncols=1, nrows=1, label_columns=False, grid_color="black", dca_color="black")

    plot_obj.plot_epsilon_map([mutat_matrix], [edges], native_idxs, [""], "mutational_pair_frustration", maxval=maxplot, dca_pairs=dca_idxs, ncols=1, nrows=1, label_columns=False, grid_color="black", dca_color="black")

def run_plot_cg_frustration(plot_obj, native_idxs, dca_idxs):

    optimized_frustration = compute_cg_frustration("optimized_data")
    initial_frustration = compute_cg_frustration("initial_data")

    edges = np.arange(0.5, 36, 1.)
    maxplot = 10
    plot_obj.plot_epsilon_map([initial_frustration], [edges], native_idxs, [""], "initial_pair_frustration", maxval=maxplot, dca_pairs=dca_idxs, ncols=1, nrows=1, label_columns=False, grid_color="black", dca_color="black")

    plot_obj.plot_epsilon_map([optimized_frustration], [edges], native_idxs, [""], "optimized_pair_frustration", maxval=maxplot, dca_pairs=dca_idxs, ncols=1, nrows=1, label_columns=False, grid_color="black", dca_color="black")

def analyze_directory(main_dir, plot_obj, native_idxs, dca_idxs, save_name, maxplot=None, reference_map=None, min_value=0.1, compute_low_bound=None, compute_high_bound=None, low_bound=-5, high_bound=5, zero_bounds=[-1, 0.78], nresidues=35, apply_shift=False, shift_value=None, dimensions=2, mark_values=None, rescale_reference=True, min_use=100, min_res_spacing=2, rescale_native=False, min_native=0.0000):
    print "Analyzing %s" % main_dir
    edges = np.arange(0.5, nresidues+1, 1.)
    aa_results = load_frustration_results(main_dir, ndim=dimensions, rescale_native=rescale_native, min_native=min_native)
    aa_frustration, shift_value = compute_frustration(aa_results["native_energy"], aa_results["decoy_avg"], aa_results["decoy_sd"], min_value=min_value, low_bound=compute_low_bound, high_bound=compute_high_bound, zero_bounds=zero_bounds, nresidues=nresidues, apply_shift=apply_shift, shift_value=shift_value, min_res_spacing=min_res_spacing)

    zscores_flat = np.copy(aa_frustration)
    zscores_flat = zscores_flat.flatten()

    sd_spread = np.copy(aa_results["decoy_sd"])
    sd_spread_flat = sd_spread.flatten()
    nonzero_indices = np.where(sd_spread_flat != 0)
    sd_spread_flat = sd_spread_flat[nonzero_indices]
    if np.shape(sd_spread_flat)[0] == 0:
        print("The standard deviation matrix has no non-zero entries")
        sd_spread_flat = [0]
    sd_bins = int(np.max(sd_spread_flat) / 0.1) + 1

    avg_e = np.copy(aa_results["native_energy"])
    energy_flat = avg_e.flatten()
    energy_flat = energy_flat[nonzero_indices]

    if maxplot is None:
        maxplot = np.max(np.abs(zscores_flat))

    print "Plotting diagnostic plots"
    plot_obj.plot_spread(zscores_flat, savename="%s_spread" % save_name, xname="zscores")
    plot_obj.plot_spread(sd_spread_flat, savename="%s_sd" % save_name, xname="zscores", bins=sd_bins)
    if np.shape(energy_flat)[0] != 0:
        plot_obj.plot_spread(energy_flat, savename="%s_energy" % save_name, xname="zscores")

    gauss_file = "%s/decoy_gaussian_reduced_chi2.dat" % (main_dir)
    counts_file = "%s/decoy_gaussian_counts.dat" % (main_dir)
    if os.path.isfile(gauss_file):
        gauss = np.loadtxt(gauss_file)
        counts = np.loadtxt(counts_file)
        # destroy the counts file
        min_used_states = np.where(counts < min_use)
        gauss[min_used_states] = 0
        for i_gs in range(nresidues):
            for j_gs in range(nresidues):
                if np.abs(i_gs - j_gs) <= min_res_spacing:
                    gauss[i_gs, j_gs] = 0

        # now take the logarithm of non-zero points
        nonzero = np.where(gauss > 0)
        gauss[nonzero] = np.log10(gauss[nonzero])

    if dimensions == 2:
        # compute aa_plot_frustration: 2-D map of frustration values
        aa_plot_frustration = np.copy(aa_frustration)
        if reference_map is not None:
            reference_copy = np.copy(reference_map)
            if rescale_reference:
                max_ref = np.max(np.abs(reference_copy))
                awsem_scale_factor = maxplot / max_ref
                print "Adjusting AWSEM scale by factor %f" % (awsem_scale_factor)
                reference_copy = reference_copy * awsem_scale_factor
            for idx in range(nresidues):
                for jdx in range(idx+1, nresidues):
                    try:
                        #assert aa_frustration[idx,jdx] == aa_frustration[jdx,idx] # check symmetry before over-write
                        #TODO : make this assertion check pass
                        aa_plot_frustration[jdx, idx] = reference_copy[jdx, idx]
                    except:
                        print aa_plot_frustration[jdx, idx]
                        raise

        print "Plotting frustration for 2 dimensions"
        plot_obj.plot_epsilon_map([aa_plot_frustration], [edges], native_idxs, [""], save_name, maxval=maxplot, dca_pairs=dca_idxs, ncols=1, nrows=1, label_columns=False, grid_color="black", dca_color="black", top_half_only=False, color_map_bounds=[low_bound, high_bound], mark_values=mark_values)

        plot_obj.plot_epsilon_map([aa_plot_frustration], [edges], native_idxs, [""], "%s_full" % save_name, maxval=maxplot, dca_pairs=dca_idxs, ncols=1, nrows=1, label_columns=False, grid_color="black", dca_color="black", top_half_only=False, mark_values=mark_values)

        # save the gaussian fit data:
        print "Plotting Gaussian Chi results for 2 dimensions"
        if os.path.isfile(gauss_file):
            try:
                plot_obj.plot_epsilon_map([gauss], [edges], native_idxs, [""], "%s_chi2" % (save_name), maxval=maxplot, dca_pairs=dca_idxs, ncols=1, nrows=1, label_columns=False, grid_color="black", dca_color="black", top_half_only=False, use_color_map="viridis", color_map_bounds=[0, np.max(gauss)])
            except:
                pass

    elif dimensions == 1:
        print "Plotting single residue plots"
        aa_plot_frustration = np.copy(aa_frustration)
        res_indices = np.arange(np.shape(aa_plot_frustration)[0])
        plt.figure()
        plot_obj.plot_multilines([res_indices, res_indices], [reference_map, aa_plot_frustration], ["Reference", "AA-Frust"], save_name, mark_values=mark_values)
        if os.path.isfile(gauss_file):
            pass

    np.savetxt("%s/%s_aa_zscores.dat" % (plot_obj.figures_dir, save_name), aa_frustration)
    np.savetxt("%s/%s_aa_zscores.dat" % (plot_obj.png_figures_dir, save_name), aa_frustration)

    all_z_values = []
    all_ref_values = []
    for idx in range(nresidues):
        for jdx in range(idx+1, nresidues):
            if aa_plot_frustration[idx, jdx] != 0:
                all_z_values.append(aa_plot_frustration[idx,jdx])
            if aa_plot_frustration[jdx, idx] != 0:
                all_ref_values.append(aa_plot_frustration[jdx,idx])

    avg_zscore = np.sum(all_z_values) / np.shape(all_z_values)[0]
    avg_refscore = np.sum(all_ref_values) / np.shape(all_ref_values)[0]

    n_positive_z_values = 0
    zvalues_pos, zvalues_mid_pos, zvalues_neg, zvalues_mid_neg = calculate_z_distribution(all_z_values)

    z_total = len(all_z_values)
    if z_total == 0:
        z_total = 1

    rvalues_pos, rvalues_mid_pos, rvalues_neg, rvalues_mid_neg = calculate_z_distribution(all_ref_values)

    r_total = len(all_ref_values)
    if r_total == 0:
        r_total = 1

    if shift_value is None:
        save_shift_value = "0"
    else:
        save_shift_value = "%f" % shift_value

    info_str = ""
    info_str += "shift value = %s\n" % save_shift_value
    info_str += "average z-score / pair = %f\n" % avg_zscore
    info_str += "Fraction z-score > 1 : %f\n" % (float(zvalues_pos) / z_total)
    info_str += "Fraction z-score < -1 : %f\n" % (float(zvalues_neg) / z_total)
    info_str += "Fraction 0 < z-score <= 1 : %f\n" % (float(zvalues_mid_pos) / z_total)
    info_str += "Fraction -1 <= z-score < 0 : %f\n" % (float(zvalues_mid_neg) / z_total)

    info_str += "\n"
    info_str += "Reference info below: \n"
    info_str += "average reference score / pair = %f\n" % avg_refscore
    info_str += "Fraction z-score > 1 : %f\n" % (float(rvalues_pos) / z_total)
    info_str += "Fraction z-score < -1 : %f\n" % (float(rvalues_neg) / z_total)
    info_str += "Fraction 0 < z-score <= 1 : %f\n" % (float(rvalues_mid_pos) / z_total)
    info_str += "Fraction -1 <= z-score < 0 : %f\n" % (float(rvalues_mid_neg) / z_total)

    save_file = "%s/%s_info.txt" % (plot_obj.figures_dir, save_name)
    f = open(save_file, "w")
    f.write(info_str)
    f.write("\n")
    f.close()

    save_file = "%s/%s_info.txt" % (plot_obj.png_figures_dir, save_name)
    f = open(save_file, "w")
    f.write(info_str)
    f.write("\n")
    f.close()

    return shift_value, aa_frustration

def analyze_multi_singleresidues(main_dir_list, label_lists, plot_obj, save_name, maxplot=None, reference_map=None, min_value=0.1, low_bound=-5, high_bound=5, zero_bounds=[-1, 0.78], nresidues=35, apply_shift=False, all_shift_values=None, mark_values=None):

    all_zscores = [reference_map]
    for main_dir,shift_value in zip(main_dir_list, all_shift_values):
        print "Analyzing %s" % main_dir
        aa_results = load_frustration_results(main_dir, ndim=1)
        aa_frustration, resultant_shift_value = compute_frustration(aa_results["native_energy"], aa_results["decoy_avg"], aa_results["decoy_sd"], min_value=min_value, low_bound=low_bound, high_bound=high_bound, zero_bounds=zero_bounds, nresidues=nresidues, apply_shift=apply_shift, shift_value=shift_value)
        all_zscores.append(aa_frustration)
    print "Found %d zscores" % (np.shape(all_zscores)[0])
    print "Plotting single residue plots"
    aa_plot_frustration = np.copy(aa_frustration)
    res_indices = np.arange(np.shape(aa_plot_frustration)[0])
    all_x = [res_indices for i in all_zscores]
    plt.figure()
    plot_obj.plot_multilines(all_x, all_zscores, label_lists, save_name, mark_values=mark_values)

    all_avgs = []
    all_sq_avgs = []
    for zs in all_zscores:
        n_stuff = np.shape(zs)[0]
        avgs = np.sum(zs) / n_stuff
        neg_idx = np.where(zs < 0)
        switches = np.ones(n_stuff)
        switches[neg_idx] = -1
        sq_avgs = np.sum(zs * zs * switches)
        all_avgs.append(avgs)
        all_sq_avgs.append(sq_avgs)

    for idx in range(len(label_lists)):
        print "For %s, computed avg = %f and sq_avg = %f" % (label_lists[idx], all_avgs[idx], all_sq_avgs[idx])

# this run method will be run for every directory
def plot_final_plots(plot_obj, native_idxs, dca_idxs, bounds=5, nresidues=1, mark_values=None, reference_map=None, min_value=0.1):
    if reference_map is not None:
        if reference_map == "a3d":
            adir, astr = get_alpha3d_awsem_strings()
            c_matrix, m_matrix = load_awsem_frustration(adir, astr, nresidues=nresidues)
            oneindex, single_frust = load_awsem_single_residue(adir, astr, nresidues=nresidues)
        elif reference_map == "fip35":
            adir, astr = get_fip35_awsem_strings()
            c_matrix, m_matrix = load_awsem_frustration(adir, astr, nresidues=nresidues)
            oneindex, single_frust = load_awsem_single_residue(adir, astr, nresidues=nresidues)
    else:
        c_matrix = None
        m_matrix = None
        single_frust = None

    this_name = "mutational_justpairwise"
    shift_value, mutation_values = analyze_directory(this_name, plot_obj, native_idxs, dca_idxs, "mutational", reference_map=m_matrix, min_value=min_value, low_bound=-bounds, high_bound=bounds, zero_bounds=[-0, 0], nresidues=nresidues, apply_shift=False, shift_value=None, dimensions=2, rescale_reference=False, min_res_spacing=3)

    true_config_values = np.zeros(np.shape(config_values))
    for i in range(73):
        for j in range(i+1, 73):
            true_config_values[i,j] = config_values[i,j]
            true_config_values[j,i] = config_values[i,j]

    shift_value, frustration_values = analyze_directory(this_name, plot_obj, native_idxs, dca_idxs, "comparison", reference_map=true_config_values, min_value=0.001, low_bound=-3, high_bound=3, zero_bounds=[-0.001, 0.001], nresidues=nresidues, apply_shift=False, shift_value=None, dimensions=2, rescale_reference=False, min_res_spacing=3)

    this_name = "config_comparison_hetero1_restricted_decoy%d_stride100" % config_decoy_cutoff
    shift_value, frustration_values = analyze_directory(this_name, plot_obj, native_idxs, dca_idxs, "comparison_flipped", reference_map=mutation_values, min_value=0.001, low_bound=-3, high_bound=3, zero_bounds=[-0.001, 0.001], nresidues=nresidues, apply_shift=False, shift_value=None, dimensions=2, rescale_reference=False, min_res_spacing=3)
    return
    all_dirs = []
    all_dirs.append("singleresidue_mutational_frustration_results_removehigh100_decoys1000")
    all_dirs.append("config_comparison_singleresidue_restricted_decoy%d_stride100_ERemove100"% config_decoy_cutoff)
    analyze_multi_singleresidues(all_dirs, ["AWSEM", "Frust_mut", "Frust_cnfg"], plot_obj, "singleresidue_confg-mutational", reference_map=single_frust, min_value=0.001, low_bound=-3, high_bound=3, zero_bounds=[-0.001, 0.001], nresidues=nresidues, mark_values=mark_values, all_shift_values=[None, None], min_res_spacing=3)

    all_dirs = []
    all_dirs.append("singleresidue_mutational_frustration_results_removehigh100_decoys1000")
    all_dirs.append("singleresidue_mutational_muta3c_frustration_results_removehigh100_decoys1000")
    analyze_multi_singleresidues(all_dirs, ["AWSEM", "a3D", "a3C"], plot_obj, "singleresidue_mutational", reference_map=single_frust, min_value=0.001, low_bound=-5, high_bound=5, zero_bounds=[-0.001, 0.001], nresidues=nresidues, mark_values=mark_values, all_shift_values=[None, None], min_res_spacing=3)

    # do the single residue plotting

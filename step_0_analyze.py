import argparse
import pickle
import random
import time
import os
import math
import matplotlib.pyplot as plt
import pylab
from MotifSearch import *
from weblogolib import *
from corebio.seq import *


def debug(message):
    """Prints message with a timestamp of the time since the last debug() call."""
    # Put the before variable in the globals so time-tracking is scopeless
    if "before" not in globals():
        global before
        before = time.time()
    # Only print if args.debug flag is set (i.e. "-debug" in commandline args)
    if args.debug:
        print "[%.2fs] %s" % (time.time() - before, message)
        before = time.time()


def parse_args():
    """Parses commandline parameters and returns a ArgumentParser().parse_args() object.
    See: http://docs.python.org/2.7/library/argparse.html"""
    parser = argparse.ArgumentParser()
    parser.add_argument("pickled_file")
    parser.add_argument("-debug", action="store_true")
    parser.add_argument("-logos", action="store_true")
    parser.add_argument("-prob_mass", action="store_true")
    parser.add_argument("-kl", action="store_true")
    return parser.parse_args()


def load_pickle(filename):
    """Loads a pickled file. For reference see:
    http://docs.python.org/2/library/pickle.html"""
    with open("data/pickles/%s.pkl" % filename, "rb") as pickle_file:
        unpickled_data = pickle.load(pickle_file)
    return unpickled_data


def seqs_from_psfm(psfm, how_many):
    """Calculates a cumulative distribution function for each column of the PSFM and randomly samples
    the distribution to generate sequences according to the PSFM frequencies for each position."""
    sequences = []
    # Generate how_many sequences
    while len(sequences) < how_many:
        sequence = ""
        for column in psfm:
            # Create the cumulative distribution function for this column
            cdf = [(i, sum(p for j, p in column.items() if j < i)) for i, _ in column.items()]
            # Sample a base according to the CDF
            base = max(i for r in [random.random()] for i, c in cdf if c <= r)
            # Add it to the sequence
            sequence += base
        sequences.append(sequence)
    return sequences


def generate_logo(psfm, folder, filename, title):
    """Create a weblogo from a PSFM by generating sequences according to the PSFM frequencies."""
    # Generate sequences from PSFM
    sequences = seqs_from_psfm(psfm, 1000)
    # weblogo options and parameters
    logo_seqs = SeqList(sequences, unambiguous_dna_alphabet)
    logo_data = LogoData.from_seqs(logo_seqs)
    logo_options = LogoOptions(color_scheme=colorscheme.nucleotide, resolution=192, logo_title=title,
        title_fontsize=6)
    logo_format = LogoFormat(logo_data, logo_options)
    logos_folder = "figures/%s/logo" % folder
    # Create the folder if it does not exist
    if not os.path.exists(logos_folder):
        os.makedirs(logos_folder)
    # Write out to file
    with open("%s/%s.png" % (logos_folder, filename), "wb") as logo_file:
        png_formatter(logo_data, logo_format, logo_file)


def plot_prob_sums(prob_sums, folder):
    """Plot the probability sums of the original motif against iterations."""
    # Create the folder if it does not exist
    if not os.path.exists("figures/%s" % folder):
        os.makedirs("figures/%s" % folder)

    # Plot on automatic axis
    figure_file = "figures/%s/probability_mass.png" % folder
    plt.title("%s | Prob Mass of Original Motif" % folder)
    plt.ylabel("Probability Mass")
    plt.xlabel("Iterations")
    plt.grid(True)
    plt.plot(prob_sums)
    pylab.savefig(figure_file, format="png")
    plt.clf()
    
    # Plot on semi-log axis
    figure_file = "figures/%s/probability_mass_semilog.png" % folder
    plt.title("%s | Prob Mass of Original Motif" % folder)
    plt.xlabel("Iterations")
    plt.ylabel("Probability Mass (semi-log)")
    plt.grid(True)
    plt.semilogy(prob_sums)
    pylab.savefig(figure_file, format="png")
    plt.clf()


def plot_KLs(KLs, folder, filename, title):
    if not os.path.exists("figures/%s" % folder):
        os.makedirs("figures/%s" % folder)

    # Plot on automatic axis
    figure_file = "figures/%s/%s.png" % (folder, filename)
    plt.title("%s | %s" % (folder, title))
    plt.ylabel("KL")
    plt.xlabel("Iterations")
    plt.grid(True)
    plt.plot(KLs)
    pylab.savefig(figure_file, format="png")
    plt.clf()
    
    # Plot on semi-log axis
    try:
        figure_file = "figures/%s/%s_semilog.png" % (folder, filename)
        plt.title("%s | %s" % (folder, title))
        plt.xlabel("Iterations")
        plt.ylabel("KL (semi-log)")
        plt.grid(True)
        plt.semilogy(KLs)
        pylab.savefig(figure_file, format="png")
    except:
        debug("Could not generate semilog plot.")
    plt.clf()


if __name__ == '__main__':
    args = parse_args()
    debug("Parsed arguments.")
    pickled_file = args.pickled_file
    data = load_pickle(pickled_file)
    debug("Loaded pickled data.")

    if args.logos:
        debug("Generating logos.")
        for i, psfm in enumerate(data["PSFMs"]):
            folder = data["filename"][:-4]
            generate_logo(psfm, folder, "logo_i%d" % i, "%s | Iter %d" % (folder, i))
        debug("Finished generating logos.")

    if args.prob_mass:
        prob_sums = [sum(iter_probs) for iter_probs in data["probabilities"]]
        plot_prob_sums(prob_sums, data["filename"][:-4])
        debug("Generated probability mass figures.")

    if args.kl:
        KLs_from_original = data["KLs"][0]
        folder = data["filename"][:-4]
        plot_KLs(KLs_from_original, folder, "KLs_from_original", "KL Divergences from Original Motif")

        KLs_from_previous = [data["KLs"][i][i - 1] for i in range(len(data["KLs"])) if i is not 0]
        plot_KLs(KLs_from_previous, folder, "KLs_from_previous", "KL Divergences from Previous Iteration")
        debug("Generated KL figures.")
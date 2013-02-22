import argparse
import csv
import math
import os
import time
import sys
import matplotlib.pyplot as plt
import pylab
from MotifSearch import *
from weblogolib import *
from corebio.seq import *

# Debugging
before = time.time()
def debug(message):
    global before
    if args.debug:
        print "[%.2fs] %s" % (time.time() - before, message)
        before = time.time()

# Commandline argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("genome_file")
parser.add_argument("motif_file")
parser.add_argument("-data")  # Data file
parser.add_argument("-score")  # Scoring method
parser.add_argument("-be")  # Binding energy
parser.add_argument("-prob")  # Probability method
parser.add_argument("-weight")  # Weighing method
parser.add_argument("-temp")  # Temperature
parser.add_argument("-cn")  # Copy number
parser.add_argument("-iterations")  # Iteration cutoff
parser.add_argument("-ct")  # Convergence threshold
parser.add_argument("-top")  # Top K sites to ouput
parser.add_argument("-sum")  # Probability sum
parser.add_argument("-plot", action="store_true")  # Output plots to /figures
parser.add_argument("-logo", action="store_true")  # Generate WebLogos
parser.add_argument("-debug", action="store_true")  # Debug mode (verbose console logging)
args = parser.parse_args()
debug("Parsed arguments.")

# Initialize variables from arguments
genome = Genome("../TFBSDB/genomes/%s.fna" % args.genome_file, FORMAT_FASTA)
original_motif = SequenceCollection("%s.txt" % args.motif_file, FORMAT_PLAIN)
iteration_cutoff = int(args.iterations) if args.iterations else 100  # Default cutoff is at 100 iterations
convergence_threshold = args.ct if args.ct else 0.001  # Default KL change threshold is 0.05
data_file = args.data + ".csv" if args.data else "data.csv"  # Adds the .csv extension, defaults to outputting to data.csv if no filename specified
scoring_method = str(args.score)
weighing_method = str(args.weight)
probability_method = str(args.prob)
debug("Initialized variables.")

# Write metadata to file
with open("data/" + data_file, "wb") as csv_file:
    writer = csv.writer(csv_file)
    writer.writerows([
        ["Genome:", genome.name],
        ["TFBS:", args.motif_file],
        ["Scoring:", scoring_method, "Weighing:", weighing_method, "Probability:", probability_method],
        ["Command:", " ".join(sys.argv)]
        ])

# Probability parameters
if args.prob:
    temperature = int(args.temp) if args.temp else 310  # Default is 310 Kelvin
    copy_number = int(args.cn) if args.cn else original_motif.size  # Default is size of the motif
    beta = 1 / (kB * temperature)

# Initialize search variables
iteration = 0
converged = False
frequencies = []
frequencies_wc = []
divergences = []

# Create matrices for the original motif
motif = original_motif  # Note: This is a reference! Use .deepcopy() to create a new instance
debug("Set motif to original_motif.")

count_matrix = PSCM(motif, STRAND_ORIGINAL)
count_matrix_wc = wc_matrix(count_matrix)
debug("Calculated count matrices for original motif.")

# Calculate the frequency matrices (PSFM)
frequency_matrix = PSFM(count_matrix)
frequency_matrix_wc = PSFM(count_matrix_wc)
debug("Calculated frequency matrices for original motif.")

# Calculate the scoring matrices (both strands)
scoring_matrix = create_scoring_matrix(scoring_method, genome, count_matrix, frequency_matrix)
scoring_matrix_wc = create_scoring_matrix(scoring_method, genome, count_matrix_wc, frequency_matrix_wc)
debug("Calculated scoring matrices for original motif.")

frequencies.append(frequency_matrix)
frequencies_wc.append(frequency_matrix_wc)
original_count_matrix = count_matrix
original_count_matrix_wc = count_matrix_wc

debug("Started iterating.")
while not converged and iteration < iteration_cutoff:
    debug("Iteration %d\n=====================" % iteration)

    # Create sliding window generator
    window_generator = sliding_window(genome, original_motif.width)

    # Weigh the scoring matrices
    if args.weight:
        # Information Content
        if weighing_method in "ic":
            scoring_matrix = [{base: (score * motif.ic[i]) for base, score in position.items()} for i, position in enumerate(scoring_matrix)]
            scoring_matrix_wc = [{base: (score * motif.ic[i]) for base, score in position.items()} for i, position in enumerate(scoring_matrix_wc)]
        # Sigmoidal function weighing
        elif weighing_method in "sigma":
            scoring_matrix = [{base: (sigma(motif.ic[i], 1, 5) * score + (1 - sigma(motif.ic[i], 1, 5)) * ((score + position[wc(base)]) / 2) )
                for base, score in position.items()} for i, position in enumerate(scoring_matrix)]
            scoring_matrix_wc = [{base: (sigma(motif.ic[i], 1, 5) * score + (1 - sigma(motif.ic[i], 1, 5)) * ((score + position[wc(base)]) / 2) )
                for base, score in position.items()} for i, position in enumerate(scoring_matrix_wc)]
        debug("Weighed the scoring matrices.")

    # Iteration 0 case; Score the original motif
    if iteration is 0:
        for i, site in original_motif.collection_generator():
            score, strand = score_sequence(site.sequence, scoring_matrix, scoring_matrix_wc)
            original_motif.update(i, {"score": score, "strand": strand})  # Update only the strand and score
        original_motif.statistics()  # Calculate statistics

    # Create a new SequenceCollection object for a motif found by using the previous motif ~~ WHERE THE MAGIC HAPPENS ~~
    motif = SequenceCollection()
    # Score every sequence in the genome and save to motif
    for position, sequence in window_generator:
        score, strand = score_sequence(sequence, scoring_matrix, scoring_matrix_wc)  # Score forward and reverse
        motif.add({"position": position, "sequence": sequence if strand is STRAND_ORIGINAL else wc(sequence), "strand": strand, "score": score})
    debug("Scored the genome and built collection (%d bytes)." % sys.getsizeof(motif.collection))

    # If we're doing probability, scale the score to binding energy
    # Note: This assumes we're using Berg and von Hippel scoring
    if args.prob or args.be:
        total_energy = 0
        if iteration is 0:
            for i, site in enumerate(original_motif):
                energy = linear_scale(site.score, -0.05, -2 * kB * temperature * original_motif.width, original_motif.stats["score"]["min"])
                original_motif.update(i, {"energy": energy})
        for i, site in motif.collection_generator():
            # Calculate the energy of binding by scaling the BvH score to [-0.05, -2kBTl]
            energy = linear_scale(site.score, -0.05, -2 * kB * temperature * original_motif.width, original_motif.stats["score"]["min"])
            motif.update(i, {"energy": energy})
            total_energy += math.exp(-beta * energy)  # Add to the total sum of energies partition function (Z)
        debug("Scaled scores to energies.")

    # Calculate probabilities for the entire motif
    if args.prob:
        debug("Beta: " + str(beta))
        debug("Z: " + str(total_energy))
        # Calculate probabilities for original motif
        if iteration is 0:
            for i, site in enumerate(original_motif):
                # Boltzmann probabilities (single TF molecule)
                if probability_method in "boltzmann":
                    original_motif.update(i, {"probability": boltzmann_probability(site.energy, total_energy, beta)})
                # Fermi-Dirac probabilities (n TF molecules; n = copy number)
                if probability_method in "fd":
                    original_motif.update(i, {"probability": fermi_dirac_probability(site.energy, total_energy, beta, copy_number)})
            original_motif.normalize("probability")
            original_motif.sort("probability")
        for i, site in motif.collection_generator():
            # Boltzmann probabilities (single TF molecule)
            if probability_method in "boltzmann":
                motif.update(i, {"probability": boltzmann_probability(site.energy, total_energy, beta)})
            # Fermi-Dirac probabilities (n TF molecules; n = copy number)
            if probability_method in "fd":
                motif.update(i, {"probability": fermi_dirac_probability(site.energy, total_energy, beta, copy_number)})
            # Normalize
        motif.normalize("probability")
        motif.sort("probability")
        debug("Calculated binding probabilities.")


    # Truncate the new motif
    if args.sum and args.prob:
        prob_sum_cutoff = float(args.sum)
        motif.truncate(prob_sum_cutoff)  # Sorts the motif by probability and truncates to only the top sites with certain probability mass

    #### FOR THE NEXT ITERATION ####
    # Calculate the count matrices (PSCM)
    count_matrix = PSCM(motif, STRAND_ORIGINAL)
    count_matrix_wc = wc_matrix(count_matrix)
    debug("Calculated count matrices.")

    # Calculate the frequency matrices (PSFM)
    frequency_matrix = PSFM(count_matrix)
    frequency_matrix_wc = PSFM(count_matrix_wc)
    debug("Calculated frequency matrices.")

    # Calculate the scoring matrices (both strands)
    scoring_matrix = create_scoring_matrix(scoring_method, genome, count_matrix, frequency_matrix)
    scoring_matrix_wc = create_scoring_matrix(scoring_method, genome, count_matrix_wc, frequency_matrix_wc)
    debug("Calculated scoring matrices.")

    frequencies.append(frequency_matrix)
    frequencies_wc.append(frequency_matrix_wc)

    # Create SequenceCollection object with new generalized motif
    motif.build_collection(motif.collection)
    debug("Initialized SequenceCollection object to motif.")

    # Calculate entropy and IC for the motif
    if iteration is 0:
        original_motif.entropies = Hxy(frequencies[0])
        original_motif.ic = RI(genome.entropy, original_motif.entropies)
    motif.entropies = Hxy(frequency_matrix)
    motif.ic = RI(genome.entropy, motif.entropies)

    # Calculate Kullback-Leibler divergences
    kl_from_original = sum([kl(p.values(), q.values()) for p, q in zip(frequency_matrix, frequencies[0])])
    kl_from_last = sum([kl(p.values(), q.values()) for p, q in zip(frequency_matrix, frequencies[-2])])
    divergences.append(kl_from_last)
    delta_kl = abs(kl_from_last - divergences[-1])
    debug("Calculated KLs. From last: %.2f, original: %.2f" % (kl_from_last, kl_from_original))

    # Output to file
    with open("data/" + data_file, "ab") as csv_file:
        writer = csv.writer(csv_file)
        # Original motif
        if iteration is 0:
            writer.writerows(
                [["Iteration:", "0"],
                ["PSCM:"],
                [""] + [n for n in range(original_motif.width)],
                ["A"] + [position["A"] for position in original_count_matrix],
                ["T"] + [position["T"] for position in original_count_matrix],
                ["C"] + [position["C"] for position in original_count_matrix],
                ["G"] + [position["G"] for position in original_count_matrix],
                ["PSFM:"],
                [""] + [n for n in range(original_motif.width)],
                ["A"] + [position["A"] for position in frequencies[0]],
                ["T"] + [position["T"] for position in frequencies[0]],
                ["C"] + [position["C"] for position in frequencies[0]],
                ["G"] + [position["G"] for position in frequencies[0]],
                ["Consensus:", consensus(original_count_matrix)],
                ["KL(0, 0):", "N/A", "KL (0, 0)", "N/A", "Delta:", "N/A"]])
            debug("Wrote PSCM and PSFM to file. Consensus: %s" % consensus(count_matrix))

            # Generate statistics from the original_motif.collection
            stats = original_motif.statistics()
            debug("Generated statistics for the original_motif.")

            # Score statistics
            writer.writerows([["Scores:"], ["Mean:", stats["score"]["mean"], "Min:", stats["score"]["min"], "Max:", stats["score"]["max"]]])
            debug("Wrote scores to file.")

            # Binding Energy statistics
            if args.be or args.prob:
                writer.writerows([["Binding Energy:"], ["Mean:", stats["energy"]["mean"], "Min:", stats["energy"]["min"], "Max:", stats["energy"]["max"], "Total:", stats["energy"]["total"]]])
                debug("Wrote energies to file.")

            # Probability statistics
            if args.prob:
                writer.writerows([["Probability:"], ["Beta:", beta, "Z:", total_energy], ["Mean:", stats["probability"]["mean"], "Min:", stats["probability"]["min"], "Max:", stats["probability"]["max"]]])
                debug("Wrote probabilities to file.")

        # PSCM and KL
        writer.writerows(
            [["Iteration:", iteration + 1],
            ["PSCM:"],
            [""] + [n for n in range(motif.width)],
            ["A"] + [position["A"] for position in count_matrix],
            ["T"] + [position["T"] for position in count_matrix],
            ["C"] + [position["C"] for position in count_matrix],
            ["G"] + [position["G"] for position in count_matrix],
            ["PSFM:"],
            [""] + [n for n in range(motif.width)],
            ["A"] + [position["A"] for position in frequency_matrix],
            ["T"] + [position["T"] for position in frequency_matrix],
            ["C"] + [position["C"] for position in frequency_matrix],
            ["G"] + [position["G"] for position in frequency_matrix],
            ["Consensus:", consensus(count_matrix)],
            ["KL(%d, 0):" % (iteration + 1), kl_from_original, "KL (%d, %d)" % (iteration + 1, iteration), kl_from_last, "Delta:", delta_kl]])
        debug("Wrote PSCM and PSFM to file. Consensus: %s" % consensus(count_matrix))

        # Generate statistics from the motif.collection
        stats = motif.statistics()
        debug("Generated statistics for the motif.")

        # Score statistics
        writer.writerows([["Scores:"], ["Mean:", stats["score"]["mean"], "Min:", stats["score"]["min"], "Max:", stats["score"]["max"]]])
        debug("Wrote scores to file.")

        # Binding Energy statistics
        if args.be or args.prob:
            writer.writerows([["Binding Energy:"], ["Mean:", stats["energy"]["mean"], "Min:", stats["energy"]["min"], "Max:", stats["energy"]["max"], "Total:", stats["energy"]["total"]]])
            debug("Wrote energies to file.")

        # Probability statistics
        if args.prob:
            writer.writerows([["Probability:"], ["Beta:", beta, "Z:", total_energy], ["Mean:", stats["probability"]["mean"], "Min:", stats["probability"]["min"], "Max:", stats["probability"]["max"]]])
            debug("Wrote probabilities to file.")

        # Top sites
        if args.top:
            top_sites = motif.top(int(args.top), "probability")
            writer.writerows([["Top Sites:"], [field for field in top_sites[0]._fields]])  # Headers
            writer.writerows([[value for value in site] for site in top_sites])  # Sites
            debug("Wrote top sites to file.")
        
        # Top sites (by probability sum)
        if args.sum:
            writer.writerows([["Top Sites", "Probability Sum:", args.sum], [field for field in motif[0]._fields]])  # Headers
            writer.writerows([[value for value in site] for site in motif])  # Sites
            debug("Wrote top sites (sum) to file.")

    # Plot figures
    if args.plot:
        folder = "figures/" + data_file[:-4]
        if not os.path.exists(folder):
            os.makedirs(folder)
        if iteration is 0:
            # Information content
            figure_file = folder + "/ic_0.png"
            plt.title("Information Content - Iteration 0")
            plt.xlabel("Position")
            plt.ylabel("Bits")
            plt.ylim(0, 2.0)
            plt.xticks(range(len(original_motif.ic)))
            plt.bar(range(len(original_motif.ic)), original_motif.ic)
            pylab.savefig(figure_file, format="png")
            plt.clf()
            
            # Probability
            figure_file = folder + "/probs_0.png"
            plt.title("Probability - Iteration 0")
            plt.xlabel("Site")
            plt.ylabel("Probability")
            plt.ylim(0, 1.0)
            plt.plot([site.probability for site in original_motif], "bo")
            pylab.savefig(figure_file, format="png")
            plt.clf()
        # Information content
        figure_file = folder + "/ic_" + str(iteration + 1) + ".png"
        plt.title("Information Content - Iteration " + str(iteration + 1))
        plt.xlabel("Position")
        plt.ylabel("Bits")
        plt.ylim(0, 2.0)
        plt.xticks(range(len(motif.ic)))
        plt.bar(range(len(motif.ic)), motif.ic)
        pylab.savefig(figure_file, format="png")
        plt.clf()
        
        # Probability

        figure_file = folder + "/probs_%d.png" % (iteration + 1)
        plt.title("Probability - Iteration %d" % (iteration + 1))
        plt.xlabel("Site")
        plt.ylabel("Probability")
        plt.ylim(0, 1.0)
        plt.plot([site.probability for site in motif], "bo")
        pylab.savefig(figure_file, format="png")
        plt.clf()
        debug("Generated figures.")

    if args.logo:
        if iteration is 0:
            logo_file = "logo_i0"
            logo_seqs = SeqList([site.sequence for site in original_motif], unambiguous_dna_alphabet)
            logo_data = LogoData.from_seqs(logo_seqs)
            logo_options = LogoOptions(color_scheme=colorscheme.nucleotide, resolution=192, logo_title="Iteration 0",
                title_fontsize=6)
            logo_format = LogoFormat(logo_data, logo_options)
            folder = "logos/" + data_file[:-4]
            if not os.path.exists(folder):
                os.makedirs(folder)
            with open("%s/%d.png" % (folder, i), "wb") as logo_file:
                png_formatter(logo_data, logo_format, logo_file)
            debug("Generated logo.")
        # Generate logo
        logo_file = "logo_i%d" % (iteration + 1)
        logo_seqs = SeqList([site.sequence for site in motif], unambiguous_dna_alphabet)
        logo_data = LogoData.from_seqs(logo_seqs)
        logo_options = LogoOptions(color_scheme=colorscheme.nucleotide, resolution=192, logo_title="Iteration %d" % (iteration + 1),
            title_fontsize=6)
        logo_format = LogoFormat(logo_data, logo_options)
        folder = "logos/" + data_file[:-4]
        if not os.path.exists(folder):
            os.makedirs(folder)
        with open("%s/%d.png" % (folder, i), "wb") as logo_file:
            png_formatter(logo_data, logo_format, logo_file)
        debug("Generated logo.")

    # Convergence check
    if iteration is not 0:
        if delta_kl < convergence_threshold and args.ct:
            converged = True
            debug("Convergence detected. dKL = %f" % delta_kl)

    # Iterate
    iteration += 1
    debug("Iteration finished.\n")

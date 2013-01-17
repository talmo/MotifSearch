import csv
import time
import sys
from bisect import bisect
from MotifSearch import *


genome = Genome("../TFBSDB/genomes/NC_000913.fna", FORMAT_FASTA)  # E. coli
#motif = SequenceCollection("crp.prodoric.txt", format=FORMAT_PLAIN)  # CRP
motif = SequenceCollection("lexa.prodoric.txt", format=FORMAT_PLAIN)  # LexA
data_file = "data/data_BvH.csv" if len(sys.argv) < 2 else sys.argv[1]
scoring_method = SCORING_BVH
temperature = 310  # Kelvin
iteration_cutoff = 100

print "Genome:", genome.name
print "Motif: %d sites, %dbp width" % (motif.collection_size, motif.width)
print "Scoring:", SCORING_STRINGS[scoring_method]
print "Data File:", data_file

with open(data_file, "wb") as csv_file:
    writer = csv.writer(csv_file)
    writer.writerows([["Genome:", genome.name], ["Motif (#, width):", len(motif.collection), motif.collection_size, motif.width], ["Scoring:", SCORING_STRINGS[scoring_method]], ["Iteration", "Sequence", "Position", "Probability"]])

beta = 1 / (Kb * temperature)
# Initial matrices (updated every iteration)
count_matrix = PSCM(motif.sequences, motif.width)
scoring_matrix = berg_von_hippel_matrix(count_matrix)

# This is needed to calculate the linear scaling parameters
motif_scores = [scorer(sequence, method=SCORING_BVH, matrix=scoring_matrix) for sequence in motif.sequences]

# Simulation variables
converged = False
iteration = 0
window_generator = sliding_window(genome.sequence, motif.width)
all_n_mers = [(window[0], window[1]) for window in window_generator]  # Store all n-mers in the genome in a list of tuples (position, sequence)
while not converged:
    b = time.time()
    # Calculate Berg and von Hippel scores
    genome_scores = [scorer(sequence, method=SCORING_BVH, matrix=scoring_matrix) for i, sequence in all_n_mers]

    # Scale scores to binding energies
    scaled_genome_scores = linear_scale(genome_scores, -0.05, -2 * Kb * temperature * motif.width, min(motif_scores))

    # Partition function (Z)
    total_genome_energy = sum([math.e ** (-beta * energy) for energy in scaled_genome_scores])

    # Calculate binding probabilities for each site
    genome_probabilities = [binding_probability(sequence, scaled_genome_scores[i], total_genome_energy, beta) for i, sequence in all_n_mers]

    # Write out the original motif's scores to a CSV file
    if iteration == 0:
        scaled_motif_scores = linear_scale(motif_scores, -0.05, -2 * Kb * temperature * motif.width, min(motif_scores))
        motif_probabilities = [binding_probability(sequence, scaled_motif_scores[i], total_genome_energy, beta) for i, sequence in enumerate(motif.sequences)]
        print "Consensus: %s | Initial motif. | Probability Sum: %.4f" % (consensus(count_matrix), sum(motif_probabilities))
        with open(data_file, "ab") as csv_file:
            writer = csv.writer(csv_file)
            writer.writerows([[iteration, sequence, None, probability] for probability, sequence in zip(motif_probabilities, motif.sequences)])

    # Iterate by re-calculating the PSCM using probabilities
    count_matrix = weighted_PSCM(all_n_mers, motif.width, genome_probabilities)
    scoring_matrix = berg_von_hippel_matrix(count_matrix)

    # Get the top k sites
    top_sites = [None] * motif.collection_size
    top_sites_probs = [None] * motif.collection_size
    for (position, site), probability in zip(all_n_mers, genome_probabilities):
        if probability > min(top_sites_probs):
            insert_pos = bisect(top_sites_probs, probability)
            top_sites.insert(insert_pos, [position, site])
            top_sites_probs.insert(insert_pos, probability)
            top_sites.pop(0)
            top_sites_probs.pop(0)

    # Write out the top k sites to a CSV file
    with open(data_file, "ab") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerows([[iteration + 1, sequence, position, probability] for (position, sequence), probability in zip(reversed(top_sites), reversed(top_sites_probs))])

    # Output status of run to console
    print "Iteration: %3d | Consensus: %s | Probability Sum: %.4f | Done in %.2f secs." % (iteration + 1, consensus(count_matrix), sum(top_sites_probs), time.time() - b)
    if iteration == iteration_cutoff:
        converged = True
    iteration += 1

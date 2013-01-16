import time
from MotifSearch import *


genome = Genome("../TFBSDB/genomes/NC_000913.fna", FORMAT_FASTA)  # E. coli
motif = SequenceCollection("crp.prodoric.txt", format=FORMAT_PLAIN)  # crp
data_file = "data/k12_crp_BvH.csv"
scoring_method = SCORING_BVH
temperature = 310  # Kelvin
iteration_cutoff = 100

print "Genome:", genome.name
print "Motif: %d sites, %dbp width" % (motif.collection_size, motif.width)
print "Scoring:", SCORING_STRINGS[scoring_method]
print "Data File:", data_file

count_matrix = PSCM(motif.sequences, motif.width)
scoring_matrix = berg_von_hippel_matrix(count_matrix)

motif_scores = [scorer(sequence, method=SCORING_BVH, matrix=scoring_matrix) for sequence in motif.sequences]
scaled_motif_scores = linear_scale(motif_scores, -0.05, -2 * Kb * temperature * motif.width, min(motif_scores))
print "Motif Scores: best = %.2f, worst = %.2f" % (max(motif_scores), min(motif_scores))
print "Motif Scaled: best = %.2f, worst = %.2f" % (min(scaled_motif_scores), max(scaled_motif_scores))
beta = 1 / (Kb * temperature)
print "Beta: %.2f" % beta

window_generator = sliding_window(genome.sequence, motif.width)
b = time.time()
genome_scores = [scorer(window[1], method=SCORING_BVH, matrix=scoring_matrix) for window in window_generator]
print "Genome Scores: best = %.2f, worst = %.2f (%.2f sec)" % (max(genome_scores), min(genome_scores), time.time() - b)

b = time.time()
scaled_genome_scores = linear_scale(genome_scores, -0.05, -2 * Kb * temperature * motif.width, min(motif_scores))
print "Genome Scaled: best = %.2f, worst = %.2f (%.2f sec)" % (min(scaled_genome_scores), max(scaled_genome_scores), time.time() - b)

b = time.time()
total_genome_energy = sum([math.e ** (-beta * energy) for energy in scaled_genome_scores])
print "Total genome energy: %d (%.2f sec)" % (total_genome_energy, time.time() - b)

b = time.time()
motif_probabilities = [binding_probability(sequence, scaled_motif_scores[i], total_genome_energy, beta) for i, sequence in enumerate(motif.sequences)]
window_generator = sliding_window(genome.sequence, motif.width)
genome_probabilities = [binding_probability(sequence, scaled_genome_scores[i], total_genome_energy, beta) for i, sequence in window_generator]
print "Motif Probabilities: best = %f, worst = %f (%.2f sec)" % (max(motif_probabilities), min(motif_probabilities), time.time() - b)
print "Genome Probabilities: best = %f, worst = %f (%.2f sec)" % (max(genome_probabilities), min(genome_probabilities), time.time() - b)

import csv
import time
import sys
from MotifSearch import *

genome = Genome("../TFBSDB/genomes/NC_000913.fna", FORMAT_FASTA)  # E. coli
motif = SequenceCollection("crp.prodoric.txt", format=FORMAT_PLAIN)  # crp
data_file = "data/data.csv" if len(sys.argv) < 2 else sys.argv[1]
scoring_method = SCORING_INFORMATION_SCORE
iteration_cutoff = 100

# Drop a third of the motif sequences
#motif.collection = motif.collection[:int(motif.collection_size * 0.66)]
#motif.sequences = motif.sequences[:int(motif.collection_size * 0.66)]
#motif.wc_sequences = motif.wc_sequences[:int(motif.collection_size * 0.66)]

print "Genome:", genome.name
print "Motif: %d/%d sites, %dbp width" % (len(motif.collection), motif.collection_size, motif.width)
print "Scoring:", SCORING_STRINGS[scoring_method]
print "Data File:", data_file
print "Iterating until convergence or %d iterations..." % iteration_cutoff

with open(data_file, "wb") as csv_file:
    writer = csv.writer(csv_file)
    writer.writerows([["Genome:", genome.name], ["Motif (#, width):", len(motif.collection), motif.collection_size, motif.width], ["Scoring:", SCORING_STRINGS[scoring_method]], ["Iteration", "Sequence", "Position", "Score"]])

converged = False
iteration = 0
while not converged:
    b = time.time()
    # Frequency matrices
    frequency_matrix = PSFM(motif.sequences, motif.width)
    wc_frequency_matrix = PSFM(motif.wc_sequences, motif.width)

    # Scoring matrices
    if scoring_method is SCORING_INFORMATION_SCORE:
        scoring_matrix = information_score_matrix(frequency_matrix, genome.entropy)
        wc_scoring_matrix = information_score_matrix(wc_frequency_matrix, genome.entropy)
    elif scoring_method is SCORING_WEIGHTED_INFORMATION_SCORE:
        scoring_matrix = weighted_information_score_matrix(frequency_matrix, genome.entropy)
        wc_scoring_matrix = weighted_information_score_matrix(wc_frequency_matrix, genome.entropy)

    # Score the initial motif in case it isn't already; this is not necessary for scoring but nice for visualization in the analysis
    if not motif.scored:
        for i, item in enumerate(motif.collection):
            score = scorer(item["sequence"], scoring_method, matrix=scoring_matrix)
            wc_score = scorer(item["sequence"], scoring_method, matrix=wc_scoring_matrix)
            motif.collection[i]["score"] = score if score > wc_score else wc_score
            motif.scored = True

    # Initialize iteration variables
    window_generator = sliding_window(genome.sequence, motif.width)
    genome_scores = []
    strandedness = []

    # Calculate the score for every position in the genome
    for position, sequence in window_generator:
        score = scorer(sequence, scoring_method, matrix=scoring_matrix)
        wc_score = scorer(sequence, scoring_method, matrix=wc_scoring_matrix)
        genome_scores.append(score if score > wc_score else wc_score)
        strandedness.append(STRAND_ORIGINAL if score > wc_score else STRAND_COMPLEMENT)  # Stores which strand the score corresponds to

    # Sort the scores
    sorted_genome_scores = sorted([[position, score] for position, score in enumerate(genome_scores)], key=lambda site: site[1], reverse=True)
    # Get the top k scores where k = the size of the original motif
    top_k = sorted_genome_scores[:len(motif.collection)]
    # Check for convergence
    if [site[0] for site in top_k] in [item["position"] for item in motif.collection] or iteration >= iteration_cutoff:
        converged = True
    # Write data to file
    with open(data_file, "ab") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerows([[iteration, site["sequence"], site["position"], site["score"]] for site in motif.collection])
    # Logging
    print "Iteration %d done in %.2f secs." % (iteration + 1, time.time() - b)
    # Update the motif and iterate
    motif = SequenceCollection([genome.sequence[site[0]:site[0] + motif.width] if strandedness[site[0]] is STRAND_ORIGINAL else wc(genome.sequence[site[0]:site[0] + motif.width]) for site in top_k], position=[site[0] for site in top_k], score=[site[1] for site in top_k], strand=[strandedness[site[0]] for site in top_k])
    iteration += 1

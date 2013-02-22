from MotifSearch import *
import argparse
import time

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
    parser.add_argument("genome")
    parser.add_argument("motif")
    parser.add_argument("-i")
    parser.add_argument("-score")
    parser.add_argument("-cn")
    parser.add_argument("-debug", action="store_true")
    parser.add_argument("-pickle", action="store_true")
    return parser.parse_args()


def run_iterations(genome, original_motif, **kwargs):
    """Bootstrap the iterate() function which does most of the heavy lifting as well as handle
    things between iterations.
        * genome must be a Genome instance (defined in MotifSearch),
        * original_motif must be a SequenceCollection instance (defined in MotifSearch).
    Iterating will stop when kwargs.how_many iterations have occurred or, if kwargs.until_converges
    is True, when convergence is detected, whichever comes first. If you want to iterate until
    convergence, set kwargs.how_many to -1.
    Returns a list of the PSCM and Z for every iteration, not including the Z for the last iteration."""

    # Figure out when to stop iterating and how to score
    until_converges = kwargs["until_converges"] if "until_converges" in kwargs else False
    how_many = kwargs["how_many"] if "how_many" in kwargs else 1
    scoring_method = kwargs["scoring_method"] if "scoring_method" in kwargs else "pssm"
    copy_number = kwargs["copy_number"] if "copy_number" in kwargs else original_motif.size
    debug("Copy Number: %d" % copy_number)

    # The very first PSCM! Every position is essentially weighted by a "probability" of 1.0
    initial_pscm = PSCM(original_motif)
    initial_scoring_matrix = create_scoring_matrix(scoring_method, genome, initial_pscm, PSFM(initial_pscm))
    debug("Calculated the initial matrices using the original motif.")

    # Iterate while we haven't met our stopping criteria
    iterations = 0
    converged = False
    PSCMs = [initial_pscm]  # The initial PSCM corresponds with the first Z
    Zs = []
    genome_energies = []
    genome_prob_masses = []
    debug("Beginning to iterate.\n")
    while (iterations < how_many or how_many is -1) and (not converged or until_converges is False):
        debug("Iteration %d\n====================" % iterations)
        # Iterate with initial matrices if this is the first iteration
        if iterations == 0:
            pscm, Z, iteration_prob_mass, iteration_energies = iterate(genome, initial_pscm, initial_scoring_matrix, copy_number)
        else:
            pscm, Z, iteration_prob_mass, iteration_energies = iterate(genome, pscm, scoring_matrix, copy_number)

        # Add the PSCM, Z, and mean score to the lists that will be returned
        PSCMs.append(pscm)
        Zs.append(Z)
        #genome_energies.append(iteration_energies)
        genome_prob_masses.append(iteration_prob_mass)

        # Re-calculate scoring matrix for the next iteration
        scoring_matrix = create_scoring_matrix(scoring_method, genome, pscm, PSFM(pscm))

        # Should we stop iterating?
        converged = check_for_convergence()
        iterations += 1
        debug("Finished iteration %d.\n" % (iterations - 1))
    debug("Finished iterating.")

    return PSCMs, Zs, genome_prob_masses

g_scores = []
g_energies = []
def iterate(genome, pscm, scoring_matrix, copy_number):
    """ ~*~ Where the magic happens! ~*~
    Returns the partition function Z and the new PSCM weighted by probability.
    The mean score is also returned as it is used as a scaling parameter for
    conversion from score to energy."""
    width = len(pscm)  # How many bases long the motif is
    temperature = 310  # Kelvin
    beta = 1 / (kB * temperature)

    # 1. Walk through the genome scoring each site and scaling to energy
    scores, genome_energies, Z = get_scores_energies((sequence for pos, sequence in sliding_window(genome.sequence, width)),
        scoring_matrix, beta, temperature, width)
    mean_score = sum(scores) / len(scores)
    mean_energy = sum(genome_energies) / len(genome_energies)
    debug("Calculated scores. Mean: %g, Min: %g, Max: %g" % (mean_score, min(scores), max(scores)))
    debug("Energies Mean: %g, Min: %g, Max: %g, Z: %g" % (mean_energy, min(genome_energies), max(genome_energies), Z))
    
    # For debugging
    #g_scores.append(scores)
    #g_energies.append(genome_energies)

    window_generator = sliding_window(genome.sequence, width)
    total_mass = 0
    genome_probabilities = []
    pscm = [{base: 0 for base in "ACTG"} for i in range(width)]

    # 2. Walk through the genome again, this time to calculate probabilities and the new PSCM
    for pos, sequence in window_generator:
        # Probability given the energy
        p = fermi_dirac_probability(genome_energies[pos], Z, beta, copy_number)

        # Weigh the counts by the probability to create a new PSCM
        for i, base in enumerate(sequence):
            pscm[i][base] += p
        genome_probabilities.append(p)

    #assert abs(1 - sum(genome_probabilities) / copy_number) < 0.05
    genome_prob_mass = sum(genome_probabilities)
    debug("Calculated probabilities and new PSCM. Genome Probability Mass: %g" % genome_prob_mass)

    return pscm, Z, genome_energies, genome_prob_mass


def get_scores_energies(sequences, scoring_matrix, beta, temperature, width):
    """Return the scores and energies for a given list sequences as scored against a scoring matrix."""
    method = SCORING_INFORMATION_SCORE
    scores = []
    energies = []
    Z = 0
    scoring_matrix_wc = wc_matrix(scoring_matrix)

    for sequence in sequences:
        score, strand = score_sequence(sequence, scoring_matrix, scoring_matrix_wc)
        scores.append(score)
        #energy = energy_scale(method, score, temperature=temperature, genomic_entropy=genome.entropy, width=width)
        energy = TRAP(score, width)
        energies.append(energy)
        Z += math.exp(-beta * energy)

    return scores, energies, Z


def get_energies(scores, beta, mean_score, width):
    """Returns a list of energy for every score in scores and the partition function Z."""
    energies = []  # Contains the energy for every site in the genome
    Z = 0  # This is the partition function, it's basically the sum of all math.exp(-beta * energy)

    # Calculate energies for every score
    for score in scores:
        energy = linear_scale(score, -0.05, -2 * (1 / beta) * width, mean_score)
        energies.append(energy)
        # Add the e ^ (-beta * energy) to the partition function
        Z += math.exp(-beta * energy)

    return energies, Z


def check_for_convergence():
    """TODO: Change in KL?"""
    return False


def PSCM_probabilities(genome, motif, pscm, genome_z, scoring_method, genome_mean_score, copy_number):
    """Returns a list of the probability for each site in the motif.
    Very important to note is that the partition function (Z) used is the one inputed, 
    not the one returned by get_energies(). This is because we need the Z for the ENTIRE GENOME.
    This is a value returned in a iterate() call."""
    # Calculate the scores
    scoring_matrix = create_scoring_matrix(scoring_method, genome, pscm, PSFM(pscm))
    scores = get_scores([site.sequence for site in motif], scoring_matrix)

    # Calculate energy for every site in the motif, but disregard the Z returned -- we need the genomic partition function
    temperature = 310 # Kelvin
    beta = 1 / (kB * temperature)
    width = len(scoring_matrix)
    energies, motif_z = get_energies(scores, beta, genome_mean_score, width)

    # Calculate probabilities for each energy (being that there is one energy value per site)
    #copy_number = motif.size
    probabilities = [fermi_dirac_probability(energy, genome_z, beta, copy_number) for energy in energies]

    return probabilities, scores, energies


def kl_matrix(PSFMs):
    """Returns a matrix of KLs between every PSFM.
    Ex: KL(PSFMs[0], PSFMs[2]) => KLs[0][2]"""
    KLs = []
    for psfm_p in PSFMs:
        row = []
        for psfm_q in PSFMs:
            row.append(sum([kl(p.values(), q.values()) for p, q in zip(psfm_p, psfm_q)]))
        KLs.append(row)

    return KLs


def save_pickle(data_to_pickle):
    import pickle
    import os
    if not os.path.exists("data/pickles"):
        os.makedirs("data/pickles")
    with open("data/pickles/%s" % data_to_pickle["filename"], "wb") as pickle_file:
        pickle.dump(data_to_pickle, pickle_file)
    debug("Pickled the data to: data/pickles/%s" % data_to_pickle["filename"])


if __name__ == '__main__':
    # Get options from the commandline parameters
    args = parse_args()
    debug("Loaded arguments.")

    # Load the genome
    genome = Genome("../TFBSDB/genomes/%s.fna" % args.genome, FORMAT_FASTA)
    debug("Loaded genome.")

    # Load the original motif
    original_motif = SequenceCollection("motifs/%s.txt" % args.motif, FORMAT_PLAIN)
    debug("Loaded original motif.")

    # Figure out how many iterations to do -- defaults to 1
    num_iterations = int(args.i) if args.i is not None else 1

    # Set the scoring method if it was passed as a parameter
    scoring_method = args.score if args.score is not None else "pssm"

    # Set the copy number if it was passed as a parameter
    copy_number = int(args.cn) if args.cn is not None else original_motif.size

    # Do it!
    PSCMs, Zs, genome_prob_masses = run_iterations(genome, original_motif, copy_number=copy_number, how_many=num_iterations, until_converges=False, scoring_method=scoring_method)

    # Get the probabilities for the original motif for each PSCM
    probabilities, scores, energies = [], [], []
    temperature = 310
    beta = 1 / (kB * temperature)
    for pscm, Z in zip(PSCMs, Zs):
        iteration_scores = []
        iteration_energies = []
        iteration_probabilities = []
        # Calculate scoring matrix based on iteration PSCM (instead of motif PSCM)
        scoring_matrix = create_scoring_matrix(scoring_method, genome, pscm, PSFM(pscm))
        scoring_matrix_wc = wc_matrix(scoring_matrix)
        # Score -> Energy -> Probability for original motif
        for sequence in (site.sequence for site in original_motif):
            score, strand = score_sequence(sequence, scoring_matrix, scoring_matrix_wc)
            #energy = energy_scale(SCORING_INFORMATION_SCORE, score, temperature=temperature, genomic_entropy=genome.entropy, width=len(original_motif[0]))
            energy = TRAP(score, original_motif.width)
            probability = fermi_dirac_probability(energy, Z, beta, copy_number)

            iteration_scores.append(score)
            iteration_energies.append(energy)
            iteration_probabilities.append(probability)
        scores.append(iteration_scores)
        energies.append(iteration_energies)
        probabilities.append(iteration_probabilities)
    debug("Calculated probabilities, scores and energies for the motif for each iteration.")

    # Calculate the probability sum for each iteration
    probability_sums = [sum(iter_probs) for iter_probs in probabilities]
    debug("Calculated probability sum for the original motifs.")

    # Calculate KLs
    PSFMs = [PSFM(pscm) for pscm in PSCMs]
    KLs = kl_matrix(PSFMs)
    debug("Calculated the KLs between each PSFM at each iteration.")

    # Pickle the data and save it to file for analysis
    if args.pickle:
        save_pickle({
            "filename": "%s_%s_%s_CN%d_I%d.pkl" % (args.genome, args.motif, scoring_method, copy_number, num_iterations),
            "genome": args.genome,
            "motif": args.motif,
            "scoring_method": scoring_method,
            "num_iterations": num_iterations,
            "copy_number": copy_number,
            "KLs": KLs,
            "PSCMs": PSCMs,
            "PSFMs": PSFMs,
            "Zs": Zs,
            "probabilities": probabilities,
            "scores": scores,
            "energies": energies,
            "genome_prob_masses": genome_prob_masses
        })

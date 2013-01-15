import math
import re
import time

# Constants
FORMAT_PLAIN = 0
FORMAT_FASTA = 1
FORMAT_GENBANK = 2
FORMAT_CSV = 3
SCORING_PSSM = 0
SCORING_INFORMATION_SCORE = 1
STRAND_ORIGINAL = 0
STRAND_COMPLEMENT = 1
Kb = 0.0019872041  # Units: kcal/mol/K; Boltzmann constant; from http://physics.nist.gov/cuu/Constants/index.html


# Functions
def log2(number):
    """Shorthand for base-2 logarithms"""
    return math.log(number, 2) if number > 0 else math.log(10E-50, 2)


def read_FASTA(FASTAFileName):
    """Returns a dictionary of names/sequences from a FASTA-formatted file"""
    f = open(FASTAFileName)
    sequences = {}
    sequenceLines = []
    for line in f.readlines():
        if line[0] == ">":
            header = line[1:-1]
            if len(sequenceLines) > 0:
                sequences[header] = "".join(sequenceLines)
        else:
            if line[-1] == "\n":
                sequenceLines.append(line[0:-1])
            else:
                sequenceLines.append(line)

    f.close()
    if len(sequenceLines) > 0:
        sequences[header] = "".join(sequenceLines)
    return sequences


def read_plain(fileName):
    """Returns a list of sequences in the format of one sequence per line, skipping any lines that don't contain only nucleotides"""
    sequences = []
    f = open(fileName)
    for line in f.readlines():
        if not re.match("[^ACTG]", line[:-1]):
            sequences.append(line.rstrip())
    return sequences


def wc(sequence):
    """Returns the reverse complement of the sequence; W.C. stands for Watson-Crick complement"""
    return ''.join([{"A": "T", "C": "G", "T": "A", "G": "C"}[base] for base in sequence[::-1]])


def PSFM(sequences, width):
    """Accepts a list of strings. Width is the length of the strings; they should all be the same length
    or width should specify the minimum length"""
    size = float(len(sequences))
    return [{base: ("".join([seq[i] for seq in sequences]).count(base) / size) for base in "ACTG"} for i in range(width)]


def PSSM(motif_PSFM, genomic_frequencies):
    """Accepts PSFM and background genomic frequencies for each base and
    returns a matrix with a score for each base at each position of the motif"""
    return [{base: log2(motif_PSFM[i][base]) - log2(genomic_frequencies[base]) for base in "ACTG"}
        for i in range(len(motif_PSFM))]


def information_score_matrix(motif_PSFM, genomic_entropy):
    return [{base: log2(motif_PSFM[i][base]) + genomic_entropy for base in "ACTG"}
        for i in range(len(motif_PSFM))]


def scorer(sequence, method=SCORING_PSSM, **kwargs):
    if method is SCORING_PSSM or method is SCORING_INFORMATION_SCORE:
        return sum([kwargs["matrix"][i][base] for i, base in enumerate(sequence)])


def sliding_window(sequence, width, step=1):
    """Generator for a sliding window walker; yields a new substring of sequence
    as well as the position of the left end of the window"""
    for position in range(0, len(sequence) - width, step):
        yield position, sequence[position: position + width]


def consensus(motif_PSFM):
    """Returns a string containing to consensus sequence given a frequency matrix"""
    return "".join([max(position, key=position.get) for position in motif_PSFM])


def binding_energy(sequence, motif_PSFM, genomic_frequencies):
    """Calculates the Maxwell-Boltzmann statistics for TF binding energy as derived by Berg and von Hippel"""
    consensus_seq = consensus(motif_PSFM)
    # Maxwell-Boltzmann binding energy (for a single TF molecule):
    #   E = ln(freq(base, motif) * freq(consensus_base, genome)) - ln(freq(consensus_base, motif) * freq(base, genome))
    #   See: Berg and von Hippel (1987), Roider et al. (2006)
    return sum([math.log((motif_PSFM[i][base] * genomic_frequencies[consensus_seq[i]]) /
        (motif_PSFM[i][consensus_seq[i]] * genomic_frequencies[base])) for i, base in enumerate(sequence)])


def Z(motif_sequences, motif_PSFM, genomic_frequencies, beta):
    """Use with binding_probability()
    - beta = 1 / (Kb * temperature) WITH TEMPERATURE IN KELVIN!!"""
    #print "E[0]:", binding_energy(motif_sequences[0], motif_PSFM, genomic_frequencies)
    #print "-beta * E[0]:", -beta * binding_energy(motif_sequences[0], motif_PSFM, genomic_frequencies)
    #print math.e ** (-beta * binding_energy(motif_sequences[0], motif_PSFM, genomic_frequencies))
    # Z = sum(e ^ (-beta * Ei)) where Ei is the binding energy for each sequence in the motif
    return sum([math.e ** (-beta * binding_energy(sequence, motif_PSFM, genomic_frequencies)) for sequence in motif_sequences])


def binding_probability(energy, beta, z):
    """Calculates the probability of binding given a binding energy.
    - binding_energy should be calculated for the sequence (Maxwell-Boltzmann or Fermi-Dirac)
    - Z comes from the Z() function
    - beta = 1 / (Kb * temperature) WITH TEMPERATURE IN KELVIN!!"""
    # P(E) = e ^ (-beta * E) / Z
    #print "e =", math.e
    #print "beta =", beta
    #print "E =",energy
    #print "Z =", z
    return math.e ** (-beta * energy) / z


def energy_weighted_PSFM(motif_sequences, temperature, genomic_frequencies, motif_PSFM, width):
    """Generates a PSFM weighted by binding energy by multiplying the count of each
    base in each sequence by the binding probability of that base at that position"""
    size = float(len(motif_sequences))
    beta = 1 / (Kb * temperature)
    z = Z(motif_sequences, motif_PSFM, genomic_frequencies, beta)
    frequency_matrix = [{base: 0 for base in "ACTG"}] * width
    for sequence in motif_sequences:
        energy = binding_energy(sequence, motif_PSFM, genomic_frequencies)
        p = binding_probability(energy, beta, z)
        for i, base in enumerate(sequence):
            frequency_matrix[i][base] += p  # Where the magic happens
    return frequency_matrix


# Classes
class Genome:
    """Provides methods and attributes for working with a genome sequence"""
    def __init__(self, filename, format=FORMAT_FASTA):
        self.load_file(filename, format)

    def __str__(self):
        """print genomeInstance will output the sequence of the instance of this class"""
        return self.sequence

    def load_file(self, filename, format=FORMAT_FASTA):
        """Formats accepted: FORMAT_FASTA, FORMAT_GENBANK, FORMAT_PLAIN"""
        if format is FORMAT_FASTA:
            FASTA = read_FASTA(filename)
            self.name, self.sequence = FASTA.keys()[0], FASTA.values()[0]
        elif format is FORMAT_GENBANK:
            self.sequence = ""
        elif format is FORMAT_PLAIN:
            self.sequence = ""
        # Store the reverse-complement
        self.wc = wc(self.sequence)
        # Calculate basic metrics
        self.length = len(self.sequence)
        self.counts = {base: self.sequence.count(base) for base in "ACTG"}
        self.frequencies = {base: float(self.counts[base]) / sum(self.counts.values()) for base in "ACTG"}
        # A priori genome entropy, or entropy before binding, calculated by:
        #     Hx = -sum(freq[b] * log_2(freq[b])) for each base b
        self.entropy = -sum([f * log2(f) for f in self.frequencies.values()])


class SequenceCollection:
    """Provides methods for working with a collection, i.e. list of sequences in a genome,
    generally for working with a motif.
    Usage:
        Creating from file:
            >>> motif = SequenceCollection("tf_binding_sites.txt", format=FORMAT_PLAIN)
        Creating from list of strings:
            >>> motif = SequenceCollection(listOfSites)
        Accessing:
            >>> motif.collection[0]
            {"name": "SITE_001", "position": 186839, "score": 0, "regulon": "DevR", "sequence": "ACTGATCGTA"}
    """
    metafields = ["name", "position", "sequence", "strand", "score", "regulon"]

    def __init__(self, sequences, **kwargs):
        """If no format is provided in the kwargs, it is assumed that sequences is a list of string sequences"""
        if "format" in kwargs:
            self.load_file(sequences, kwargs["format"])
        else:
            self.build_collection(sequences, **kwargs)
        self.scored = False

    def __str__(self):
        return "\n".join(["%-8s %s %s" % (item["position"] if item["position"] is not None else "N/A", item["sequence"], item["score"]) for item in self.collection])

    def load_file(self, filename, format=FORMAT_PLAIN):
        """Formats accepted: FORMAT_FASTA, FORMAT_GENBANK, FORMAT_PLAIN"""
        # FASTA: One name per line prefixed with ">"; one sequence per line
        if format == FORMAT_FASTA:
            FASTA = read_FASTA(filename)
            self.build_collection(FASTA.keys(), name=FASTA.values())
        # Plain: One sequence per line
        elif format == FORMAT_PLAIN:
            sequences = read_plain(filename)
            self.build_collection(sequences)
        # TODO: CSV Parsing

    def build_collection(self, sequences, **kwargs):
        """Constructs the self.collection attribute with default value for the metafields if not provided in the arguments."""
        self.collection_size = len(sequences)
        self.width = min([len(seq) for seq in sequences])
        self.collection = []
        # For each sequence in sequences, construct a dictionary with the metafields including available data
        for i, sequence in enumerate(sequences):
            self.collection.append(dict({field: kwargs[field][i] if field in kwargs else None for field in SequenceCollection.metafields}.items() +
                {"sequence": sequence}.items()))
        # Adjust for strandedness
        for i, item in enumerate(self.collection):
            if item["strand"] is STRAND_COMPLEMENT:
                self.collection[i]["sequence"] = wc(item["sequence"])
        # For convenience: expose just the sequences
        self.sequences = [item["sequence"] for item in self.collection]
        self.wc_sequences = [wc(item["sequence"]) for item in self.collection]

    def csv(self):
        return ",".join(SequenceCollection.metafields) + "\n" + "\n".join([",".join([str(item[field]) if item[field] is not None else "" for field in item]) for item in self.collection])


# TODO: class ESA


genome = Genome("../TFBSDB/genomes/NC_000913.fna", FORMAT_FASTA)  # E. coli
motif = SequenceCollection("crp.prodoric.txt", format=FORMAT_PLAIN)  # crp
scoring_method = SCORING_INFORMATION_SCORE
weigh_by_binding_energy = True
temperature = 310  # Kelvin, E. coli optimally grows at 37 C = 310 K
iteration_cutoff = 100

converged = False
iteration = 0
while not converged:
    # Frequency matrices
    frequency_matrix = PSFM(motif.sequences, motif.width)
    wc_frequency_matrix = PSFM(motif.wc_sequences, motif.width)
    # Re-calculate frequencies using binding energy
    if weigh_by_binding_energy:
        frequency_matrix = energy_weighted_PSFM(motif.sequences, temperature, genome.frequencies, frequency_matrix, motif.width)
        wc_frequency_matrix = energy_weighted_PSFM(motif.wc_sequences, temperature, genome.frequencies, wc_frequency_matrix, motif.width)

    # Scoring matrices
    if scoring_method is SCORING_PSSM:
        scoring_matrix = PSSM(frequency_matrix, genome.frequencies)
        wc_scoring_matrix = PSSM(wc_frequency_matrix, genome.frequencies)
    elif scoring_method is SCORING_INFORMATION_SCORE:
        scoring_matrix = information_score_matrix(frequency_matrix, genome.entropy)
        wc_scoring_matrix = information_score_matrix(wc_frequency_matrix, genome.entropy)

    # Score the initial motif in case it isn't already; this is not 100% necessary by the way
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
    # Get the top k scores
    top_k = sorted_genome_scores[:len(motif.collection)]
    # Check for convergence
    if [site[0] for site in top_k] in [item["position"] for item in motif.collection] or iteration >= iteration_cutoff:
        converged = True
    # Logging
    print "\n".join([",".join([str(iteration), site["sequence"], str(site["position"]), str(site["score"])]) for site in motif.collection])
    # Update the motif and iterate
    motif = SequenceCollection([genome.sequence[site[0]:site[0] + motif.width] if strandedness[site[0]] is STRAND_ORIGINAL else wc(genome.sequence[site[0]:site[0] + motif.width]) for site in top_k], position=[site[0] for site in top_k], score=[site[1] for site in top_k], strand=[strandedness[site[0]] for site in top_k])
    iteration += 1

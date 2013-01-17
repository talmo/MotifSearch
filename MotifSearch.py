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
SCORING_WEIGHTED_INFORMATION_SCORE = 2
SCORING_BVH = 3
STRAND_ORIGINAL = 0
STRAND_COMPLEMENT = 1
SCORING_STRINGS = ["PSSM", "Information Scoring", "Weighted Information Scoring", "Berg and von Hippel"]
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


def PSCM(sequences, width):
    """Retuns the count of each base at each position.
    Width is the length of the shortest sequence at most and at minimum 1."""
    return [{base: "".join([seq[i] for seq in sequences]).count(base) for base in "ACTG"} for i in range(width)]


def weighted_PSCM(sequence_pos, width, weights):
    """Retuns the weighted count of each base at each position.
    sequence_pos should be a list of (position, sequence) tuples.
    Width is the length of the shortest sequence at most and at minimum 1.
    Multiplies the counts for the n-th sequence with the n-th weight."""
    count_matrix = [{base: 0 for base in "ACTG"} for i in range(width)]
    for position, sequence in sequence_pos:
        for i, base in enumerate(sequence):
            count_matrix[i][base] += weights[position]
    return count_matrix


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


def weighted_information_score_matrix(motif_PSFM, genomic_entropy):
    return [{base: log2(motif_PSFM[i][base]) * (genomic_entropy + sum([motif_PSFM[i][b] * log2(motif_PSFM[i][b]) for b in "ACTG"])) + genomic_entropy
        for base in "ACTG"} for i in range(len(motif_PSFM))]


def scorer(sequence, method=SCORING_PSSM, **kwargs):
    return sum([kwargs["matrix"][i][base] for i, base in enumerate(sequence)])


def sliding_window(sequence, width, step=1):
    """Generator for a sliding window walker; yields a new substring of sequence
    as well as the position of the left end of the window"""
    for position in range(0, len(sequence) - width, step):
        yield position, sequence[position: position + width]


def consensus(matrix):
    """Returns a string containing to consensus sequence given a frequency or count matrix"""
    return "".join([max(position, key=position.get) for position in matrix])


def berg_von_hippel_matrix(motif_PSCM):
    """Returns a scoring matrix based on the Berg and von Hippel method of site scoring
    to be used for calculation of binding energy and probability."""
    cons_seq = consensus(motif_PSCM)
    return [{base: math.log((motif_PSCM[i][base] + 0.5) / (motif_PSCM[i][cons_seq[i]] + 0.5)) for base in "ACTG"}
        for i in range(len(motif_PSCM))]


def linear_scale(values, worst, best, worst_unscaled):
    return [((worst - best) / worst_unscaled) * value + best for value in values]


def binding_probability(sequence, binding_energy, genome_total_energy, beta):
    """Returns a probability given the binding energy of a site, the binding energies for every site in the genome
    and the beta value, where binding energy is the scaled Berg and von Hippel site score."""
    return math.exp(-beta * binding_energy) / genome_total_energy


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

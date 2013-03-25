"""
A collection of utility functions relating to computational genomics.
"""
import math
import re
import time
import heapq
from collections import namedtuple

# Constants
FORMAT_PLAIN = 0
FORMAT_FASTA = 1
FORMAT_GENBANK = 2
FORMAT_CSV = 3
FORMAT_LIST = 4
FORMAT_DICT = 5
PROBABILITY_METHODS = {"boltzmann": "Boltzmann Probability", "fd": "Fermi-Dirac Probability"}
SCORING_PSSM = 0
SCORING_INFORMATION_SCORE = 1
SCORING_WEIGHTED_INFORMATION_SCORE = 2
SCORING_BVH = 3
SCORING_METHODS = {"pssm": "PSSM", "is": "Information Scoring", "wis": "Weighted Information Scoring", "bvh": "Berg and von Hippel", "wbvh": "Weighted Berg and von Hippel"}
SCORING_STRINGS = ["PSSM", "Information Scoring", "Weighted Information Scoring", "Berg and von Hippel"]
STRAND_ORIGINAL = 0
STRAND_COMPLEMENT = 1
kB = 0.0019872041  # Units: kcal/mol/K; Boltzmann constant; from http://physics.nist.gov/cuu/Constants/index.html


# Functions
def log2(number):
    """Shorthand for base-2 logarithms"""
    return math.log(number, 2) if number > 0 else math.log(10E-50, 2)


def read_FASTA(FASTAFileName):
    """Returns a dictionary of names/sequences from a FASTA-formatted file"""
    # TODO: Double-check some edge cases, i.e. Rosalind -> "Finding a Shared Motif"
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


def wc_matrix(matrix):
    """Returns the reverse complement of a position-weight matrix."""
    return [{"A": position["T"], "T": position["A"], "C": position["G"], "G": position["C"]} for position in matrix[::-1]]


def PSCM(sequences, strand=STRAND_ORIGINAL):
    """Returns the count of each base at each position.
    If STRAND_COMPLEMENT is specified as the second parameter, the reverse-complement of each sequence in sequences is used.
    The counts at position i is multiplied by weights[i] if it exists."""
    width = len(sequences[0])
    # Support for SequenceCollection instances
    if isinstance(sequences, SequenceCollection):
        collection = sequences
        sequences = (site.sequence for site in collection)
        width = len(collection[0].sequence)

    # Initialize count matrix with zeroes
    count_matrix = [{base: 0 for base in "ACTG"} for i in xrange(width)]
    for sequence in sequences:
        # Apply strandedness
        if strand is STRAND_COMPLEMENT:
            sequence = wc(sequence)
        # Count
        for position, base in enumerate(sequence):
            count_matrix[position][base] += 1

    return count_matrix


def PSFM(count_matrix, total=-1):
    """Accepts a PSCM (count matrix) which should be a list of dictionaries with the
    counts for each base at each position.
    Returns a matrix with the frequency of each base at each position."""
    if total == -1:
        total = sum(count_matrix[0].values())  # The total for any given column should be the same
    return [{base: float(count) / total for base, count in position.items()} for position in count_matrix]


def PSSM(motif_PSFM, genomic_frequencies):
    """Accepts PSFM and background genomic frequencies for each base and
    returns a matrix with a score for each base at each position of the motif"""
    return [{base: log2(motif_PSFM[i][base]) - log2(genomic_frequencies[base]) for base in "ACTG"}
        for i in range(len(motif_PSFM))]


def information_score_matrix(motif_PSFM, genomic_entropy):
    """Information scoring is done by doing log(freq) + H(X) <- Genome a priori entropy"""
    return [{base: log2(motif_PSFM[i][base]) + genomic_entropy for base in "ACTG"}
        for i in range(len(motif_PSFM))]


def weighted_information_score_matrix(motif_PSFM, genomic_entropy):
    """Naive weighing method for information scoring -> multiply each position frequency by its information content"""
    return [{base: log2(motif_PSFM[i][base]) * (genomic_entropy + sum([motif_PSFM[i][b] * log2(motif_PSFM[i][b]) for b in "ACTG"])) + genomic_entropy
        for base in "ACTG"} for i in range(len(motif_PSFM))]


def create_scoring_matrix(method, genome, count_matrix, frequency_matrix):
    """Helper function to generalize scoring methods; not very robust
    Returns the scoring matrix given a scoring method"""
    if method in "pssm":
        return PSSM(frequency_matrix, genome.frequencies)
    elif method in "is":
        return information_score_matrix(frequency_matrix, genome.entropy)
    elif method in "bvh":
        return berg_von_hippel_matrix(count_matrix)


def score_sequence(sequence, matrix, matrix_wc):
    """Returns the sum of the score in the matrix for the base at each position of the sequence.
    The sequence (as given) is scored by the scoring matrix given and then by the reverse-complemented
    matrix given and the greater of the two is returned."""
    original = sum([matrix[i][base] for i, base in enumerate(sequence)])
    reverse = sum([matrix_wc[i][base] for i, base in enumerate(sequence)])
    return (original, STRAND_ORIGINAL) if original > reverse else (reverse, STRAND_COMPLEMENT)


def sliding_window(sequence, width=2, step=1):
    """Generator for a sliding window walker; yields a new substring of sequence
    as well as the position of the left end of the window."""
    for position in range(0, len(sequence) - width, step):
        yield position, sequence[position: position + width]


def consensus(matrix):
    """Returns a string containing the consensus sequence of given a frequency or count matrix."""
    return "".join([max(position, key=position.get) for position in matrix])


def berg_von_hippel_matrix(motif_PSCM):
    """Returns a scoring matrix based on the Berg and von Hippel method of site scoring
    to be used for calculation of binding energy and probability.
    See: Berg OG and von Hippel PH (1987) "Selection of DNA binding sites by regulatory proteins. Statistical-mechanical theory and application to operators and promoters.". """
    cons_seq = consensus(motif_PSCM)
    return [{base: math.log((motif_PSCM[i][cons_seq[i]] + 0.5) / (motif_PSCM[i][base] + 0.5)) for base in "ACTG"}
        for i in range(len(motif_PSCM))]


def IS_to_energy(method, score, **kwargs):
    """A theoretical conversion of information score to energy. Upper bounded by 0 => -4kBT.
    kwargs should include temperature, width of the sequence (length), and genomic a priori entropy (H(X))."""
    m = ((2 * kB * kwargs["temperature"] * (2 - kwargs["width"])) / (kwargs["genomic_entropy"] * kwargs["width"]))
    b = -4 * kB * kwargs["temperature"]
    energy = m * score + b
    return energy if score >= 0 else (-4 * kB * kwargs["temperature"])


def TRAP(bvh_score, width):
    """Converts a Berg and von Hippel score to binding energy based on biologically-derived parameters (TRAP).
    See: Manke (2008) and Roider (2007)."""
    lamb = 0.7
    ln_R_0 = 0.585 * width - 5.66
    return (1 / lamb) * bvh_score + ln_R_0

def boltzmann_probability(binding_energy, Z, beta):
    """Returns a binding probability given the binding energy of a site, the partition function(Z)
    and the beta value, where binding energy is measured in kBT."""
    return math.exp(-beta * binding_energy) / Z


def fermi_dirac_probability(binding_energy, Z, beta, copy_number):
    """Returns a binding probability given a binding energy of a site, beta, Z (partition function),
    and copy number (number of TF molecules present). Keep in mind that this based on an APPROXIMATION of mu."""
    mu = (1 / beta) * (math.log(copy_number) - math.log(Z))
    return 1 / (1 + math.exp(beta * (binding_energy - mu)))

def fd_exact(binding_energy, Z, beta, copy_number):
    pass

def kl(p, q):
    """Calculates the Kullback-Leibler diverence for two distributions p and q.
    This is essentially the information gain between two probability distributions."""
    if max(p) > 1 or sum(p) > 1:
        p = normalize(p)
    if max(q) > 1 or sum(q) > 1:
        q = normalize(q)

    # Add pseudo-count to prevent division by 0 errors
    for i in range(len(q)):
        if q[i] == 0:
            q[i] += 10E-50

    return sum([log2(pi / qi) * pi for pi, qi in zip(p, q)])


def normalize(distribution):
    """Normalize the distribution to the range 0-1.0."""
    total = sum(distribution)
    return [p / total for p in distribution]


def top_k_sites(k, sites):
    """Returns the top scoring sites of a collection."""
    return heapq.nlargest(k, scores)


def IC(genomic_entropy, psfm):
    """Returns information content of each position of the PSFM.
        Information content == Mutual information == Motif redundancy index
        IC = H(X) - H(X|Y)  => the amount of information (in bits) required to encode Y given X
        H(X) = a priori entropy => genomic entropy
        H(X|Y) = a posteriori entropy => frequency * log2(frequency)
        """
    return [genomic_entropy + sum([base * log2(base) for base in position]) for position in psfm]


# Classes
class Genome:
    """Provides methods and attributes for working with a genome sequence"""
    def __init__(self, filename, format=FORMAT_FASTA):
        self.load_file(filename, format)

    def __str__(self):
        """print genomeInstance will output the sequence of the instance of this class"""
        return self.sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, k):
        return self.sequence.__getitem__(k)

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


SEQ_FIELDS = ["sequence", "position", "strand", "score", "energy", "probability"]
Seq = namedtuple("Seq", SEQ_FIELDS)


class SequenceCollection:
    """Provides methods for working with a collection, i.e. list of sequences in a genome,
    generally for working with a motif.
    Usage:
        Creating from file:
            >>> motif = SequenceCollection("tf_binding_sites.txt", FORMAT_PLAIN)
        Creating from list of strings:
            >>> motif = SequenceCollection(listOfSites)
        Accessing:
            >>> motif.collection[0]
            {"name": "SITE_001", "position": 186839, "score": 0, "regulon": "DevR", "sequence": "ACTGATCGTA"}
    """

    def __init__(self, collection=[], format=FORMAT_LIST):
        """If no format is provided, it is assumed that sequences is a list"""
        if len(collection) > 0:
            if format is FORMAT_LIST:
                self.build_collection(collection)
            elif format is FORMAT_DICT:
                self.build_collection(collection)
            else:
                self.load_file(collection, format)
        else:
            self.collection = []

    def __str__(self):
        return "\n".join([" ".join([item for item in site]) for site in self.collection])

    def __getitem__(self, k):
        return self.collection.__getitem__(k)

    def load_file(self, filename, format):
        """Formats accepted: FORMAT_FASTA, FORMAT_PLAIN"""
        # FASTA: One name per line prefixed with ">"; one sequence per line
        if format is FORMAT_FASTA:
            FASTA = read_FASTA(filename)
            self.build_collection([{"name": name, "sequence": sequence} for name, sequence in FASTA.items()])
        # Plain: One sequence per line
        elif format is FORMAT_PLAIN:
            sequences = read_plain(filename)
            self.build_collection([[sequence] + [None] * (len(SEQ_FIELDS) - 1) for sequence in sequences])
        # TODO: CSV Parsing

    def build_collection(self, collection):
        """Sets the self.collection attribute and calculates additional information about the collection.
        Important note on strandedness: The ["sequence"] attribute of each sequence must be the STRAND_ORIGINAL
        version! Meaning, if the sequence is 3'-ACTG-5' and it is in the STRAND_COMPLEMENT, then on the actual
        strand, the sequence would be 5'-TGAC-3'.
        **Basically** Make sure to reverse complement any sequence that is on the other strand BEFORE it gets to
        this function."""
        self.collection = [Seq(*site) for site in collection]
        self.size = len(self.collection)
        self.width = len(self.collection[0].sequence)
        self.stats = self.statistics()

    def add(self, site):
        """Converts a dictionary to a namedtuple and stores in the collection."""
        # Populate all fields in case some of them don't exist
        for field in SEQ_FIELDS:
            try:
                site[field]
            except:
                site[field] = None
        self.collection.append(Seq(*[site[field] for field in SEQ_FIELDS]))

    def update(self, i, site):
        """Converts a dictionary to a namedtuple and replaces the element at i in the collection with it."""
        # Populate all fields in case some of them don't exist
        for field in SEQ_FIELDS:
            try:
                site[field]
            except:
                site[field] = getattr(self.collection[i], field)
        self.collection[i] = Seq(*[site[field] for field in SEQ_FIELDS])

    def statistics(self):
        """Calculate some basic statistics on the quantifiable properties of the SequenceCollection."""
        stats = {}
        for category in ["score", "energy", "probability"]:
            if getattr(self.collection[0], category) is not None:
                stats[category] = {subcat: 0 for subcat in ["min", "max", "total", "mean"]}
                stats[category]["min"] = min([getattr(site, category) for site in self.collection])
                stats[category]["max"] = max([getattr(site, category) for site in self.collection])
                stats[category]["total"] = sum([getattr(site, category) for site in self.collection])
                stats[category]["mean"] = float(stats[category]["total"]) / self.size

        # Save as instance variable
        self.stats = stats
        return stats

    def top(self, k, sort_by="score"):
        """Returns the top sites of a collection. sort_by should be the key to sort by, defaults to "score"."""
        return heapq.nlargest(k, self.collection, lambda site: getattr(site, sort_by))

    def collection_generator(self):
        """A generator expression used to avoid loading large collection objects into memory."""
        for i, site in enumerate(self.collection):
            yield i, site

    def sort(self, category, reverse=True):
        """Sorts the collection by the specified category. Generally 'score' or 'probability'."""
        self.collection.sort(reverse=reverse, key=lambda site: getattr(site, category))

    def normalize(self, category):
        """Normalizes the collection by the category specified so that the attribute of all sites add up to 1.0.
        Generally used with 'probability'."""
        total = sum([getattr(site, category) for site in self.collection])
        for i, site in enumerate(self.collection):
            self.update(i, {category: getattr(site, category) / total})

    def truncate(self, prob_sum):
        # Add the probability of the site to the total probability mass until the total is greater than prob_sum
        total, i = 0, 0
        for site in self.collection:
            total += site.probability
            if total > prob_sum:
                break
            i += 1
        # Truncate the collection list
        self.collection = self.collection[:i]
        # Re-calculate statistics if they had been calculated before
        if hasattr(self, "stats"):
            self.statistics()



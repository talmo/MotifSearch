MotifSearch
===========

A simple set of classes and functions to search for search for motifs in a
genome based on positional scoring, i.e. PSSM, as well as iteration for
convergence upon an optimized motif.

### Features:
- Scoring matrix-based scoring:
	- Position Specific Frequency Matrix (PSFM)
	- Position Specific Scoring Matrix (PSSM) scoring (aka PWM)
	- Information Scoring
	- Berg and von Hippel (1987)
- Binding energy and probability based on Berg and von Hippel
- Information statistics (entropy, etc.)
- Sliding window search
- FASTA file parsing

### Soon* to feature:
- Extended Suffix Array sequence searching
- GenBank parsing (likely using BioPython)

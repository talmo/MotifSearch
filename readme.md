MotifSearch
===========

Transcription factors (TFs) bind exhibit both specific and non-specific binding to DNA sequences. In spite of significant research efforts, a precise understanding of the mechanism of TF binding remains elusive.

In an effort to better understand TFs, computational models have been developed to simulate binding. A well-known binding simulation model was described by Marko Djordjevic et al. in a 2003 seminal paper, _A Biophysical Approach to Transcription Factor Binding Site Discovery_. This model allows for calculation of the probability of binding given a few physical parameters and the binding energy of the specific sequence of DNA.

One property of these models is what we have termed "self-consistency". Self-consistency refers to the relative quantity of TF molecules binding to the original starting motif of sites. Essentially, starting off with a set of binding sites, how many of them are recovered and with what probability?

In addition, what will happen if we "evolve" our model by refining the starting motif with the probability that is recovered? In theory, the model should narrow down the probability to a set of sites within the genome that most realistically represents the binding motif that describes a TF.

This program implements this biophysical model of TF binding and iteratively evolves the binding motif by calculating binding probability.

For more information, check out the project's [wiki page](http://erilllab.biosci.umbc.edu/wiki/index.php5/MotifSearch).
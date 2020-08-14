# COVID-19 Preprocess

Preprocess COVID-19 strains data by the algorithm proposed by *[Tempel](https://pubmed.ncbi.nlm.nih.gov/31999330/)*. This is a relatively naive implementation of the original one aiming to test whether time-series training samples construction (explained in *Tempel* Fig.2) can be applied to COVID-19.

# Data Files

`src.paths.file_paths` stores paths to the data which were not uploaded due to the file size limitation of GitHub. The files missing are listed below:
* `covid_amino_acids_path` contains virus types and their corresponding amino acids sequence.
* `covid_metadata_path` contains virus types and their corresponding discovery time.
* `covid_dna_sequences_path` contains virus types and their corresponding DNA sequence.
* `covid_protvec_path` is the [protVec](https://arxiv.org/abs/1503.05140) of amino acids trigrams.
* `covid_overlapping_trigrams_path` contains the embedded sequences represented by numeric vectors.

# Usage

Run `src.execution.setup`.
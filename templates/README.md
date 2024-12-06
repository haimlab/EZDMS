
This program is designed to search for amino acids in specified variable sites and return a list of the combinations of expressed amino acids.

The program accepts a single-line FASTA file and a FASTQ file.

The input FASTA file must be a single-line FASTA file that denotes the variable regions using nonstandard nucleotide symbols. Gaps will be interpreted as variable sites, so please ensure that no gaps are used in the FASTA file unless they denote the variable sites.

The FASTQ file is converted into a FASTA sequence based on the selected Phred score (the default Phred score is 20).

The program searches for variable sites by comparing the nucleotides up to a specified distance on the three-prime and five-prime ends. In cases of a perfect match, it counts the corresponding amino acid. In sequences where all specified amino variable sites are found, the count is added to the output list.

The program returns a CSV file containing each possible amino acid combination and the number of times it was detected.
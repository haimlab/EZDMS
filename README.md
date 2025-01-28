Requirements: 

Python >= 3.7 

 small change

This code searches for variable sites within fastq files and returns the frequency of each combination of amino acids from all variable sites within each sequence. 

Inputs: 

-fasta: The path to a reference fasta file that should contain the complete nucleotide sequence. The expected variable sites must be marked in this file as any non standard nucleotide including N,X,-,?. Anything in the fasta Sequence that is not an A,C,T,G will be detected as a variable site. 

-fastq: The path to the fastq file to be analyzed.  

This function will take in a single dna sequence in fasta format and return the begging site and end site in a variable region.  

Command line options: 

-h, --help show this help message and exit 

-o OUTPUT, --output OUTPUT 

Path to the output file: default="output" 

-p PHREAD, --phread PHREAD 

Threshold value to filliter phread score: default=20 

-5 FIVE_PRIME, --five_prime FIVE_PRIME 

distance from variability sites to the 5-prime end of the required match: default=8 

-3 THREE_PRIME, --three_prime THREE_PRIME 

distance from variability sites to the 3-prime end of the required match: default=8 

-v VARIABLE, --variable VARIABLE 

number of variable sites: default=2 

 

 

 To measure how the selective pressure from Temsavir affects amino acid preference in HIV Env, AD8 plasmids were mutated to contain variable sites at site 375 and 426. These sites were chosen due to their high rate of mutation in response to Temsavir. Overlap PCR was used to generate the two NNK sites, and a homologous barcode sequence was inserted to differentiate WT plasmids from Mutant plasmids. The resulting mutant plasmid library was sequenced using Oxford nano pore technology to generate a fastq file. The fastq library was interpreted using EZDMS Amino Count tool which matches the nucleotides flanking a variable site in barcoded sequences to confirm an even distribution of amino acids at the expected NNK sites.  

293 T Cells were transformed using mutant plasmid library to act as a source of HIV virions. These virions are produced from multiple plasmids and cannot be used to test the fitness of any induvial plasmid strain but were sequence and analyzed with EZDMS Amino ACID count to ensure there was an adequate number of mutant plasmids with evenly distributed amino acids. The first round of transfection has an even distribution of amino acids as multiple plasmids are in each cell allowing unfit genes to be passed on in fit capsids.  

Round two of transfection was performed with a Low MOI viral passage was on A3R5.7 T-cells to dilute the plasmid so that each cell would only being infected once. Each cell only being infected once the produced capsid should match the genotype, they contain allowing for selective pressures to eliminate unfit mutations. Four days after infection the HIV RNA was extracted and sequenced with Oxford nano pore technology.  

Amino Acid Preferences are then calculated using EZDMS utilizing the fastq sequences from round 1 as the preselection library and the fastq sequences from round 2 as the post selection library section library. The Enrichment Ratio  
𝜑𝑝
𝜑
p
 
 for amino acid ‘
𝑦 
y
 
 
’ and a wildtype amino acid of ‘
𝑤𝑡 
w
t
 
 
’ is calculated as the ratio of ‘
𝑦𝑚
y
m
 
’ to ‘
𝑤𝑡𝑚
w
t
m
 
’ after infection divided by the ratio of ‘
𝑦𝑝
y
p
 
’ to ‘
𝑤𝑡𝑝
w
t
p
 
’ before infection.  

The Enrichment Ratio for an amino acid at a single site is rescaled to one based on the sum of all Enrichment Ratios to generate the amino acid preference 
𝜋𝑦
𝜋
y
 
.  

 

 

 

usage: find_variable_sites.py [-h] [-o OUTPUT] [-p PHREAD] [-5 FIVE_PRIME] [-3 THREE_PRIME] [-v VARIABLE] fasta fastq 

Example command line: 

python find_variable_sites.py Ref_375X BNKWKD_3_Library_375X.fastq -o "out.txt" -p 20 -5 8 -3 8 -v 2 

. .venv/bin/activate

 pyinstaller --onefile --windowed --add-data "templates/index.html:templates" --add-data "templates/preference.html:templates" app.py

 pyinstaller --onefile --add-data "templates/:templates" app.py

 pyinstaller --onefile --add-data "templates:templates"  --add-data "modules:modules" app.py 

 pyinstaller --onefile --windowed --add-data "templates:templates" --icon=FGV.icns app.py 
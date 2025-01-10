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

 

 

 

 

 

usage: find_variable_sites.py [-h] [-o OUTPUT] [-p PHREAD] [-5 FIVE_PRIME] [-3 THREE_PRIME] [-v VARIABLE] fasta fastq 

Example command line: 

python find_variable_sites.py Ref_375X BNKWKD_3_Library_375X.fastq -o "out.txt" -p 20 -5 8 -3 8 -v 2 

. .venv/bin/activate

 pyinstaller --onefile --windowed --add-data "templates/index.html:templates" --add-data "templates/preference.html:templates" app.py

 pyinstaller --onefile --add-data "templates/:templates" app.py

 pyinstaller --onefile --add-data "templates:templates"  --add-data "modules:modules" app.py 

 pyinstaller --onefile --windowed --add-data "templates:templates" --icon=FGV.icns app.py 
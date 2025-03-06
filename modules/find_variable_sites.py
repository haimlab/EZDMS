"""
Created on fri Nov 1 2024

@author: Samuel McCarthy-Potter
"""

'''
This code function is to find variable regions in a fasta file. 
This function will take in a single dna sequence in fasta format and return the begging site and end site in a variable region. 

input - .txt .fasta 
a single sequnce fasta file 

output - list of list of sites 
the begning and end site of each varable region 

'''
# Import necessary libraries
from argparse import ArgumentParser # For command-line argument parsing
from re import findall # For regular expression searching
from sys import argv,stderr # For handling command-line arguments and errors
import os

from modules.customErrors import *

def init_argparse():
    """
    Sets up command-line argument parsing for the program.
    
    Allows users to provide input FASTA and FASTQ files, along with output and various configuration options.
    :
    Unkown Error in init_argparse.
    """


    # Define command-line arguments
    parser = ArgumentParser(description ='This program creates a list of the amino acid groups found from specified variability sites')
    parser.add_argument('fasta', type=str, help='Path to the fasta file:\n', default='input.fasta')
    parser.add_argument('fastq', type=str, help='Path to the fastq file:\n', default='fastq.fastq')
    parser.add_argument('-o', '--output', type=str, help='Path to the output file: \tdefault="output"\n', default='output')
    parser.add_argument('-p', '--phread', type=int, help='threshold value to filliter phread score: \tdefault=20\n', default=20)
    parser.add_argument('-5', '--five_prime', type=int, help='distance from variability sites to the 5 prime end of the required match: \tdefault=8\n', default=8)
    parser.add_argument('-3', '--three_prime', type=int, help='distance from variability sites to the 3 prime end of the required match: \tdefault=8\n', default=8)
    parser.add_argument('-v', '--variable', type=int, help='number of variable sites: \tdefault=2\n', default=2)
    
     # If no arguments are passed, show the help message
    if len(argv)==1:
        parser.print_help(stderr)
        
    return parser

def fasta_to_single_line_string(input_fasta,is_file = True):
    """
    Converts a FASTA file into a single-line DNA sequence string.
    
    The sequence is read and non-sequence lines (e.g., headers) are discarded.
    :
    Errors reading the input file. Check that it is the correct fasta format.
    """

    if is_file:
        with open(input_fasta,"r") as open_fasta_file:

            line_list = list(open_fasta_file)

            open_fasta_file.close()
    else:
        line_list = input_fasta.split("\n")

    fasta_sequence = ''

    # Iterate through lines and build the sequence string
    for line in line_list:
        strip_line = line.strip()
        if strip_line [0] == ">":
            # Ensure there is only one sequence in the file
            if len(fasta_sequence) > 0:
                global error_detected
                error_detected = 'There are more then one Sequences in this file'
                print('There are more then one Sequences in this file')
                raise FastaSequenceError('There are more then one Sequences in this file')
        elif strip_line.isalpha():
            fasta_sequence += strip_line.upper() # Only keep alphabetic characters (A, T, C, G)
        else:
            print(f"Unknown character in input fasta {line}")
            raise FastaSequenceError(f"Unknown character in input fasta ({line})")

    return fasta_sequence


def translate_codon(codon):
    """
    Translates a 3-nucleotide codon to its corresponding amino acid.
    
    Uses a predefined genetic code dictionary.
    """

    genetic_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '.', 'TAG': '.',
        'TGT': 'C', 'TGC': 'C', 'TGA': '.', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    return genetic_code.get(codon.upper(), '-')


def reverse_complement(nucleotide):
    """
    Returns the reverse complement of a given nucleotide sequence.
    
    A -> T, T -> A, C -> G, G -> C.
    """

    reverse_complement = []
    for x in nucleotide:
        y = '-'
        if str(x).upper() == "A":
            y = "T"
        if str(x).upper() == "T":
            y = "A"
        if str(x).upper() == "C":
            y = "G"
        if str(x).upper() == "G":
            y = "C"

        reverse_complement.append(y)

    reverse_complement.reverse() # Reverse the order of the sequence
    return(reverse_complement)


def reverse_complement_nucleotide(nucleotide_list,variable_sites_number = 2):
    """
    Reverses and complements a list of nucleotide sequences.
    
    The function handles a list of variable-length regions.
    """

    reverse_nucleotide_list = []

    for item in nucleotide_list:
        if len(item) == variable_sites_number:
            reverse_nucleotide = []
            for i in item:
                reverse_nucleotide += (''.join(reverse_complement(i)),)
            reverse_nucleotide.reverse()
            reverse_nucleotide_list.append(reverse_nucleotide)
    
    return reverse_nucleotide_list

def reverse_complement_region_marker(region_marker_list):
    #Reverses the region marker sequences 

    reverse_complement_region_marker_list = []

    for item in region_marker_list:
        reverse_region_marker = (''.join(reverse_complement(item[1])),''.join(reverse_complement(item[0])),item[-1])
        reverse_complement_region_marker_list.append(reverse_region_marker)

    reverse_complement_region_marker_list.reverse()

    return reverse_complement_region_marker_list

def find_variable_sites(fasta_sequence,variable_sites_number = 2):
    """
    input - .txt .fasta 
    a single sequnce fasta file 

    output - list of list of sites 
    all sites of each varable region 
    :
    Error Variable Sites requested not equal to Variable Sites found.
    """

    #returns a list of all continuous regions with nonstandard nucleotides

    list_non_standard_nucleotide_region = []

    non_standard_nucleotide_site = False

    # Iterate through the sequence to find non-standard nucleotides
    for index,site in enumerate(fasta_sequence):
        if site in ["A","T","C","G"]:
            non_standard_nucleotide_site = False
        else:
            if non_standard_nucleotide_site:
                list_non_standard_nucleotide_region[-1].append(index)
            else:
                list_non_standard_nucleotide_region.append([index])
                non_standard_nucleotide_site = True

    # Check if the number of variable regions matches the expected number
    if len(list_non_standard_nucleotide_region) != variable_sites_number:
        print(list_non_standard_nucleotide_region,variable_sites_number)
        raise VariableSites(f'Variable Sites requested {variable_sites_number} not equal to Variable Sites found {len(list_non_standard_nucleotide_region)}')
 
    return list_non_standard_nucleotide_region

Default_distance_from_region = 8 # Default distance from variability sites to the 5' and 3' regions

def find_guide_sequences(fasta_sequence,list_non_standard_nucleotide_region,distance_5_prime = Default_distance_from_region,distance_3_prime = Default_distance_from_region):
    """
    Identifies regions surrounding variable regions in the given DNA sequence.
    
    This function is used to extract sequences upstream and downstream of variable sites.
    :
    Error Error in Error_detection
    """

    region_marker_list = []

    for region in list_non_standard_nucleotide_region:
        five_prime = fasta_sequence[(region[0] - (distance_5_prime )):region[0]]
        three_prime = fasta_sequence[(region[-1] +1 ):(region[-1] + distance_3_prime+1)]
        length_region = len(region)
        region_marker = (five_prime,three_prime,length_region)
        region_marker_list.append(region_marker)

    return region_marker_list

def find_codon_list(region_marker_list,fasta_line_list,barcode=""):
    """
    :
    Error Number of matching nucleotides invalid.
    """
    #Find guide sequences

    codon_list = []
    for index,fasta_line in enumerate(fasta_line_list):
        sub_codon_list = []
        for region_marker in region_marker_list:
            variable_region = "." * region_marker[-1]
            search_template = f"{region_marker[0] + variable_region + region_marker[1]}"

            look_for = findall(search_template ,fasta_line)

            barcode_found = True
            if barcode != "":
                print('barcode != ""',barcode,findall(barcode ,fasta_line))
                barcode_found = False
                if len(findall(barcode ,fasta_line)) > 0:
                    print(barcode,'barcode_found')
                    barcode_found = True

            

            if len(look_for) > 0 and barcode_found:
                three_prime_end = -len(region_marker[1])
                codon = look_for[0][len(region_marker[0]):three_prime_end]
                
                if three_prime_end == 0:
                    codon = look_for[0][len(region_marker[0]):]

                sub_codon_list.append(codon)



        codon_list.append(sub_codon_list)

    return codon_list

def get_fastq_sequence_list(in_put_fastq,is_file = True):
    """
    :
    Error check if fastq file is valid.
    """
    #Preprocess fastq file into individual list of sequences

    sequence_list  = []
    if is_file:
        with open(in_put_fastq, 'r') as file:
            file_list = list(file)
            file.close()
    else:
        file_list = list(in_put_fastq.split("\n"))

    for site,line in enumerate(file_list):
        line = line.strip()
        if site % 4 == 0:
            sequence_list.append([])
            sequence_list[-1].append(line)
        else:
            sequence_list[-1].append(line)

    if sequence_list[-1] == ['']:
        sequence_list.pop()

    return sequence_list


def fastq_to_fasta_sequence(sequence_list,phread_score = 10):
    """
    :
    Error check if fastq file is valid.
    """
    #Converts the each fastq sequence into a fasta sequence based on phread score 

    list_fasta_sequence = []

    for sequence in sequence_list:
        fasta_string = list(sequence[1])
        for site,char in enumerate(sequence[-1]):
            if ord(char) < 21 + phread_score:
                fasta_string[site] = "-"
        list_fasta_sequence.append([">" + str(sequence[0][1:]),''.join(fasta_string)])

    return [n[1] for n in list_fasta_sequence]


def convert_codons_to_amino_acid_list(codon_list,variable_sites_number):
    """
    :
    Error check if variable sites number and fasta file is valid.
    """
    #This transforms the codon list to a list of amino acids

    amino_acid_list = []

    for codons in codon_list:
        if len(codons) == variable_sites_number:
            amino_pair = []

            for codon in codons:
                amino_pair.append(translate_codon(codon))

            if len(amino_pair) == variable_sites_number and "-" not in amino_pair:
                amino_acid_list.append(amino_pair)

    return amino_acid_list


def load_amino_dic(variable_sites_number=2):
    """
    :
    Error check if variable sites number is valid.
    """
    #Adds the amino acids counts to the dictionary  

    amino_acids = [".",1],["A",2],["C",1],["D",1],["E",1],["F",1],["G",2],["H",1],["I",1],["K",1],["L",3],["M",1],["N",1],["P",2],["Q",1],["R",3],["S",3],["T",2],["V",2],["W",1],["Y",1]

    amino_dic = {}


    loop = range(len(amino_acids) ** variable_sites_number)

    for i in loop:
        key = ""
        count = i

        for j in range(variable_sites_number):

            amino_len = len(amino_acids)
            m = ((count) % (amino_len ** (j+1))) / (amino_len ** (j))
            key += f"{amino_acids[int(m)][0]}\t"
            count -= m
            
        amino_dic[key] = 0

    

    return amino_dic


def populate_amino_dic(amino_acid_list,amino_dic):
    """
    :
    Error in populate_amino_dic
    """
    #Initializes the amino acid dictionary with empty values 

    for amino_pair in amino_acid_list:
        try:
            key = ""
            for amino in amino_pair:
                key += amino + "\t"
            amino_dic [key] += 1
        except:
            raise UnknownError(f"Unknown character in populate_amino_dic {amino_pair}")


def write_out_file(out_file,amino_dic_list,region_marker_list):
    """
    :
    Error in write_out_file
    """
    #Print out file as a excel or CSV

    try:
        assert False
        import pandas as pd
        print('module pandas run')

        out_file = f'{out_file}.xlsx'

        print("Out:",out_file)


        import pandas as pd

        # Create an Excel file with one empty sheet
        with pd.ExcelWriter(out_file, engine='openpyxl', mode='w') as writer:
            # Create an empty DataFrame to add a default sheet
            pd.DataFrame().to_excel(writer, sheet_name="ERROR")
        with pd.ExcelWriter(out_file, engine='openpyxl', mode='a', if_sheet_exists='new') as writer:
            for amino_index,amino_dic in enumerate(amino_dic_list):
                columns=[]

                # Build the column names for the DataFrame
                for index,item in enumerate(region_marker_list):
                    columns.append(f"Site {index+1}:" + str(item[0]))
                columns.append("Count")
                amino = []
                for key in amino_dic.keys():
                    key_list = key.split("\t")
                    key_list[-1] = amino_dic[key]
                    amino.append(key_list)

                df = pd.DataFrame(list(amino), columns=columns)
                df.to_excel(writer, sheet_name=f'barcode_{amino_index+1}', index=False)
            try:
                writer.book.remove(writer.book["ERROR"])
            except:
                pass

        return out_file

    except:
          #Print out file as a excel or CSV
        #amino_dic = amino_dic_list[0]

        out_file = out_file + ".csv"

        print(out_file)

        with open(out_file,"w") as f:
            for index,amino_dic in enumerate(amino_dic_list):
                f.write(f"BarCode_{index+1},")
                for index,item in enumerate(region_marker_list):
                    f.write(f"Site_{index+1}:" + str(item[0]))
                    f.write(",")

                    f.write("Count")

                f.write("\n")

                for key in amino_dic.keys():
                    f.write(",".join(key.split("\t")) +str(amino_dic[key]) +"\n")

                f.write("\n")

            f.close

        return out_file

def Error_out():
    """
    Provides error output if something goes wrong during processing.
    
    This function writes errors to a file called "ERROR.txt".
    :
    Error in Error_out
    """
    
    global error_detected
    out_file = "ERROR.txt"
    print(out_file)
    with open(out_file,"w") as f:
        f.write(error_detected)
    f.close
    return out_file


def main(input_fasta_file,in_put_fastq,out_file = "",is_file_fasta = True,is_file_fastq = True,phread_score = 20,distance_5_prime = 8,distance_3_prime = 8,variable_sites_number = 2,barcode_list=[""]):
    #main driver code 
    """
    Provides error output if something goes wrong during processing.
    
    This function writes errors to a file called "ERROR.txt".
    :
    Error in main driver code
    """

    print("main driver code ")
    print(len(barcode_list))
    print(barcode_list)

    amino_dic_list = []


    for barcode in barcode_list:
        barcode = str(barcode).upper()
        print("here",barcode)


        ref_fasta_sequence = fasta_to_single_line_string(input_fasta_file,is_file_fasta) #Converts fasta file into single line

        list_non_standard_nucleotide_region = find_variable_sites(ref_fasta_sequence,variable_sites_number) #returns a list of all continuous regions with nonstandard nucleotides

        region_marker_list = find_guide_sequences(ref_fasta_sequence,list_non_standard_nucleotide_region,distance_5_prime,distance_3_prime) #Find guide sequences

        sequence_list = get_fastq_sequence_list(in_put_fastq,is_file_fastq) #Preprocess fastq file into individual list of sequences

        fasta_line_list = fastq_to_fasta_sequence(sequence_list,phread_score) #Converts the each fastq sequence into a fasta sequence based on phread score 

        codon_list = find_codon_list(region_marker_list,fasta_line_list,barcode) #Find guide sequences

        amino_acid_list = convert_codons_to_amino_acid_list(codon_list,variable_sites_number) #This transforms the codon list to a list of amino acids

        amino_dic = load_amino_dic(variable_sites_number)  #Initializes the amino acid dictionary with empty values 

        populate_amino_dic(amino_acid_list,amino_dic) #Adds the amino acids counts to the dictionary

        reverse_complement_region_marker_list = reverse_complement_region_marker(region_marker_list) #Reverses the region marker sequences

        reverse_codon_list = find_codon_list(reverse_complement_region_marker_list,fasta_line_list,''.join(reverse_complement(list(barcode)))) #Find guide sequences

        reverse_codon_list = reverse_complement_nucleotide(reverse_codon_list,variable_sites_number) #Reverses the nucleotide sequences

        reverse_amino_acid_list = convert_codons_to_amino_acid_list(reverse_codon_list,variable_sites_number) #This transforms the codon list to a list of amino acids

        populate_amino_dic(reverse_amino_acid_list,amino_dic) #Adds the amino acids counts to the dictionary 

        amino_dic_list.append(amino_dic)

    if len(out_file) > 0:
        print("write_out_file")
        return write_out_file(out_file,amino_dic_list,list_non_standard_nucleotide_region)
    else:
        print("return amino_dic")
        return amino_dic


if __name__ == '__main__':

    parser = init_argparse()
    args = parser.parse_args()

    # Inputs ================================================================================================
    input_fasta_file = args.fasta
    in_put_fastq = args.fastq
    out_file = args.output
    phread_score = args.phread
    five_prime = args.five_prime
    three_prime = args.three_prime
    variable_sites = args.variable
    # ========================================================================================================




    # Check Inputs are valid =================================================================================
    
    assert type(input_fasta_file) == str
    assert type(in_put_fastq) == str
    assert type(out_file) == str
    assert type(phread_score) == int
    assert type(five_prime) == int
    assert type(three_prime) == int
    assert type(variable_sites) == int

    # ========================================================================================================
    

    print(f"Settings selected \nInput fasta file:\t{input_fasta_file}")
    print(f"Input fastq file:\t{in_put_fastq}")
    print(f"Output file name:\t{out_file}")
    print(f"Phread score:\t{phread_score}")
    print(f"Five prime base pair match:\t{five_prime}")
    print(f"Three prime base pair match:\t{three_prime}")
    print(f"Variable sites:\t{variable_sites}")
    print(main(input_fasta_file,in_put_fastq,out_file,True,True,phread_score,five_prime,three_prime,variable_sites))
    
from argparse import ArgumentParser
from re import findall
from sys import argv, stderr
from modules.customErrors import *

# Default distance from variability sites to the 5' and 3' regions
Default_distance_from_region = 8

def init_argparse():
    """Sets up command-line argument parsing for the program."""
    parser = ArgumentParser(description='This program creates a list of the amino acid groups found from specified variability sites')
    parser.add_argument('fasta', type=str, help='Path to the fasta file:', default='input.fasta')
    parser.add_argument('fastq', type=str, help='Path to the fastq file:', default='fastq.fastq')
    parser.add_argument('-o', '--output', type=str, help='Path to the output file (default: "output")', default='output')
    parser.add_argument('-p', '--phread', type=int, help='Phread score threshold (default: 20)', default=20)
    parser.add_argument('-5', '--five_prime', type=int, help='Distance from variability sites to the 5\' end (default: 8)', default=8)
    parser.add_argument('-3', '--three_prime', type=int, help='Distance from variability sites to the 3\' end (default: 8)', default=8)
    parser.add_argument('-v', '--variable', type=int, help='Number of variable sites (default: 2)', default=2)

    if len(argv) == 1:
        parser.print_help(stderr)
    return parser

def fasta_to_single_line_string(input_fasta, is_file=True):
    """Converts a FASTA file into a single-line DNA sequence string."""
    if is_file:
        with open(input_fasta, "r") as open_fasta_file:
            line_list = list(open_fasta_file)
    else:
        line_list = input_fasta.split("\n")

    fasta_sequence = ''
    for line in line_list:
        strip_line = line.strip()
        if strip_line.startswith(">"):
            if len(fasta_sequence) > 0:
                global error_detected
                error_detected = 'There are more than one sequences in this file'
                print(error_detected)
                raise FastaSequenceError(error_detected)
        elif strip_line.isalpha():
            fasta_sequence += strip_line.upper()
        else:
            print(f"Unknown character in input fasta {line}")
            raise FastaSequenceError(f"Unknown character in input fasta ({line})")
    return fasta_sequence

def translate_codon(codon):
    """Translates a 3-nucleotide codon to its corresponding amino acid."""
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
    """Returns the reverse complement of a given nucleotide sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    rev_comp = [complement.get(base.upper(), '-') for base in nucleotide]
    rev_comp.reverse()
    return rev_comp

def reverse_complement_nucleotide(nucleotide_list, variable_sites_number=2):
    """Reverses and complements a list of nucleotide sequences."""
    reverse_nucleotide_list = []
    for item in nucleotide_list:
        if len(item) == variable_sites_number:
            reverse_nucleotide = []
            for i in item:
                reverse_nucleotide.append(''.join(reverse_complement(i)))
            reverse_nucleotide.reverse()
            reverse_nucleotide_list.append(reverse_nucleotide)
    return reverse_nucleotide_list

def reverse_complement_region_marker(region_marker_list):
    """Reverses the region marker sequences."""
    reverse_complement_region_marker_list = []
    for item in region_marker_list:
        reverse_marker = (
            ''.join(reverse_complement(item[1])),
            ''.join(reverse_complement(item[0])),
            item[-1]
        )
        reverse_complement_region_marker_list.append(reverse_marker)
    reverse_complement_region_marker_list.reverse()
    return reverse_complement_region_marker_list

def find_variable_sites(fasta_sequence, variable_sites_number=2):
    """Returns a list of all continuous regions with nonstandard nucleotides."""
    regions = []
    in_region = False
    for idx, base in enumerate(fasta_sequence):
        if base in ["A", "T", "C", "G"]:
            in_region = False
        else:
            if in_region:
                regions[-1].append(idx)
            else:
                regions.append([idx])
                in_region = True
    if len(regions) != variable_sites_number:
        print(regions, variable_sites_number)
        raise VariableSites(f'Variable Sites requested {variable_sites_number} not equal to Variable Sites found {len(regions)}')
    return regions

def find_guide_sequences(fasta_sequence, variable_regions, distance_5_prime=Default_distance_from_region, distance_3_prime=Default_distance_from_region):
    """Identifies regions surrounding variable regions in the given DNA sequence."""
    region_marker_list = []
    for region in variable_regions:
        five_prime = fasta_sequence[(region[0] - distance_5_prime):region[0]]
        three_prime = fasta_sequence[(region[-1] + 1):(region[-1] + distance_3_prime + 1)]
        length_region = len(region)
        region_marker_list.append((five_prime, three_prime, length_region))
    return region_marker_list

def find_codon_list(region_marker_list, fasta_line_list, barcode=""):
    """Finds guide sequences."""
    codon_list = []
    for fasta_line in fasta_line_list:
        sub_codon_list = []
        for region_marker in region_marker_list:
            variable_region = "." * region_marker[-1]
            search_template = f"{region_marker[0]}{variable_region}{region_marker[1]}"
            look_for = findall(search_template, fasta_line)
            barcode_found = True
            if barcode:
                barcode_found = bool(findall(barcode, fasta_line))
            if look_for and barcode_found:
                three_prime_end = -len(region_marker[1])
                codon = look_for[0][len(region_marker[0]):three_prime_end] if three_prime_end != 0 else look_for[0][len(region_marker[0]):]
                sub_codon_list.append(codon)
        codon_list.append(sub_codon_list)
    return codon_list

def get_fastq_sequence_list(in_put_fastq, is_file=True):
    """Preprocesses fastq file into individual list of sequences."""
    sequence_list = []
    if is_file:
        with open(in_put_fastq, 'r') as file:
            file_list = list(file)
    else:
        file_list = in_put_fastq.split("\n")
    for i, line in enumerate(file_list):
        line = line.strip()
        if i % 4 == 0:
            sequence_list.append([line])
        else:
            sequence_list[-1].append(line)
    if sequence_list and sequence_list[-1] == ['']:
        sequence_list.pop()
    return sequence_list

def fastq_to_fasta_sequence(sequence_list, phread_score=10):
    """Converts each fastq sequence into a fasta sequence based on phread score."""
    list_fasta_sequence = []
    for sequence in sequence_list:
        fasta_string = list(sequence[1])
        for site, char in enumerate(sequence[-1]):
            if ord(char) < 21 + phread_score:
                fasta_string[site] = "-"
        list_fasta_sequence.append([">" + str(sequence[0][1:]), ''.join(fasta_string)])
    return [n[1] for n in list_fasta_sequence]

def convert_codons_to_amino_acid_list(codon_list, variable_sites_number):
    """Transforms the codon list to a list of amino acids."""
    amino_acid_list = []
    for codons in codon_list:
        if len(codons) == variable_sites_number:
            amino_pair = [translate_codon(codon) for codon in codons]
            if len(amino_pair) == variable_sites_number and "-" not in amino_pair:
                amino_acid_list.append(amino_pair)
    return amino_acid_list

def load_amino_dic(variable_sites_number=2):
    """Initializes the amino acid dictionary with all possible combinations."""
    amino_acids = [(".", 1), ("A", 2), ("C", 1), ("D", 1), ("E", 1), ("F", 1), ("G", 2), ("H", 1), ("I", 1), ("K", 1), ("L", 3), ("M", 1), ("N", 1), ("P", 2), ("Q", 1), ("R", 3), ("S", 3), ("T", 2), ("V", 2), ("W", 1), ("Y", 1)]
    amino_dic = {}
    amino_len = len(amino_acids)
    for i in range(amino_len ** variable_sites_number):
        key = ""
        count = i
        for j in range(variable_sites_number):
            m = ((count) % (amino_len ** (j + 1))) // (amino_len ** j)
            key += f"{amino_acids[int(m)][0]}\t"
            count -= m * (amino_len ** j)
        amino_dic[key] = 0
    return amino_dic

def populate_amino_dic(amino_acid_list, amino_dic):
    """Adds the amino acids counts to the dictionary."""
    for amino_pair in amino_acid_list:
        try:
            key = "".join([amino + "\t" for amino in amino_pair])
            amino_dic[key] += 1
        except Exception:
            raise UnknownError(f"Unknown character in populate_amino_dic {amino_pair}")

def write_out_file(out_file, amino_dic_list, region_marker_list):
    """Writes output as Excel or CSV."""
    try:
        import pandas as pd
        out_file = f'{out_file}.xlsx'
        with pd.ExcelWriter(out_file, engine='openpyxl', mode='w') as writer:
            pd.DataFrame().to_excel(writer, sheet_name="ERROR")
        with pd.ExcelWriter(out_file, engine='openpyxl', mode='a', if_sheet_exists='new') as writer:
            for amino_index, amino_dic in enumerate(amino_dic_list):
                columns = [f"Amino Acid {i+1}" for i, item in enumerate(region_marker_list)]
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
            except Exception:
                pass
        return out_file
    except Exception:
        out_file = out_file + ".csv"
        with open(out_file, "w") as f:
            for index, amino_dic in enumerate(amino_dic_list):
                f.write(f"Amino Acid_{index+1},")
                for idx, item in enumerate(region_marker_list):
                    f.write(f"Count {idx+1}:" + str(item[0]) + ",")
                f.write("\n")
                for key in amino_dic.keys():
                    f.write(",".join(key.split("\t")) + str(amino_dic[key]) + "\n")
                f.write("\n")
        return out_file

def Error_out():
    """Provides error output if something goes wrong during processing."""
    global error_detected
    out_file = "ERROR.txt"
    with open(out_file, "w") as f:
        f.write(error_detected)
    return out_file

def main(input_fasta_file, in_put_fastq, out_file="", is_file_fasta=True, is_file_fastq=True, phread_score=20, distance_5_prime=8, distance_3_prime=8, variable_sites_number=2, barcode_list=[""]):
    """Main driver code."""
    amino_dic_list = []
    for barcode in barcode_list:
        barcode = str(barcode).upper()
        ref_fasta_sequence = fasta_to_single_line_string(input_fasta_file, is_file_fasta)
        list_non_standard_nucleotide_region = find_variable_sites(ref_fasta_sequence, variable_sites_number)
        region_marker_list = find_guide_sequences(ref_fasta_sequence, list_non_standard_nucleotide_region, distance_5_prime, distance_3_prime)
        sequence_list = get_fastq_sequence_list(in_put_fastq, is_file_fastq)
        fasta_line_list = fastq_to_fasta_sequence(sequence_list, phread_score)
        codon_list = find_codon_list(region_marker_list, fasta_line_list, barcode)
        amino_acid_list = convert_codons_to_amino_acid_list(codon_list, variable_sites_number)
        amino_dic = load_amino_dic(variable_sites_number)
        populate_amino_dic(amino_acid_list, amino_dic)
        reverse_complement_region_marker_list = reverse_complement_region_marker(region_marker_list)
        reverse_codon_list = find_codon_list(reverse_complement_region_marker_list, fasta_line_list, ''.join(reverse_complement(list(barcode))))
        reverse_codon_list = reverse_complement_nucleotide(reverse_codon_list, variable_sites_number)
        reverse_amino_acid_list = convert_codons_to_amino_acid_list(reverse_codon_list, variable_sites_number)
        populate_amino_dic(reverse_amino_acid_list, amino_dic)
        amino_dic_list.append(amino_dic)
    if out_file:
        return write_out_file(out_file, amino_dic_list, list_non_standard_nucleotide_region)
    else:
        return amino_dic

if __name__ == '__main__':
    parser = init_argparse()
    args = parser.parse_args()
    input_fasta_file = args.fasta
    in_put_fastq = args.fastq
    out_file = args.output
    phread_score = args.phread
    five_prime = args.five_prime
    three_prime = args.three_prime
    variable_sites = args.variable

    assert type(input_fasta_file) == str
    assert type(in_put_fastq) == str
    assert type(out_file) == str
    assert type(phread_score) == int
    assert type(five_prime) == int
    assert type(three_prime) == int
    assert type(variable_sites) == int

    print(f"Settings selected \nInput fasta file:\t{input_fasta_file}")
    print(f"Input fastq file:\t{in_put_fastq}")
    print(f"Output file name:\t{out_file}")
    print(f"Phread score:\t{phread_score}")
    print(f"Five prime base pair match:\t{five_prime}")
    print(f"Three prime base pair match:\t{three_prime}")
    print(f"Variable sites:\t{variable_sites}")
    print(main(input_fasta_file, in_put_fastq, out_file, True, True, phread_score, five_prime, three_prime, variable_sites))

import pytest

import find_variable_sites as FVS

#test functions 

global error_detected
error_detected = None


def test_Fasta_string_5():
    assert "FASTA" == FVS.fasta_to_single_line_string("Fasta",False)

def test_find_variable_sites_string_6_1():
    in_put_string = "AAABBBCCCGGGLLLCCC"

    assert [[3, 4, 5], [12, 13, 14]] == FVS.find_variable_sites(in_put_string)

def test_find_variable_sites_string_6_2():
    in_put_string = "AAABBBCCCGGGLLLCCCTTTMMMGGG"

    assert [[3, 4, 5], [12, 13, 14],[21,22,23]] == FVS.find_variable_sites(in_put_string,3)

def test_find_variable_sites_string_6_3():
    in_put_string = "AAABBBCCCGGGCCC"

    assert [[3, 4, 5]] == FVS.find_variable_sites(in_put_string,1)

def test_find_guide_sequences_7_1():
    fasta_sequence = "AAABBBCCCGGGLLLCCC"
    list_non_standard_nucleotide_region = [[3, 4, 5], [12, 13, 14]]
    distance_from_region = 3

    assert [('AAA','CCC',3),('GGG','CCC',3)] == FVS.find_guide_sequences(fasta_sequence,list_non_standard_nucleotide_region, distance_from_region,distance_from_region)

def test_find_guide_sequences_7_2():
    fasta_sequence = "AAABBBCCCGGGLLLCCCTTTMMMGGG"
    list_non_standard_nucleotide_region = [[3, 4, 5], [12, 13, 14]]
    distance_from_region = 3

    assert [('AAA','CCC',3),('GGG','CCC',3)] == FVS.find_guide_sequences(fasta_sequence,list_non_standard_nucleotide_region, distance_from_region,distance_from_region)

def test_find_guide_sequences_7_3():
    fasta_sequence = "AAABBBCCCGGGLLLCCCTTTMMMGGG"
    list_non_standard_nucleotide_region = [[3, 4, 5], [12, 13, 14],[21,22,23]]
    distance_from_region = 3

    assert [('AAA','CCC',3),('GGG','CCC',3),('TTT','GGG',3)] == FVS.find_guide_sequences(fasta_sequence,list_non_standard_nucleotide_region, distance_from_region,distance_from_region)

def test_find_guide_sequences_7_point_4():
    fasta_sequence = "AAABBBCCCGGGLLLCCCT"
    list_non_standard_nucleotide_region = [[3, 4, 5], [12, 13, 14]]
    distance_5_prime = 3
    distance_3_prime = 4

    assert [('AAA','CCCG',3),('GGG','CCCT',3)] == FVS.find_guide_sequences(fasta_sequence,list_non_standard_nucleotide_region, distance_5_prime,distance_3_prime)

def test_find_guide_sequences_7_point_5():
    fasta_sequence = "AAABBBCCCGGGLLLCCCT"
    list_non_standard_nucleotide_region = [[3, 4, 5]]
    distance_5_prime = 3
    distance_3_prime = 3

    assert [('AAA','CCC',3)] == FVS.find_guide_sequences(fasta_sequence,list_non_standard_nucleotide_region, distance_5_prime,distance_3_prime)

def test_find_guide_sequences_7_point_6():
    fasta_sequence = "AAABBBCCCGGGLLLCCCT"
    list_non_standard_nucleotide_region = [[3, 4, 5]]
    distance_5_prime = 3
    distance_3_prime = 2

    assert [('AAA','CC',3)] == FVS.find_guide_sequences(fasta_sequence,list_non_standard_nucleotide_region, distance_5_prime,distance_3_prime)

def test_find_guide_sequences_7_point_7():
    fasta_sequence = "AAABBBCCCGGGLLLCCCT"
    list_non_standard_nucleotide_region = [[3, 4, 5]]
    distance_5_prime = 3
    distance_3_prime = 1

    assert [('AAA','C',3)] == FVS.find_guide_sequences(fasta_sequence,list_non_standard_nucleotide_region, distance_5_prime,distance_3_prime)

def test_find_guide_sequences_7_point_8():
    fasta_sequence = "AAABBBCCCGGGLLLCCCT"
    list_non_standard_nucleotide_region = [[3, 4, 5]]
    distance_5_prime = 3
    distance_3_prime = 0

    assert [('AAA','',3)] == FVS.find_guide_sequences(fasta_sequence,list_non_standard_nucleotide_region, distance_5_prime,distance_3_prime)


def test_find_codon_list_8():
    region_marker_list = [('TTTAAA','CCCGGG',3),('CCCGGG','CCCTTT',3)]
    fasta_line = ["TTTAAAGGGCCCGGGAAACCCTTT","TTTAAATTTCCCGGGTTTCCCTTT"]

    assert [["GGG","AAA"],["TTT","TTT"]] == FVS.find_codon_list(region_marker_list,fasta_line)

def test_find_codon_list_8_1():
    region_marker_list = [('TTTAAA','',3),('CCCGGG','',3)]
    fasta_line = ["TTTAAAGGGCCCGGGAAACCCTTT","TTTAAATTTCCCGGGTTTCCCTTT"]

    assert [["GGG","AAA"],["TTT","TTT"]] == FVS.find_codon_list(region_marker_list,fasta_line)



def test_get_fastq_sequence_list_9():
    in_put_fastq = "@0_GGZBQ6_1\nAAACCCTTTGGG\n+\n((()))((()))\n@0_GGZBQ6_1\nAAACCCTTTGGG\n+\n((()))((()))\n"

    assert [["@0_GGZBQ6_1","AAACCCTTTGGG","+","((()))((()))"],["@0_GGZBQ6_1","AAACCCTTTGGG","+","((()))((()))"]] == FVS.get_fastq_sequence_list(in_put_fastq,False)

def test_fastq_to_fasta_sequence_10():
    sequence_list = [["@0_GGZBQ6_1","AAACCCTTTGGG","+","((()))((()))"],["@0_GGZBQ6_1","AAACCCTTTGGG","+","((()))((()))"]]

    assert ['---CCC---GGG','---CCC---GGG'] == FVS.fastq_to_fasta_sequence(sequence_list,20)

def test_populate_amino_dic_11():
    amino_dic = {".\t.\t":0 , "A\tA\t" : 0}

    amino_acid_list = [[".","."],["A","A"]]

    FVS.populate_amino_dic(amino_acid_list,amino_dic)

    print(amino_dic)

    assert amino_dic == {".\t.\t":1 , "A\tA\t" :1}

def test_reverse_complement_region_marker_12():
    region_marker_list = [('TTTAAA','CCCGGG',3),('CCCGGG','CCCTTT',3)]
    
    assert [('AAAGGG', 'CCCGGG', 3),('CCCGGG', 'TTTAAA', 3)] == FVS.reverse_complement_region_marker(region_marker_list)


def test_convert_codons_to_amino_acid_list_13():
    reverse_codon_list = [['CGA', 'ATT'], [], [], ['CAT', 'CTT']]

    assert [['R', 'I'],['H', 'L']] == FVS.convert_codons_to_amino_acid_list(reverse_codon_list,2)

def test_convert_codons_to_amino_acid_list_13_1():
    reverse_codon_list = [['CGA'], [], [], [ 'CTT']]

    assert [['R'],[ 'L']] == FVS.convert_codons_to_amino_acid_list(reverse_codon_list,1)

def test_convert_codons_to_amino_acid_list_13_2():
    reverse_codon_list = [['CGA', 'ATT','TAG'], [], ['CGA', 'ATT'], ['CAT', 'CTT','TAG']]

    assert [['R', 'I','.'],['H', 'L','.']] == FVS.convert_codons_to_amino_acid_list(reverse_codon_list,3)

def test_translate_codon_14():
    assert "R" == FVS.translate_codon("CGA")

def test_fasta_to_single_line_string_15():
    input_fasta_file = ">DMS_Test_1\ncactttaaatcaagtagctacaaaattaaaagaacaatttgggaataataaaacaatagtctttaatcaatcctcaggaggggacccagaaattgtaatgcacNNKtttaattgtggTggCgaatttttctactgtaattcaacacaactgtttaatagtacttggaattttaatggtacttggaatttaacacaatcgaatggtactgaaggaaatgacactatcacactcccCtgtagaatTaaacaaattataaacNNKtggcaagaagtaggaaaagcaatgtatgcccctcccatcagaggacaaattagatgttcatcaaatattacagggctgatattaacaagagatggtggaaataaccacaataatgataccgagacctttagacctggaggaggagatatgagggacaattggagaagtgaattatataaatataaagtagtaaaaattgaaccattaggagtagcacccaccaaggcaaagagaagagttgtgcagagagaaaaaag"
    FVS.fasta_to_single_line_string(input_fasta_file,False)

def test_find_codon_list_16():

    region_marker_list = [('TCTTGCCA', 'GTTTATAA', 3), ('CAATTAAA', 'GTGCATTA', 3),('CGGTTAAA', 'GTGCACCA', 3)]

    fasta_line_list = ['TCTTGCCACGAGTTTATAACAATTAAAATTGTGCATTACGGTTAAAGGGGTGCACCA']

    Distance_from_region = 8

    assert [["CGA","ATT","GGG"]] == FVS.find_codon_list(region_marker_list,fasta_line_list)

def test_find_codon_list_16_1():

    region_marker_list = [('TTGCCA', 'GTTTATAA', 3), ('CAATTAAA', 'GTGCATTA', 3),('CGGTTAAA', 'GTGCA', 3)]

    fasta_line_list = ['TTGCCACGAGTTTATAACAATTAAAATTGTGCATTACGGTTAAAGGGGTGCA']

    Distance_from_region = 8

    assert [["CGA","ATT","GGG"]] == FVS.find_codon_list(region_marker_list,fasta_line_list)

def test_find_codon_list_16_2():

    region_marker_list = [('TTGCCA', 'GTTTATAA', 3)]

    fasta_line_list = ['TTGCCACGAGTTTATAACAATTAAAATTGTGCATTACGGTTAAAGGGGTGCA']

    Distance_from_region = 8

    assert [["CGA"]] == FVS.find_codon_list(region_marker_list,fasta_line_list)

def test_find_codon_list_16_3():

    region_marker_list = [('', 'GTTTATAA', 3)]

    fasta_line_list = ['TTGCCACGAGTTTATAACAATTAAAATTGTGCATTACGGTTAAAGGGGTGCA']

    Distance_from_region = 8

    assert [["CGA"]] == FVS.find_codon_list(region_marker_list,fasta_line_list)

def test_find_codon_list_16_3():

    region_marker_list = [('TTGCCA', '', 3)]

    fasta_line_list = ['TTGCCACGAGTTTATAACAATTAAAATTGTGCATTACGGTTAAAGGGGTGCA']

    Distance_from_region = 8

    assert [["CGA"]] == FVS.find_codon_list(region_marker_list,fasta_line_list)

def test_convert_codons_to_amino_acid_list_17():
    codon_list = [["CGA","ATT"]]

    assert [['R', 'I']] == FVS.convert_codons_to_amino_acid_list(codon_list,2)

def test_convert_codons_to_amino_acid_list_17_1():
    codon_list = [["TAG"]]

    assert [['.']] == FVS.convert_codons_to_amino_acid_list(codon_list,1)
    
def test_convert_codons_to_amino_acid_list_17_2():
    codon_list = [["CGA","ATT","TAG"]]

    assert [['R', 'I','.']] == FVS.convert_codons_to_amino_acid_list(codon_list,3)

def test_reverse_complement_nucleotide_18():
    reverse_codon_list = [['CGA', 'ATT']]

    assert [['AAT', 'TCG']] == FVS.reverse_complement_nucleotide(reverse_codon_list)

def test_reverse_complement_nucleotide_18_1():
    reverse_codon_list = [['CGA']]

    assert [['TCG']] == FVS.reverse_complement_nucleotide(reverse_codon_list,1)

def test_reverse_complement_nucleotide_18_2():
    reverse_codon_list = [['ATT','CGA', 'ATT']]

    assert [['AAT', 'TCG', 'AAT']] == FVS.reverse_complement_nucleotide(reverse_codon_list,3)
    
def test_load_amino_dic_19():
    for x in range(5):
        variable_sites_number = x
        
        assert 21**x == len(FVS.load_amino_dic(variable_sites_number))

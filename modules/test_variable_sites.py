import pytest

import modules.find_variable_sites as FVS

import modules.Re_scaling as RS

#test functions 

global error_detected
error_detected = None


def test_Fasta_string_5():
    assert "FASTA" == FVS.fasta_to_single_line_string("Fasta",False)

def test_find_variable_sites_string_6_1():
    in_put_string = "AAABBBCCCGGGLLLCCC"

    assert [[3, 4, 5], [12, 13, 14]] == FVS.find_variable_sites(in_put_string,2)

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

def test_main_20():
    input_fasta_file = "sample/Ref_375X426X" 
    in_put_fastq = "static/BNKWKD_3_Library_375X.fastq"
    out_file = "" 
    phread_score = 20 
    five_prime = 8
    three_prime = 8
    variable_sites = 2 

    dic = {'.\t.\t': 0, 'A\t.\t': 0, 'C\t.\t': 0, 'D\t.\t': 0, 'E\t.\t': 0, 'F\t.\t': 0, 'G\t.\t': 0, 'H\t.\t': 0, 'I\t.\t': 0, 'K\t.\t': 0, 'L\t.\t': 0, 'M\t.\t': 0, 'N\t.\t': 0, 'P\t.\t': 0, 'Q\t.\t': 0, 'R\t.\t': 0, 'S\t.\t': 0, 'T\t.\t': 0, 'V\t.\t': 0, 'W\t.\t': 0, 'Y\t.\t': 0, '.\tA\t': 0, 'A\tA\t': 0, 'C\tA\t': 0, 'D\tA\t': 0, 'E\tA\t': 0, 'F\tA\t': 0, 'G\tA\t': 0, 'H\tA\t': 0, 'I\tA\t': 0, 'K\tA\t': 0, 'L\tA\t': 0, 'M\tA\t': 0, 'N\tA\t': 0, 'P\tA\t': 0, 'Q\tA\t': 0, 'R\tA\t': 0, 'S\tA\t': 0, 'T\tA\t': 0, 'V\tA\t': 0, 'W\tA\t': 0, 'Y\tA\t': 0, '.\tC\t': 0, 'A\tC\t': 0, 'C\tC\t': 0, 'D\tC\t': 0, 'E\tC\t': 0, 'F\tC\t': 0, 'G\tC\t': 0, 'H\tC\t': 0, 'I\tC\t': 0, 'K\tC\t': 0, 'L\tC\t': 0, 'M\tC\t': 0, 'N\tC\t': 0, 'P\tC\t': 0, 'Q\tC\t': 0, 'R\tC\t': 0, 'S\tC\t': 0, 'T\tC\t': 0, 'V\tC\t': 0, 'W\tC\t': 0, 'Y\tC\t': 0, '.\tD\t': 0, 'A\tD\t': 0, 'C\tD\t': 0, 'D\tD\t': 0, 'E\tD\t': 0, 'F\tD\t': 0, 'G\tD\t': 0, 'H\tD\t': 0, 'I\tD\t': 0, 'K\tD\t': 0, 'L\tD\t': 0, 'M\tD\t': 0, 'N\tD\t': 0, 'P\tD\t': 0, 'Q\tD\t': 0, 'R\tD\t': 0, 'S\tD\t': 0, 'T\tD\t': 0, 'V\tD\t': 0, 'W\tD\t': 0, 'Y\tD\t': 0, '.\tE\t': 0, 'A\tE\t': 0, 'C\tE\t': 0, 'D\tE\t': 0, 'E\tE\t': 0, 'F\tE\t': 0, 'G\tE\t': 0, 'H\tE\t': 0, 'I\tE\t': 0, 'K\tE\t': 0, 'L\tE\t': 0, 'M\tE\t': 0, 'N\tE\t': 0, 'P\tE\t': 0, 'Q\tE\t': 0, 'R\tE\t': 0, 'S\tE\t': 0, 'T\tE\t': 0, 'V\tE\t': 0, 'W\tE\t': 0, 'Y\tE\t': 0, '.\tF\t': 0, 'A\tF\t': 0, 'C\tF\t': 0, 'D\tF\t': 0, 'E\tF\t': 0, 'F\tF\t': 0, 'G\tF\t': 0, 'H\tF\t': 0, 'I\tF\t': 0, 'K\tF\t': 0, 'L\tF\t': 0, 'M\tF\t': 0, 'N\tF\t': 0, 'P\tF\t': 0, 'Q\tF\t': 0, 'R\tF\t': 0, 'S\tF\t': 0, 'T\tF\t': 0, 'V\tF\t': 0, 'W\tF\t': 0, 'Y\tF\t': 0, '.\tG\t': 0, 'A\tG\t': 0, 'C\tG\t': 0, 'D\tG\t': 0, 'E\tG\t': 0, 'F\tG\t': 0, 'G\tG\t': 0, 'H\tG\t': 0, 'I\tG\t': 0, 'K\tG\t': 0, 'L\tG\t': 0, 'M\tG\t': 0, 'N\tG\t': 0, 'P\tG\t': 0, 'Q\tG\t': 0, 'R\tG\t': 0, 'S\tG\t': 0, 'T\tG\t': 0, 'V\tG\t': 0, 'W\tG\t': 0, 'Y\tG\t': 0, '.\tH\t': 0, 'A\tH\t': 0, 'C\tH\t': 0, 'D\tH\t': 0, 'E\tH\t': 0, 'F\tH\t': 0, 'G\tH\t': 0, 'H\tH\t': 0, 'I\tH\t': 0, 'K\tH\t': 0, 'L\tH\t': 0, 'M\tH\t': 0, 'N\tH\t': 0, 'P\tH\t': 0, 'Q\tH\t': 0, 'R\tH\t': 0, 'S\tH\t': 0, 'T\tH\t': 0, 'V\tH\t': 0, 'W\tH\t': 0, 'Y\tH\t': 0, '.\tI\t': 0, 'A\tI\t': 0, 'C\tI\t': 0, 'D\tI\t': 0, 'E\tI\t': 0, 'F\tI\t': 0, 'G\tI\t': 0, 'H\tI\t': 0, 'I\tI\t': 0, 'K\tI\t': 0, 'L\tI\t': 0, 'M\tI\t': 0, 'N\tI\t': 0, 'P\tI\t': 0, 'Q\tI\t': 0, 'R\tI\t': 0, 'S\tI\t': 0, 'T\tI\t': 0, 'V\tI\t': 0, 'W\tI\t': 0, 'Y\tI\t': 0, '.\tK\t': 0, 'A\tK\t': 0, 'C\tK\t': 0, 'D\tK\t': 0, 'E\tK\t': 0, 'F\tK\t': 0, 'G\tK\t': 0, 'H\tK\t': 0, 'I\tK\t': 0, 'K\tK\t': 0, 'L\tK\t': 0, 'M\tK\t': 0, 'N\tK\t': 0, 'P\tK\t': 0, 'Q\tK\t': 0, 'R\tK\t': 0, 'S\tK\t': 0, 'T\tK\t': 0, 'V\tK\t': 0, 'W\tK\t': 0, 'Y\tK\t': 0, '.\tL\t': 0, 'A\tL\t': 0, 'C\tL\t': 0, 'D\tL\t': 0, 'E\tL\t': 0, 'F\tL\t': 0, 'G\tL\t': 0, 'H\tL\t': 0, 'I\tL\t': 0, 'K\tL\t': 0, 'L\tL\t': 0, 'M\tL\t': 0, 'N\tL\t': 0, 'P\tL\t': 0, 'Q\tL\t': 0, 'R\tL\t': 0, 'S\tL\t': 0, 'T\tL\t': 0, 'V\tL\t': 0, 'W\tL\t': 0, 'Y\tL\t': 0, '.\tM\t': 191, 'A\tM\t': 75, 'C\tM\t': 90, 'D\tM\t': 184, 'E\tM\t': 130, 'F\tM\t': 124, 'G\tM\t': 134, 'H\tM\t': 147, 'I\tM\t': 255, 'K\tM\t': 244, 'L\tM\t': 386, 'M\tM\t': 112, 'N\tM\t': 300, 'P\tM\t': 246, 'Q\tM\t': 178, 'R\tM\t': 264, 'S\tM\t': 347, 'T\tM\t': 507, 'V\tM\t': 162, 'W\tM\t': 73, 'Y\tM\t': 289, '.\tN\t': 0, 'A\tN\t': 0, 'C\tN\t': 0, 'D\tN\t': 0, 'E\tN\t': 0, 'F\tN\t': 0, 'G\tN\t': 0, 'H\tN\t': 0, 'I\tN\t': 0, 'K\tN\t': 0, 'L\tN\t': 0, 'M\tN\t': 0, 'N\tN\t': 0, 'P\tN\t': 0, 'Q\tN\t': 0, 'R\tN\t': 0, 'S\tN\t': 0, 'T\tN\t': 0, 'V\tN\t': 0, 'W\tN\t': 0, 'Y\tN\t': 0, '.\tP\t': 0, 'A\tP\t': 0, 'C\tP\t': 0, 'D\tP\t': 0, 'E\tP\t': 0, 'F\tP\t': 0, 'G\tP\t': 0, 'H\tP\t': 0, 'I\tP\t': 0, 'K\tP\t': 0, 'L\tP\t': 0, 'M\tP\t': 0, 'N\tP\t': 0, 'P\tP\t': 0, 'Q\tP\t': 0, 'R\tP\t': 0, 'S\tP\t': 0, 'T\tP\t': 0, 'V\tP\t': 0, 'W\tP\t': 0, 'Y\tP\t': 0, '.\tQ\t': 0, 'A\tQ\t': 0, 'C\tQ\t': 0, 'D\tQ\t': 0, 'E\tQ\t': 0, 'F\tQ\t': 0, 'G\tQ\t': 0, 'H\tQ\t': 0, 'I\tQ\t': 0, 'K\tQ\t': 0, 'L\tQ\t': 0, 'M\tQ\t': 0, 'N\tQ\t': 0, 'P\tQ\t': 0, 'Q\tQ\t': 0, 'R\tQ\t': 0, 'S\tQ\t': 0, 'T\tQ\t': 0, 'V\tQ\t': 0, 'W\tQ\t': 0, 'Y\tQ\t': 0, '.\tR\t': 0, 'A\tR\t': 0, 'C\tR\t': 0, 'D\tR\t': 0, 'E\tR\t': 0, 'F\tR\t': 0, 'G\tR\t': 0, 'H\tR\t': 0, 'I\tR\t': 0, 'K\tR\t': 0, 'L\tR\t': 0, 'M\tR\t': 0, 'N\tR\t': 0, 'P\tR\t': 0, 'Q\tR\t': 0, 'R\tR\t': 0, 'S\tR\t': 0, 'T\tR\t': 0, 'V\tR\t': 0, 'W\tR\t': 0, 'Y\tR\t': 0, '.\tS\t': 0, 'A\tS\t': 0, 'C\tS\t': 0, 'D\tS\t': 0, 'E\tS\t': 0, 'F\tS\t': 0, 'G\tS\t': 0, 'H\tS\t': 0, 'I\tS\t': 0, 'K\tS\t': 0, 'L\tS\t': 0, 'M\tS\t': 0, 'N\tS\t': 0, 'P\tS\t': 0, 'Q\tS\t': 0, 'R\tS\t': 0, 'S\tS\t': 0, 'T\tS\t': 0, 'V\tS\t': 0, 'W\tS\t': 0, 'Y\tS\t': 0, '.\tT\t': 0, 'A\tT\t': 0, 'C\tT\t': 0, 'D\tT\t': 0, 'E\tT\t': 0, 'F\tT\t': 0, 'G\tT\t': 0, 'H\tT\t': 0, 'I\tT\t': 0, 'K\tT\t': 0, 'L\tT\t': 0, 'M\tT\t': 0, 'N\tT\t': 0, 'P\tT\t': 0, 'Q\tT\t': 0, 'R\tT\t': 0, 'S\tT\t': 0, 'T\tT\t': 0, 'V\tT\t': 0, 'W\tT\t': 0, 'Y\tT\t': 0, '.\tV\t': 0, 'A\tV\t': 0, 'C\tV\t': 0, 'D\tV\t': 0, 'E\tV\t': 0, 'F\tV\t': 0, 'G\tV\t': 0, 'H\tV\t': 0, 'I\tV\t': 0, 'K\tV\t': 0, 'L\tV\t': 0, 'M\tV\t': 0, 'N\tV\t': 0, 'P\tV\t': 0, 'Q\tV\t': 0, 'R\tV\t': 0, 'S\tV\t': 0, 'T\tV\t': 0, 'V\tV\t': 0, 'W\tV\t': 0, 'Y\tV\t': 0, '.\tW\t': 0, 'A\tW\t': 0, 'C\tW\t': 0, 'D\tW\t': 0, 'E\tW\t': 0, 'F\tW\t': 0, 'G\tW\t': 0, 'H\tW\t': 0, 'I\tW\t': 0, 'K\tW\t': 0, 'L\tW\t': 0, 'M\tW\t': 0, 'N\tW\t': 0, 'P\tW\t': 0, 'Q\tW\t': 0, 'R\tW\t': 0, 'S\tW\t': 0, 'T\tW\t': 0, 'V\tW\t': 0, 'W\tW\t': 0, 'Y\tW\t': 0, '.\tY\t': 0, 'A\tY\t': 0, 'C\tY\t': 0, 'D\tY\t': 0, 'E\tY\t': 0, 'F\tY\t': 0, 'G\tY\t': 0, 'H\tY\t': 0, 'I\tY\t': 0, 'K\tY\t': 0, 'L\tY\t': 0, 'M\tY\t': 0, 'N\tY\t': 0, 'P\tY\t': 0, 'Q\tY\t': 0, 'R\tY\t': 0, 'S\tY\t': 0, 'T\tY\t': 0, 'V\tY\t': 0, 'W\tY\t': 0, 'Y\tY\t': 0}

    pre_amino_dict = FVS.main(input_fasta_file,in_put_fastq,out_file,True,True,phread_score,five_prime,three_prime,variable_sites)
    assert pre_amino_dict == dic

def test_find_variable_sites_1():

    #the enrichment ratio, which is the relative frequency of mutations to x after selection versus before selection

    urX = 100 #pre mutation
    urWT = 100 #pre mutation

    frX = 10 #post mutation
    frWT = 50 #post mutation

    #assert negative selection
    assert .2 == RS.calculateEnrichmentRatio(frX,frWT,urX,urWT)

def test_find_variable_sites_2():
    urX = 100 #pre mutation
    urWT = 100 #pre mutation

    frX = 50 #post mutation
    frWT = 50 #post mutation


    #assert netrual selection
    assert 1 == RS.calculateEnrichmentRatio(frX,frWT,urX,urWT)

def test_find_variable_sites_3():
    urX = 100 #pre mutation
    urWT = 100 #pre mutation

    frX = 50 #post mutation
    frWT = 10 #post mutation

    #assert positive selection
    assert 5 == RS.calculateEnrichmentRatio(frX,frWT,urX,urWT)

def test_find_variable_sites_4():
    EnrichmentRatioDict = {1:1}
    PreferenceSite = 1
    assert 1 == RS.calculatePreference(EnrichmentRatioDict,PreferenceSite)

def test_find_variable_sites_5():
    EnrichmentRatioDict = {1:1,2:1,3:1,4:1}
    PreferenceSite = 1
    assert .25 == RS.calculatePreference(EnrichmentRatioDict,PreferenceSite)

def test_find_variable_sites_6():
    WT ="1" 
    pre_amino_dict = {"1\t":100,"2\t":100,"3\t":100,"4\t":100,"5\t":100,"6\t":100,"7\t":100,"8\t":100}
    post_amino_dict = {"1\t":50,"2\t":10,"3\t":80,"4\t":100,"5\t":00,"6\t":30,"7\t":5,"8\t":100}
    expected_dict = {'1\t': 0.13333333333333333, '2\t': 0.02666666666666667, '3\t': 0.21333333333333335, '4\t': 0.26666666666666666, '5\t': 0.0, '6\t': 0.08, '7\t': 0.013333333333333334, '8\t': 0.26666666666666666}
    print(RS.main(WT,pre_amino_dict,post_amino_dict,""))



def calculateReScaledEnrichmentRatio(SiReported,SMedianSynonymous,SMedianBottom):
    SiScaled = ((SiReported - SMedianSynonymous)/-SMedianBottom)+1
    return SiScaled

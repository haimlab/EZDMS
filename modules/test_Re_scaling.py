from modules.Re_scaling import main

import modules.find_variable_sites as FVS

if __name__ == '__main__':
    input_fasta_file = 'static/Ref_375X.fa' 
    in_put_fastq = "static/S375X_Lib_1_Viral Stock_RT_cDNA_Extra PCR.fastq"
    in_put_fastq2 = "static/H2MBPL_1_Lib1_S375X_A3R5.7_Non-Int_0.001_PCRx2.fastq"
    out_file = "" 
    phread_score = 20 
    five_prime = 8
    three_prime = 8
    variable_sites = 1

    pre_amino_dict = FVS.main(input_fasta_file,in_put_fastq,out_file,True,True,phread_score,five_prime,three_prime,variable_sites,["ggTggCga"])
    post_amino_dict = FVS.main(input_fasta_file,in_put_fastq2,out_file,True,True,phread_score,five_prime,three_prime,variable_sites,["ggTggCga"])
    print(main("S",pre_amino_dict,post_amino_dict,"here",1))
    '''
    print()
    WT ="1" 
    pre_amino_dict = {"1\t":100,"2\t":100,"3\t":100,"4\t":100,"5\t":100,"6\t":100,"7\t":100,"8\t":100}
    post_amino_dict = {"1\t":50,"2\t":10,"3\t":80,"4\t":100,"5\t":00,"6\t":30,"7\t":5,"8\t":100}
    expected_dict = {'1\t': 0.13333333333333333, '2\t': 0.02666666666666667, '3\t': 0.21333333333333335, '4\t': 0.26666666666666666, '5\t': 0.0, '6\t': 0.08, '7\t': 0.013333333333333334, '8\t': 0.26666666666666666}
    assert expected_dict == main(WT,pre_amino_dict,post_amino_dict,"")
    '''


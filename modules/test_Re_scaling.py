from modules.Re_scaling import main

import modules.find_variable_sites as FVS

if __name__ == '__main__':
    input_fasta_file = 'sample/Ref_375X.fa' 
    in_put_fastq = "/Volumes/rdss_hhaim/LAB PROJECTS/Sam/FQ_VSearch/sample/BNKWKD_3_Library_375X.fastq"
    in_put_fastq2 = "sample/GGZBQ6_1_DMS_DX_1.fastq"
    out_file = "" 
    phread_score = 20 
    five_prime = 8
    three_prime = 8
    variable_sites = 1

    #pre_amino_dict = FVS.main(input_fasta_file,in_put_fastq,out_file,True,True,phread_score,five_prime,three_prime,variable_sites)
    #post_amino_dict = FVS.main(input_fasta_file,in_put_fastq2,out_file,True,True,phread_score,five_prime,three_prime,variable_sites)
    #print(main("S",pre_amino_dict,post_amino_dict,"here"))
    print()
    WT ="1" 
    pre_amino_dict = {"1\t":100,"2\t":100,"3\t":100,"4\t":100,"5\t":100,"6\t":100,"7\t":100,"8\t":100}
    post_amino_dict = {"1\t":50,"2\t":10,"3\t":80,"4\t":100,"5\t":00,"6\t":30,"7\t":5,"8\t":100}
    expected_dict = {'1\t': 0.13333333333333333, '2\t': 0.02666666666666667, '3\t': 0.21333333333333335, '4\t': 0.26666666666666666, '5\t': 0.0, '6\t': 0.08, '7\t': 0.013333333333333334, '8\t': 0.26666666666666666}
    assert expected_dict == main(WT,pre_amino_dict,post_amino_dict,"")


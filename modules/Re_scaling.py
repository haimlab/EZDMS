import modules.find_variable_sites as FVS

import numpy as np

def calculateReScaledEnrichmentRatio(SiReported,SMedianSynonymous,SMedianBottom):
    SiScaled = ((SiReported - SMedianSynonymous)/-SMedianBottom)+1
    return SiScaled

def calculateEnrichmentRatio(frX,frWT,urX,urWT):
    EnrichmentRatio = ((frX/frWT)/(urX/urWT))
    return EnrichmentRatio

def calculatePreference(EnrichmentRatioDict,PreferenceSite):
    EnrichmentRatioList = []
    for key in EnrichmentRatioDict.keys():
        EnrichmentRatioList.append(float(EnrichmentRatioDict[key]))
    Preference = EnrichmentRatioDict[PreferenceSite]/sum(EnrichmentRatioList)
    return Preference

def write_out_file(out_file,amino_dic,region_marker_list):
    """
    :
    Error in write_out_file
    """
    #Print out file as a excel or CSV

    out_file = out_file + ".csv"

    print(out_file)

    with open(out_file,"w") as f:
        for index,item in enumerate(region_marker_list):
            f.write(f"Site_{index+1}:" + str(item[0]))
            f.write(",")

        f.write("amino_acid_counter")
        f.write("\n")

        for key in amino_dic .keys():
            f.write(",".join(key.split("\t")) +str(amino_dic[key]) +"\n")

    f.close

    return out_file

def main(pre_amino_dict,post_amino_dict,out_file = "1"):
    for key in pre_amino_dict.keys():
        pre_amino_dict[key] = pre_amino_dict[key]/2
    pre_amino_dict['S\t'] = 9000

    pre_amino_dict_WT = {}
    for key in pre_amino_dict.keys():
        pre_amino_dict_WT[key] = 100
    pre_amino_dict_WT['S\t'] = 9000

    post_amino_dict_WT = {}
    for key in pre_amino_dict.keys():
        post_amino_dict_WT[key] = 11
    post_amino_dict_WT['S\t'] = 2000

    WT = 'S\t'

    pre_amino_list = []
    for key in pre_amino_dict.keys():
        pre_amino_list.append(float(pre_amino_dict[key]))
    pre_amino_sum = sum(pre_amino_list)

    post_amino_list = []
    for key in post_amino_dict.keys():

        post_amino_list.append(float(post_amino_dict[key]))
    post_amino_sum = sum(post_amino_list)

    pre_amino_WT_list = []
    for key in post_amino_dict_WT.keys():
        pre_amino_WT_list.append(float(pre_amino_dict_WT[key]))
    pre_amino_WT_sum = sum(pre_amino_WT_list)

    post_amino_WT_list = []
    for key in post_amino_dict_WT.keys():
        post_amino_WT_list.append(float(post_amino_dict_WT[key]))
    post_amino_WT_sum = sum(post_amino_WT_list)

    Enrichment_dict = {}
    for key in pre_amino_dict.keys():
        Enrichment_dict[key] = float(calculateEnrichmentRatio(pre_amino_dict[key]/pre_amino_sum,pre_amino_dict_WT[WT]/pre_amino_WT_sum ,post_amino_dict[key]/post_amino_sum,post_amino_dict_WT[WT]/post_amino_WT_sum))
        
    Preference_dict = {}
    for key in Enrichment_dict.keys():
        Preference_dict[key] = calculatePreference(Enrichment_dict,key)

    post_amino_WT_list = []
    for key in post_amino_dict_WT.keys():
        post_amino_WT_list.append(float(post_amino_dict_WT[key]))

    ReScaledEnrichmentRatio = {}
    for key in Preference_dict.keys():

        SiReported = Preference_dict[key]

        SMedianSynonymous = Preference_dict[WT]

        SMedianBottom = np.percentile(post_amino_WT_list, 1)

        ReScaledEnrichmentRatio[key] = calculateReScaledEnrichmentRatio(SiReported,SMedianSynonymous,SMedianBottom)

    if len(out_file) > 0:
        print("run")
        return write_out_file(out_file, Preference_dict,[[1]])
    else:
        return Preference_dict

if __name__ == '__main__':
    input_fasta_file = '/Volumes/rdss_hhaim/LAB PROJECTS/Sam/FQ_VSearch/sample/Ref_"375X"426X' #"ttaatcaatcctcaggaggggacccagaaattgtaatgcacNNKtttaattgtggTggCgaatttttctactgtaattcaacacaact" 
    in_put_fastq = "BNKWKD_2_A3R5.7_375X426X.fastq"
    out_file = "" 
    phread_score = 20 
    five_prime = 8
    three_prime = 8
    variable_sites = 1

    pre_amino_dict = FVS.main(input_fasta_file,in_put_fastq,out_file,True,True,phread_score,five_prime,three_prime,variable_sites)
    post_amino_dict = FVS.main(input_fasta_file,in_put_fastq,out_file,True,True,phread_score,five_prime,three_prime,variable_sites)
    print(main(pre_amino_dict,post_amino_dict,"here"))



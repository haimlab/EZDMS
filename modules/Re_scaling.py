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

def main(WT_list,pre_amino_dict,post_amino_dict,out_file = "1"):

    WT = "\t".join(WT_list.split(","))+"\t"

    for key in pre_amino_dict.keys():
        pre_amino_dict[key] = pre_amino_dict[key]/2

    pre_amino_list = []
    for key in pre_amino_dict.keys():
        pre_amino_list.append(float(pre_amino_dict[key]))
    pre_amino_sum = sum(pre_amino_list)

    post_amino_list = []
    for key in post_amino_dict.keys():
        post_amino_list.append(float(post_amino_dict[key]))
    post_amino_sum = sum(post_amino_list)

    Enrichment_dict = {}
    print(pre_amino_dict[WT])
    pre_amino_wild_type = pre_amino_dict[WT]/pre_amino_sum
    if pre_amino_wild_type == 0:
        pre_amino_wild_type = 1
    post_amino_wild_type = post_amino_dict[WT]/post_amino_sum
    if post_amino_wild_type == 0:
        post_amino_wild_type= 1
    for key in pre_amino_dict.keys():
        try:

            pre_amino_ratio = pre_amino_dict[key]/pre_amino_sum
            post_amino_ratio = post_amino_dict[key]/post_amino_sum
            
            Enrichment_dict[key] = float(calculateEnrichmentRatio(post_amino_ratio,post_amino_wild_type,pre_amino_ratio,pre_amino_wild_type))
            #print(key,Enrichment_dict[key])

        except ZeroDivisionError:
            Enrichment_dict[key] = 0

    Preference_dict = {}
    for key in Enrichment_dict.keys():
        Preference_dict[key] = calculatePreference(Enrichment_dict,key)



    """ReScaledEnrichmentRatio = {}

    for key in Preference_dict.keys():

        SiReported = Preference_dict[key]

        SMedianSynonymous = Preference_dict[WT]

        SMedianBottom = np.percentile(post_amino_list, 1)

        ReScaledEnrichmentRatio[key] = calculateReScaledEnrichmentRatio(SiReported,SMedianSynonymous,SMedianBottom) 
        """

    print(out_file,len(out_file))
    if len(out_file) > 0:
        print("run")
        return write_out_file(out_file, Preference_dict,[[1]])
    else:
        return Preference_dict

if __name__ == '__main__':
    input_fasta_file = 'Work_flow/Reference.fa' 
    in_put_fastq = "/Users/smccarthypotter/FQ_VSearch/Work_flow/Plasmid Library_S375X.fastq"
    in_put_fastq2 = "sample/GGZBQ6_1_DMS_DX_1.fastq"
    out_file = "" 
    phread_score = 20 
    five_prime = 8
    three_prime = 8
    variable_sites = 1

    WT = 'S'


    pre_amino_dict = FVS.main(input_fasta_file,in_put_fastq,out_file,True,True,phread_score,five_prime,three_prime,variable_sites)
    post_amino_dict = FVS.main(input_fasta_file,in_put_fastq2,out_file,True,True,phread_score,five_prime,three_prime,variable_sites)
    print(main(WT,pre_amino_dict,post_amino_dict,"here"))





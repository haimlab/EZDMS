import modules.find_variable_sites as FVS
import numpy as np
import pandas as pd

def calculateReScaledEnrichmentRatio(SiReported, SMedianSynonymous, SMedianBottom):
    """Calculates the rescaled enrichment ratio."""
    SiScaled = ((SiReported - SMedianSynonymous) / -SMedianBottom) + 1
    return SiScaled

def calculateEnrichmentRatio(frX, frWT, urX, urWT):
    """Calculates the enrichment ratio."""
    return (frX / frWT) / (urX / urWT)

def calculatePreference(EnrichmentRatioDict, PreferenceSite):
    """Calculates the preference for a given site."""
    EnrichmentRatioList = [float(val) for val in EnrichmentRatioDict.values()]
    Preference = EnrichmentRatioDict[PreferenceSite] / sum(EnrichmentRatioList)
    return Preference

def write_out_file(out_file, amino_dic, region_marker_list, add_column=False):
    """
    Writes the output to an Excel (.xlsx) file with proper headers.
    - If three columns: Amino Acid 1, Amino Acid 2, Preference
    - If two columns: Amino Acid 1, Preference
    """
    out_file = out_file + ".xlsx"

    # Prepare data for DataFrame
    data = []
    for key in amino_dic.keys():
        split_key = key.split("\t")
        split_key = [aa for aa in split_key if aa]  # Remove empty strings
        if add_column:
            row = split_key + [amino_dic[key][0], amino_dic[key][1]]
        else:
            row = split_key + [amino_dic[key]]
        data.append(row)

    # Determine columns based on number of amino acids/sites
    n_sites = len(data[0]) - (2 if add_column else 1)
    if n_sites == 2:
        columns = ["Amino Acid 1", "Amino Acid 2", "Preference"]
    elif n_sites == 1:
        columns = ["Amino Acid 1", "Preference"]
    else:
        columns = [f"Amino Acid {i+1}" for i in range(n_sites)] + ["Preference"]

    df = pd.DataFrame(data, columns=columns)

    # Write to Excel
    df.to_excel(out_file, index=False)
    print(f"Output written to {out_file}")
    return out_file

def main(WT_list, pre_amino_dict, post_amino_dict, out_file="1", add_column=False):
    """Main driver for enrichment and preference calculations."""
    WT = "\t".join(WT_list.split(",")) + "\t"

    # Halve all pre-amino values
    for key in pre_amino_dict.keys():
        pre_amino_dict[key] = pre_amino_dict[key] / 2

    pre_amino_list = [float(val) for val in pre_amino_dict.values()]
    pre_amino_sum = sum(pre_amino_list)

    post_amino_list = [float(val) for val in post_amino_dict.values()]
    post_amino_sum = sum(post_amino_list)

    Enrichment_dict = {}
    pre_amino_wild_type = pre_amino_dict.get(WT, 0) / pre_amino_sum if pre_amino_sum else 1
    if pre_amino_wild_type == 0:
        pre_amino_wild_type = 1
    post_amino_wild_type = post_amino_dict.get(WT, 0) / post_amino_sum if post_amino_sum else 1
    if post_amino_wild_type == 0:
        post_amino_wild_type = 1

    for key in pre_amino_dict.keys():
        try:
            pre_amino_ratio = pre_amino_dict[key] / pre_amino_sum if pre_amino_sum else 0
            post_amino_ratio = post_amino_dict[key] / post_amino_sum if post_amino_sum else 0
            Enrichment_dict[key] = float(calculateEnrichmentRatio(post_amino_ratio, post_amino_wild_type, pre_amino_ratio, pre_amino_wild_type))
        except ZeroDivisionError:
            Enrichment_dict[key] = 0

    Preference_dict = {}
    if add_column:
        for key in Enrichment_dict.keys():
            preference = calculatePreference(Enrichment_dict, key)
            Preference_dict[key] = [preference, preference * (len(Enrichment_dict) / 1)]
    else:
        for key in Enrichment_dict.keys():
            Preference_dict[key] = calculatePreference(Enrichment_dict, key)

    if out_file:
        return write_out_file(out_file, Preference_dict, [[1]], add_column)
    else:
        return Preference_dict

if __name__ == '__main__':
    input_fasta_file = 'Work_flow/Reference.fa'
    in_put_fastq = "/Users/smccarthypotter/FQ_VSearch/Work_flow/Plasmid Library_S375X.fastq"
    in_put_fastq2 = "sample/GGZBQ6_1_DMS_DX_1.fastq"
    out_file = "output"
    phread_score = 20
    five_prime = 8
    three_prime = 8
    variable_sites = 1

    WT = 'S'

    pre_amino_dict = FVS.main(input_fasta_file, in_put_fastq, out_file, True, True, phread_score, five_prime, three_prime, variable_sites)
    post_amino_dict = FVS.main(input_fasta_file, in_put_fastq2, out_file, True, True, phread_score, five_prime, three_prime, variable_sites)
    print(main(WT, pre_amino_dict, post_amino_dict, "here"))

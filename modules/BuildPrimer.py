from modules.find_variable_sites import fasta_to_single_line_string,reverse_complement,find_variable_sites

def main(input_fasta_file,Out_path):

    Out_path = Out_path + ".txt"


    is_file_fasta = True

    print()

    var_site = (find_variable_sites(fasta_to_single_line_string(input_fasta_file,True),1))

    print(var_site)

    var_site_s = var_site[0][0]
    var_site_e = var_site[0][-1]+1

    ref_fasta_sequence = fasta_to_single_line_string(input_fasta_file,is_file_fasta)

    print(ref_fasta_sequence[var_site_s-40:var_site_s-30],ref_fasta_sequence[var_site_s-30:var_site_s],ref_fasta_sequence[var_site_s:var_site_e],ref_fasta_sequence[var_site_e:var_site_e+30])

    print(" "*10,''.join(reverse_complement(ref_fasta_sequence[var_site_e:var_site_e+30])),"MNN",''.join(reverse_complement(ref_fasta_sequence[var_site_s-30:var_site_s])),''.join(reverse_complement(ref_fasta_sequence[var_site_s-40:var_site_s-30])))

    #get forward primer 1
    for end_index_f_p1, char in enumerate(reversed(ref_fasta_sequence[0:30])):
        if char == "G" or char == "C":
            break

    FP1 = ref_fasta_sequence[0:30-end_index_f_p1]



    #get Reverse primer 2
    for end_index_R_p2, char in enumerate((ref_fasta_sequence[-30:])):
        if char == "G" or char == "C":
            break

    RP2 = ''.join(reverse_complement(ref_fasta_sequence[-30:]))

    #get reverse primer 1
    for end_index_RP1, char in enumerate(ref_fasta_sequence[var_site_s-(27):var_site_s]):
        if char == "G" or char == "C":
            break

    f_reverse_primer_1 = ref_fasta_sequence[var_site_s-(27-end_index_RP1):var_site_s]


    for start_index_RP1, char in enumerate(reversed(ref_fasta_sequence[var_site_e:var_site_e+30+end_index_RP1])):
        #print(char,start_index_RP1)
        if char == "G" or char == "C":
            break

    RP1 = ''.join([''.join(reverse_complement(ref_fasta_sequence[var_site_e:var_site_e+(27-end_index_RP1)])),"MNN",''.join(reverse_complement(ref_fasta_sequence[var_site_s-30-end_index_RP1+start_index_RP1:var_site_s]))])

    for start_index_FP2, char in enumerate(reversed(ref_fasta_sequence[var_site_s:var_site_s+30])):
        if char == "G" or char == "C":
            break

    for end_index_FP2, char in enumerate(reversed(ref_fasta_sequence[var_site_e:var_site_e+30])):
        if char == "G" or char == "C":
            break
    FP2 = ref_fasta_sequence[var_site_e:var_site_e+30-end_index_FP2]

    with open(Out_path, "w") as f:

        f.write("FP1"+" "+FP1+" "+str(len(FP1))+"\n")

        f.write("FP2"+" "+FP2+" "+str(len(FP2))+"\n")

        f.write("RP1"+" "+RP1+" "+str(len(RP1))+"\n")

        f.write("RP2"+" "+RP2+" "+str(len(RP2))+"\n")

        FP1_list_problem = []

        for index, char in enumerate(FP1[:-8]):
            FP1_list_problem.append(FP1[index:index+8])

        f.write("\n\n")

        f.write("FP1_list_problem:\n")
        for index,x in enumerate(FP1_list_problem):
            if ''.join(reverse_complement(x)) in FP1_list_problem:
                f.write(str(x)+" "+str(index) +" "+ str(FP1_list_problem.index(''.join(reverse_complement(x)))+"\n"))

        f.write("\n\n")

        FP2_list_problem = []

        for index, char in enumerate(FP2[:-8]):
            FP2_list_problem.append(FP2[index:index+8])

        f.write("RP2_list_problem:\n")
        for index,x in enumerate(FP2_list_problem):
            if ''.join(reverse_complement(x)) in FP2_list_problem:
                f.write(" "+str(x)+" "+str(index)+" "+str(FP2_list_problem.index(''.join(reverse_complement(x)))))

        f.write("\n")

        RP1_list_problem = []

        for index, char in enumerate(RP1[:-8]):
            RP1_list_problem.append(RP1[index:index+8])

        f.write("\n")


        f.write("RP1_list_problem:\n")
        for index,x in enumerate(RP1_list_problem):
            if ''.join(reverse_complement(x)) in RP1_list_problem:
                f.write(str(x)+" "+str(index)+" "+str(RP1_list_problem.index(''.join(reverse_complement(x)))))

        f.write("\n\n")

        RP2_list_problem = []

        for index, char in enumerate(RP2[:-8]):
            RP1_list_problem.append(RP2[index:index+8])

        f.write("RP2_list_problem:\n")
        for index,x in enumerate(RP2_list_problem):
            if ''.join(reverse_complement(x)) in RP2_list_problem:
                f.write(str(x)+" "+str(index)+" "+str(RP2_list_problem.index(''.join(reverse_complement(x)))))

    return Out_path
if __name__ == '__main__':

    Out_path = "1.txt"

    main("/Volumes/rdss_hhaim/LAB PROJECTS/Sam/FQ_VSearch/uploads/test_fasta.fa",Out_path)
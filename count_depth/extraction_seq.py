import os


def extraction_seq(seqname_path, seqtaxidaccession_path):
    seq2list_dict = {}
    with open(seqtaxidaccession_path, 'r') as files:
        for line in files:
            line_list = line.strip("\n").split("\t")
            seq2list_dict[line_list[0]] = line_list
    seq_line_list = []
    with open(seqname_path, 'r') as pf:
        for lin in pf:
            lin = lin.strip("\n")
            try:
                print(seq_line_list[lin])
                print(lin)
                seq_line_list.append(seq_line_list[lin])
            except:
                pass
    with open("./extraction_sequence/11.21/new_seq_taxid_accession.txt", 'w') as f:
        for seq_line in seq_line_list:
            f.write("\t".join(seq_line) + "\n")


if __name__ == '__main__':
    seqname_path = "/share/data5/hegh/project1/5.17/toolkit/Sample/count_depth/extraction_sequence/11.21/seq_list.txt"
    seqtaxidaccession_path = "/share/data5/hegh/project1/5.17/toolkit/Sample/count_depth/extraction_sequence/11.21/seq_taxid_accession.txt"
    extraction_seq(seqname_path, seqtaxidaccession_path)




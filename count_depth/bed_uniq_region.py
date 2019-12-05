import argparse
import base64
import hashlib
import json
import logging
import os
from urllib import parse

import requests

logger = logging.getLogger(__name__)


class UniqBed(object):
    def __init__(self, inputdir, inputbed, outdir, sign):
        self.inputdir = inputdir
        self.inputbed = inputbed
        self.outdir = outdir
        self.sign = sign

        # 读取汇总后bed文件中的数据
    def read_all_bed(self):
        all_bed_dict = {}
        with open(self.inputbed, 'r')as files:
            for line in files:
                line_list = line.strip("\n").split("\t")
                line_list[1] = int(line_list[1])
                line_list[2] = int(line_list[2])
                try:
                    locus_list = all_bed_dict[line_list[0]]
                    flag = True
                    for locus in locus_list:
                        if line_list[1] == locus[1]:
                            locus[1] = line_list[2]
                            flag = False
                            break
                        elif line_list[1] < locus[1]:
                            if line_list[1] < locus[0]:
                                locus[0] = line_list[1]
                            if line_list[2] > locus[1]:
                                locus[1] = line_list[2]
                            flag = False
                            break
                    if flag:
                        locus_list.append([line_list[1], line_list[2]])
                except:
                    all_bed_dict[line_list[0]] = [[line_list[1], line_list[2]]]
        print(len(all_bed_dict))
        return all_bed_dict

    def write_bed(self, all_bed_dict):
        bed_path = os.path.join(self.outdir, 'uniq.bed')
        with open(bed_path, 'w') as f:
            for key, value in all_bed_dict.items():
                for locus in value:
                    f.write(key + "\t" + str(locus[0]) + "\t" + str(locus[1]) + "\n")

    def run(self):
        all_bed_dict = self.read_all_bed()
        self.write_bed(all_bed_dict)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='从库文件中提取序列')
    parser.add_argument('-i', '--inputdir', help='存放输入文件list的路径',
                        default="/share/data2/xujm/docker/test_pathogenic_smk/analysis/qc_analysis/taxid_*/paired.taxid_*.q60.runid_gt3.bed")
    parser.add_argument('-f', '--inputbed', help='all_depth.bed',
                        default="/share/data2/xujm/docker/test_pathogenic_smk/staging/190813_MN00302_0308_A000H2W35Y/B1900731/mapping/bwa/GCF_000847245.1/B1900731.paired.GCF_000847245.1_uniq.bed")
    parser.add_argument('-o', '--outdir', help='输出文件的路径',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/count_depth/outfiles")
    parser.add_argument('-s', '--sign', help='是否需要生成中间文件seq_taxid_accession.txt 0代表不生成， 1代表生成', default=0)
    parser.add_argument('-d', '--debug', help='是否打开调试模式', default=False)
    args = parser.parse_args()
    logging_level = logging.DEBUG if args.debug else logging.INFO
    log_path = os.path.join(args.outdir, 'root.log')
    logging.basicConfig(level=logging_level, filename=log_path,
                        format='[%(asctime)s] p%(process)s {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s')
    # logger.debug(args.tsvfile, args.txtpath, args.outfna)
    cb = UniqBed(args.inputdir, args.inputbed, args.outdir, args.sign)
    cb.run()
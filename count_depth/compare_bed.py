import argparse
import base64
import hashlib
import json
import logging
import os
from urllib import parse

import requests

logger = logging.getLogger(__name__)


class CompareBed(object):
    def __init__(self, inputdir, inputbed, outdir, sign):
        self.inputdir = inputdir
        self.inputbed = inputbed
        self.outdir = outdir
        self.sign = sign

    # 通过输入路径获取bed文件的路径
    def get_bed_path(self):
        bed_path_list = []
        bed_path_files = os.path.join(self.outdir, 'bed_path.txt')
        code = os.system(f"ls {self.inputdir} > {bed_path_files}")
        if code == 0:
            with open(bed_path_files, 'r') as files:
                for line in files:
                    bed_path_list.append(line.strip("\n"))
        print(bed_path_list)
        return bed_path_list

    # 根据路径把所有的bed文件汇总到一起
    def collect_all_bed(self, bed_path_list):
        str_path = " ".join(bed_path_list)
        all_bed_path = os.path.join(self.outdir, 'all_bed_data.bed')
        code = os.system(f"cat {str_path} >> {all_bed_path}")
        if code == 0:
            pass

    # 读取汇总后bed文件中的数据
    def read_all_bed(self):
        all_bed_dict = {}
        all_bed_path = os.path.join(self.outdir, 'all_bed_data.bed')
        with open(all_bed_path, 'r')as files:
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

    # 读取输入文件中另一个bed文件中的数据
    def read_other_bed(self):
        other_bed_dict = {}
        with open(self.inputbed, 'r')as files:
            for line in files:
                line_list = line.strip("\n").split("\t")
                line_list[1] = int(line_list[1])
                line_list[2] = int(line_list[2])
                try:
                    locus_list = other_bed_dict[line_list[0]]
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
                    other_bed_dict[line_list[0]] = [[line_list[1], line_list[2]]]
        print(len(other_bed_dict))
        return other_bed_dict

    # 计算2个bed 文件的位点数，区间数
    def count_bed(self, all_bed_dict, other_bed_dict):
        count_bed_list = []
        common_bed_list = []
        for key, value in all_bed_dict.items():
            row_list = []
            row_list.append(key)
            length = 0
            for locus in value:
                length += locus[1] - locus[0]
            row_list.append(str(length))
            row_list.append(str(len(value)))
            try:
                other_row = other_bed_dict[key]
                olength = 0
                for olocus in other_row:
                    olength += olocus[1] - olocus[0]
                row_list.append(str(olength))
                row_list.append(str(len(other_row)))
                common_bed_list.append(row_list)
            except:
                row_list.append("-")
                row_list.append("-")
                count_bed_list.append(row_list)
        return count_bed_list, common_bed_list

    # 把计算后的结果写入到输出文件中
    def wrirte_bed_data(self, count_bed_list, common_bed_list):
        outfile = os.path.join(self.outdir, 'count_bed.tsv')
        with open(outfile, 'w') as f:
            for rowitem in common_bed_list:
                f.write("\t".join(rowitem) + "\n")
            for rowlist in count_bed_list:
                f.write("\t".join(rowlist) + "\n")
        print(f"{outfile} 文件写入完毕")
        logger.info(f"{outfile} 文件写入完毕")

    def run(self):
        bed_path_list = self.get_bed_path()
        self.collect_all_bed(bed_path_list)
        all_bed_dict = self.read_all_bed()
        other_bed_dict = self.read_other_bed()
        count_bed_list, common_bed_list = self.count_bed(all_bed_dict, other_bed_dict)
        self.wrirte_bed_data(count_bed_list, common_bed_list)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='从库文件中提取序列')
    parser.add_argument('-i', '--inputdir', help='存放输入文件list的路径',
                        default="/share/data2/xujm/docker/test_pathogenic_smk/analysis/qc_analysis/taxid_*/paired.taxid_*.q60.runid_gt3.bed")
    parser.add_argument('-f', '--inputbed', help='all_depth.bed',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/count_depth/all_depth.bed")
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
    cb = CompareBed(args.inputdir, args.inputbed, args.outdir, args.sign)
    cb.run()
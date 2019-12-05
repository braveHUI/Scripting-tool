import argparse
import logging
import os
import time
from logging.handlers import TimedRotatingFileHandler
import pysam
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class ComparisonBwa(object):
    def __init__(self, input_list, inputfna, inputdir, bwa, outdir):
        self.input_list = input_list
        self.inputfna = inputfna
        self.inputdir = inputdir
        self.bwa = bwa
        self.outdir = outdir

    def combine_path_files(self):
        bash_path = os.path.join(self.outdir, 'step1_combine.sh')
        fa_list = os.path.join(self.outdir, "fa_list.txt")
        code = os.system("echo 'cat {} {} > {} '>>{}".format(self.input_list, self.inputfna, fa_list, bash_path))
        if code == 0:
            os.system("bash {}".format(bash_path))
        return fa_list

    def combine_fa_files(self, fa_list):
        # print("11111111111")
        all_path = os.path.join(self.outdir, 'data')
        all_fa_file = os.path.join(all_path, 'all.fa')
        fa_path_list = []
        with open(fa_list, 'r') as files:
            for line in files:
                fa_path_list.append(line.strip("\n"))
        fa_path_str = " ".join(fa_path_list)
        step2_bash = os.path.join(self.outdir, 'step2_combine_fa.sh')
        shell = 'cat {} >> {}'.format(fa_path_str, all_fa_file)
        with open(step2_bash, 'w') as f:
            f.write(shell)
        # code = os.system("echo 'cat {} > {}' >>{}".format(fa_path_str, all_fa, step2_bash))
        os.system("bash {}".format(step2_bash))

    def bwa_build(self):
        all_fa = os.path.join(self.outdir, 'data')
        all_fa = os.path.join(all_fa, 'all.fa')
        step3_bash = os.path.join(self.outdir, 'step3_build_bwa.sh')
        # data_path = os.path.join(self.outdir, 'data')
        code = os.system("echo 'bwa index {}' > {}".format(all_fa, step3_bash))
        if code == 0:
            os.system("bash {}".format(step3_bash))

    def get_reads_path(self):
        files = os.listdir(self.inputdir)
        dir_list = [os.path.join(self.inputdir, dis) for dis in files if
                    os.path.isdir(os.path.join(self.inputdir, dis)) and dis.startswith("19")]
        bk_dir_list = []
        for path in dir_list:
            bk_files = os.listdir(path)
            for fil in bk_files:
                if os.path.isdir(os.path.join(path, fil)):
                    if fil.startswith("BK") or fil.startswith("ZK"):
                        bk_dir_list.append(os.path.join(path, fil))
        files_path_list = []
        name_item_list = []
        for file_pa in bk_dir_list:
            filepath = os.path.join(file_pa, "clean")
            files_list = os.listdir(filepath)
            name_list = [name for name in files_list if name.endswith("filHuman.fastq") and not name.startswith("single")]
            for name in name_list:
                name_item_list.append(name)
            if len(name_list) == 2:
                strna = str(name_list[0].split(".")[0]).split("_")[-1]
                if strna == "R1":
                    files_path = [os.path.join(filepath, name) for name in name_list]
                    files_path_list.append(files_path)
                else:
                    files_path = [os.path.join(filepath, name_list[1]),  os.path.join(filepath, name_list[0])]
                    files_path_list.append(files_path)
            else:
                print(name_list)
        print(len(name_item_list))
        name_item_list = {}.fromkeys(name_item_list).keys()
        print(len(name_item_list))
        return files_path_list

    def get_sam_files(self, files_path_list):
        sam_path = []
        sam_path_dir = os.path.join(self.outdir, 'sam_dirs')
        step4_bash = os.path.join(self.outdir, 'step4_bwa_alignment.sh')
        # all_fa_path = os.path.join(self.outdir, 'data/all.fa')
        all_fa_path = "/share/data5/hegh/project1/5.17/toolkit/Sample/bwa/most_files_bwa/data/all.fa"
        print(len(files_path_list))
        with open(step4_bash, 'w') as f:
            for path_liat in files_path_list:
                gz = str(path_liat[0].split("/")[-1]).split(".")[0].split("_")[0]
                gz = str(gz) + ".sam"
                outfile = os.path.join(sam_path_dir, gz)
                sam_path.append(outfile)
                strshell = "/share/data2/hegh/miniconda3/envs/ncbi/bin/bwa mem -t 4 {} {} {} >{}".format(all_fa_path, path_liat[0], path_liat[1], outfile)
                f.write(strshell + "\n")
        # code = os.system("perl /share/data2/xujm/bin/submit-jobs.sh.SGE.pl  -i {} -n 80 -c 4 -m 5G -r".format(step4_bash))
        # code = os.system("bash {}".format(step4_bash))
        # if code == 0:
        return sam_path

    def samtools_sort_file(self, sam_path):
        sam_path_dir = os.path.join(self.outdir, 'sam_dirs')
        # sam_files = os.listdir(sam_path_dir)
        # sam_files = [name for name in sam_files if name.endswith(".sam")]
        sort_sam_list = []
        step5_bash = os.path.join(self.outdir, 'step5_samtools_sort.sh')
        with open(step5_bash, 'w') as f:
            for samp in sam_path:
                sort_name = samp.split("/")[-1].split(".")[0]
                sampath = os.path.join(sam_path_dir, samp)
                sort_sam_path = os.path.join(sam_path_dir, sort_name + "_sort.bam")
                sort_sam_list.append(sort_sam_path)
                strshell = "/share/data2/hegh/miniconda3/envs/ncbi/bin/samtools sort {} > {}".format(sampath, sort_sam_path)
                f.write(strshell + "\n")
        # code = os.system("bash {}".format(step5_bash))
        # if code == 0:
        return sort_sam_list

    def samtools_index_files(self, sort_sam_list):
        step6_bash = os.path.join(self.outdir, 'step6_samtools_index.sh')
        with open(step6_bash, 'w') as f:
            for sortpa in sort_sam_list:
                strshell = "/share/data2/hegh/miniconda3/envs/ncbi/bin/samtools index {}".format(sortpa)
                f.write(strshell + "\n")
        # code = os.system("bash {}".format(step6_bash))
        # if code == 0:
        return sort_sam_list

    def samtools_depth_fles(self, sort_sam_list):
        step7_bash = os.path.join(self.outdir, 'step7_samtools_depth.sh')
        depth_path = os.path.join(self.outdir, 'txt_dirs')
        with open(step7_bash, 'w') as f:
            for samname in sort_sam_list:
                name = samname.split("/")[-1].split("_")[0]
                depth_name = name + "_depth_Q60.txt"
                depth_file = os.path.join(depth_path, depth_name)
                strshell = "/share/data2/hegh/miniconda3/envs/ncbi/bin/samtools depth -Q 60 {} > {}".format(samname, depth_file)
                f.write(strshell + "\n")
        # code = os.system("bash {}".format(step7_bash))
        # if code == 0:
        sam_path_dir = os.path.join(self.outdir, 'sam_dirs')
        sam_files = os.listdir(sam_path_dir)
        sam_path_files = [os.path.join(sam_path_dir, name) for name in sam_files if name.endswith(".sam")]
        step8_bash = os.path.join(self.outdir, 'step8_detele_sam.sh')
        with open(step8_bash, 'w')as f:
            for path in sam_path_files:
                strshell= "rm {}".format(path)
                f.write(strshell + "\n")

    def run(self):
        # fa_list = self.combine_path_files()
        fa_list = "/share/data5/hegh/project1/5.17/toolkit/Sample/bwa/most_files_bwa/fa_list.txt"
        # self.combine_fa_files(fa_list)
        # self.bwa_build()
        files_path_list = self.get_reads_path()
        files_path_list = [files_path_list[i] for i in range(0, 3)]
        sam_path = self.get_sam_files(files_path_list)
        sort_sam_list = self.samtools_sort_file(sam_path)
        sort_sam_list = self.samtools_index_files(sort_sam_list)
        self.samtools_depth_fles(sort_sam_list)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Scan Illumina Runs and send status to Promegene BMS')
    parser.add_argument('-i', '--input_list', help='存放输入文件map的路径',
                        default="/share/data5/zhangh/database/Pathogen/pan_genome_reference.list")
    parser.add_argument('-f', '--inputfna', help='存放输入文件map的路径',
                        default="/share/data5/zhangh/database/Pathogen/uniq_refseq_taxid_prokka_ffn.list")
    parser.add_argument('-p', '--inputdir', help='存放下机数据的路径',
                        default="/share/data7/pmid/2019")
    parser.add_argument('-n', '--bwa', help='bwa建库路径',
                        default="/share/data6/PMiD/references/bwa_db/seq5_7/all.fna")
    parser.add_argument('-o', '--outdir', help='存放输出文件的目录',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/bwa/most_files_bwa/test")
    parser.add_argument('-d', '--debug', help='是否打开调试模式', default=False)
    args = parser.parse_args()
    logging_level = logging.DEBUG if args.debug else logging.INFO
    log_path = args.outfile + ".log"
    handler = TimedRotatingFileHandler(log_path,
                                       when="midnight",
                                       interval=1,
                                       backupCount=5)
    logging.basicConfig(level=logging_level,
                        format='[%(asctime)s] p%(process)s {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
                        handlers=[handler])
    start = time.process_time()
    cb = ComparisonBwa(args.input_list, args.inputfna, args.inputdir, args.bwa, args.outdir)
    cb.run()
    end = time.process_time()
    logger.debug("main()开始时间是{}， 结束时间是{} ，总运行时间是{}".format(start, end, end-start))

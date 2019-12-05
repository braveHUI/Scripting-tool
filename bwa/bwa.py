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
    def __init__(self, inputfile, inputdir, outfile, bwa, outdir):
        self.inputfile = inputfile
        self.outfile = outfile
        self.path = inputdir
        self.bwa = bwa
        self.outdir = outdir

    def read_excel_data(self):
        datas = pd.read_excel(self.inputfile)
        data_array = np.array(datas)
        data_list = data_array.tolist()
        path_list = []
        for line in data_list:
            item = {}
            list_p = line[-1].split("_")
            if len(list_p) > 2:
                list_r = list_p[1]
                if "MN00302" == list_r:
                    item["sequencer"] = 3
                    item["runid"] = line[-1]
                    item["gz"] = line[0]
                    path_list.append(item)
                    logger.debug("增加%s成功" % item)
                elif "M04840" == list_r:
                    item["sequencer"] = 2
                    item["runid"] = line[-1]
                    item["gz"] = line[0]
                    path_list.append(item)
        # data_dict = {line[0]: os.path.join(self.inputdir, line[-1]) for line in data_list}
        print(path_list)
        return path_list

    # 寻找存放fastq文件的最新的Alignment文件夹
    def find_alignment_dir(self, runid):
        logger.debug("进入find_alignment_dir方法")
        air_dirs = []
        num = 0
        try:
            ali_path = os.path.join(self.path, runid)
            dire_al = os.listdir(ali_path)
            for ali in dire_al:
                if ali.startswith("Alignment"):
                    air_dirs.append(ali)
            if len(air_dirs):
                for ae in air_dirs:
                    ali_num = int(ae.split("_")[1])
                    if ali_num > num:
                        num = ali_num
        except Exception as e:
            logger.error(e)
        finally:
            ali_dirs_name = "Alignment" + "_" + str(num)
            logger.debug("最新的文件夹是%s" % ali_dirs_name)
            logger.debug("退出find_alignment_dir方法")
            return ali_dirs_name

    # 在文件夹中找出fastq文件
    def get_fastq_file(self, diritem):
        logger.debug("进入get_fastq_file方法")
        fastq_path_list = []
        if diritem["sequencer"] == 3:
            fastqtt_path = self.find_alignment_dir(diritem["runid"])
            fastq_path = os.path.join(self.path, diritem["runid"], fastqtt_path)
            if os.path.exists(fastq_path):
                try:
                    temp = os.listdir(fastq_path)[0]
                    fastq_path = os.path.join(fastq_path, temp, 'Fastq')
                    fastq_list = os.listdir(fastq_path)
                    if len(fastq_list):
                        fastq_list = [faq for faq in fastq_list if faq.endswith("fastq.gz") and diritem["gz"] in faq]
                        fastq_list = sorted(fastq_list)
                        for fastq in fastq_list:
                            fastq = os.path.join(fastq_path, fastq)
                            fastq_path_list.append(fastq)
                except Exception as e:
                    logger.error("%s查找失败,错误信息：%s" % (diritem["runid"], e))
            logger.debug("退出get_fastq_file方法")
        elif diritem["sequencer"] == 2:
            fastq_path = os.path.join(self.path, diritem["runid"], 'Data/Intensities/BaseCalls')
            if os.path.exists(fastq_path):
                try:
                    fastq_list = os.listdir(fastq_path)
                    if len(fastq_list):
                        fastq_list = [faq for faq in fastq_list if faq.endswith("fastq.gz") and diritem["gz"] in faq]
                        fastq_list = sorted(fastq_list)
                        fastq_list = sorted(fastq_list)
                        for fastq in fastq_list:
                            fastq = os.path.join(fastq_path, fastq)
                            fastq_path_list.append(fastq)
                except Exception as e:
                    logger.error("%s查找失败,错误信息：%s" % (diritem["runid"], e))
            logger.debug("退出get_fastq_file方法")
        return fastq_path_list

    def get_sam_files(self, path_list):
        sam_path = []
        for path in path_list:
            fastq_path_list = self.get_fastq_file(path)
            if len(fastq_path_list) == 2:
                # print(fastq_path_list)
                gz = str(fastq_path_list[0].split("/")[-1]).split(".")[0]
                gz = str(gz) + "_aln-se.sam"
                outfile = os.path.join(self.outdir, gz)
                # os.system("bwa mem {} {} | gzip -3 > {}".format(self.bwa, fastq, outfile))
                logger.info("bwa mem {} {} {} > {}".format(self.bwa, fastq_path_list[0], fastq_path_list[1], outfile))
                # code = os.system("bwa mem {} {} {} > {}".format(self.bwa, fastq_path_list[0], fastq_path_list[1], outfile))
                # if code == 0:
                #     sam_path.append(outfile)
        return sam_path

    def get_bam_files(self, sam_path):
        bam_files = []
        for path in sam_path:
            sam = path.split("/")[-1].split(".")[0]
            sam += ".bam"
            bam_path = os.path.join(self.outdir, sam)
            code = os.system("samtools view -b -S {} > {} ".format(path, bam_path))
            if code == 0:
                bam_files.append(bam_path)
        return bam_files

    def parse_bam_file(self):
        files_list = os.listdir(self.outdir)
        bam_files = [os.path.join(self.outdir, name) for name in files_list if name.endswith(".bam")]
        for path in bam_files:
            name = str(str(path.split("/")[-1]).split(".")[0])
            prefix_path = os.path.join(self.outdir, name + "Prefix.bam")
            print(prefix_path)
            print(path)
            code = os.system("samtools sort  {} {}".format(path, prefix_path))
            if code == 0:
                code1 = os.system("samtools index {}".format(path))
                if code1 == 0:
                    seq_5 = 0
                    seq_7 = 0
                    samfile = pysam.AlignmentFile(path, 'rb')
                    for read in samfile.fetch():
                        if read.is_proper_pair:
                            if read.reference_name == "seq_7" and read.mapq >= 60:
                                seq_7 += 1
                            elif read.reference_name == "seq_5" and read.mapq >= 60:
                                seq_5 += 1
                    parse_path = name + "parse.txt"
                    txt_path = os.path.join(self.outdir, parse_path)
                    with open(txt_path, 'w') as f:
                        f.write("seq_5" + "\t" + "seq_7" + "\n")
                        f.write(str(seq_5) + "\t" + str(seq_7) + "\n")

    def parse_bam_file1(self):
        files_list = os.listdir(self.outdir)
        bam_files = [os.path.join(self.outdir, name) for name in files_list if name.endswith(".bam")]
        for path in bam_files:
            name = str(str(path.split("/")[-1]).split(".")[0])
            seq_5 = 0
            seq_7 = 0
            samfile = pysam.AlignmentFile(path, 'rb')
            for read in samfile.fetch(until_eof=True):
                if read.is_proper_pair:
                    if read.reference_name == "seq_7" and read.mapq >= 60:
                        seq_7 += 1
                    elif read.reference_name == "seq_5" and read.mapq >= 60:
                        seq_5 += 1
            parse_path = name + "parse.txt"
            txt_path = os.path.join(self.outdir, parse_path)
            with open(txt_path, 'w') as f:
                f.write("seq_5" + "\t" + "seq_7" + "\n")
                f.write(str(seq_5) + "\t" + str(seq_7) + "\n")



    def count_bam_reads(self, bam_files):
        for path in bam_files:
            name = path.split("/")[-1].split(".")[0]
            idxstat_path = name + "_idxstat.txt"
            txt_path = os.path.join(self.outdir, idxstat_path)

            code1 = os.system("samtools index {}".format(path))
            if code1 == 0:
                code2 = os.system("samtools idxstats {} > {}".format(path, txt_path))
                if code2 == 0:
                    items = []
                    with open(txt_path, 'r') as files:
                        for line in files:
                            line_list = line.strip("\n").split("\t")
                            items.append(line_list[2])
                    items.pop()
                    mapple_name = name + "mapped.txt"
                    mapped_path = os.path.join(self.outdir, mapple_name)
                    with open(mapped_path, 'w') as f:
                        f.write("seq_5" + "\t" + "seq_7" + "\n")
                        f.write("\t".join(items) + "\n")

    def count_bamreads(self):
        files_list = os.listdir(self.outdir)
        bam_files = [os.path.join(self.outdir, name) for name in files_list if name.endswith(".bam")]
        print(bam_files)
        for path in bam_files:
            name = str(path.split("/")[-1].split(".")[0])
            idxstat_path = name + "_idxstat.txt"
            txt_path = os.path.join(self.outdir, idxstat_path)
            code1 = os.system("samtools index {}".format(path))
            if code1 == 0:
                code2 = os.system("samtools idxstats {} > {}".format(path, txt_path))
                if code2 == 0:
                    items = []
                    with open(txt_path, 'r') as files:
                        for line in files:
                            line_list = line.strip("\n").split("\t")
                            items.append(line_list[2])
                    items.pop()
                    mapple_name = name + "mapped.txt"
                    mapped_path = os.path.join(self.outdir, mapple_name)
                    with open(mapped_path, 'w') as f:
                        f.write("seq_5" + "\t" + "seq_7" + "\n")
                        f.write("\t".join(items) + "\n")


    def run(self):
        path_list = self.read_excel_data()
        # print(path_list)
        sam_files = self.get_sam_files(path_list)
        # bam_files = self.get_bam_files(sam_files)
        # # bam_files = ["/share/data5/hegh/project1/5.17/toolkit/Sample/bwa/test/test.bam"]
        # self.count_bam_reads(bam_files)
        # self.count_bamreads()
        # self.parse_bam_file1()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Scan Illumina Runs and send status to Promegene BMS')
    parser.add_argument('-i', '--inputfile', help='存放输入文件map的路径',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/bwa/etiology.xlsx")
    parser.add_argument('-p', '--inputdir', help='存放下机数据的路径',
                        default="/share/data4/illumina/")
    parser.add_argument('-n', '--bwa', help='bwa建库路径',
                        default="/share/data6/PMiD/references/bwa_db/seq5_7/all.fna")
    parser.add_argument('-', '--outdir', help='存放输出文件的目录',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/bwa/bam_mapped_files")
    parser.add_argument('-o', '--outfile', help='存放输出文件的路径',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/bwa/bam_mapped_files")
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
    cb = ComparisonBwa(args.inputfile, args.inputdir, args.outfile, args.bwa, args.outdir)
    cb.run()
    end = time.process_time()
    logger.debug("main()开始时间是{}， 结束时间是{} ，总运行时间是{}".format(start, end, end-start))

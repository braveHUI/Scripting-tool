import argparse
import base64
import hashlib
import json
import logging
import os
from urllib import parse
import pysam
import requests

logger = logging.getLogger(__name__)


class ExtractAlignmentSequence(object):
    def __init__(self, inputrefseqfile, inputdir, outdir, sign):
        self.inputrefseqfile = inputrefseqfile
        self.inputdir = inputdir
        self.outdir = outdir
        self.sign = int(sign)

    def read_seq_data(self):
        with open(self.inputrefseqfile, 'r') as files:
            seq_list = []
            for line in files:
                line_list = line.strip("\n").split("\t")
                if line_list[1] == "1031542":
                    seq_list.append(line_list[0])
        return seq_list

    def get_sam_path(self):
        sam_path = os.listdir(self.inputdir)
        sam_file_path = [os.path.join(self.inputdir, name) for name in sam_path if name.endswith(".sam")]
        outbam_path = os.path.join(self.outdir, '1031542.bam')
        samfile1 = pysam.AlignmentFile("/share/data5/hegh/project1/5.17/toolkit/Sample/bwa/most_files_bwa/sam_dirs/ZK190814_R1.sam", "rb")
        pairedreads = pysam.AlignmentFile(outbam_path, "ab", template=samfile1)
        for path in sam_file_path:
            samfile = pysam.AlignmentFile(path, "rb")
            for read in samfile.fetch():
                try:
                    value = seq_list[read.reference_name]
                    pairedreads.write(read)
                except:
                    pass
            samfile.close()
        pairedreads.close()






    def run(self):
        seq_list = self.read_seq_data()
        self.get_sam_path()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='从库文件中提取序列')
    parser.add_argument('-f', '--inputrefseqfile', help='存放输入文件的路径seq_taxid_accession.txt',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/count_depth/outfiles/seq_taxid_accession.txt")
    parser.add_argument('-i', '--inputdir', help='存放输入文件的目录',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/bwa/most_files_bwa/sam_dirs")
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
    eas = ExtractAlignmentSequence(args.inputrefseqfile, args.inputdir, args.outdir, args.sign)
    eas.run()
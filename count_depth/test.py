import argparse
import logging
import os
import time
from logging.handlers import TimedRotatingFileHandler

logger = logging.getLogger(__name__)


class GetSpeciesTaxid(object):
    def __init__(self, inputfile, inputdir):
        self.inputfile = inputfile
        self.inputdir = inputdir

    def write_data(self):
        speciestaxid2path = {}
        with open(self.inputfile, 'r') as files:
            for line in files:
                if not line.startswith("#"):
                    line_list = line.strip("\n").split("\t")
                    try:
                        path_list = speciestaxid2path[line_list[6]]
                        ftp_path = line_list[-3].split("//")[-1].split("/")[-4:]
                        strftppath = "/".join(ftp_path) + "/" + ftp_path[-1] + "_genomic.fna.gz"
                        gcf_path = os.path.join(self.inputdir, strftppath)
                        path_list.append(gcf_path)
                    except:
                        ftp_path = line_list[-3].split("//")[-1].split("/")[-4:]
                        strftppath = "/".join(ftp_path) + "/" + ftp_path[-1] + "_genomic.fna.gz"
                        gcf_path = os.path.join(self.inputdir, strftppath)
                        speciestaxid2path[line_list[6]] = [gcf_path]
        return speciestaxid2path

    def get_species_path(self, speciestaxid2path):
        for key, value in speciestaxid2path.items():
            if len(value) == 3:
                print("*" * 200)
                print(key)
                print(value)

    def run(self):
        speciestaxid2path = self.write_data()
        self.get_species_path(speciestaxid2path)





if __name__ == '__main__':
    '''
           1.根据assembly_summary_refseq.txt文件中的assembly_accession的ftp路径结合下载的GCF的路径去获取report的路径，
           2. 根据获取的report的路径去读取其中RefSeq-Accn和Sequence-Length这两列的数据，
            保存到对应的assembly_summary_refseq这一列后，生成asembly_sumary_refseq_update.txt文件， 
            同时生成包含RefSeq-Accn，taxid， species_taxid， assembly_accession这四列的reseq_taxid_sepcie_accession.txt文件
       '''
    parser = argparse.ArgumentParser(description='Scan Illumina Runs and send status to Promegene BMS')
    parser.add_argument('-i', '--inputfile', help='存放输入文件的路径assembly_summary_refseq.txt',
                        default="/share/data6/PMiD/library/refseq/ascp/assembly_summary_refseq.txt")
    parser.add_argument('-p', '--inputdir', help='存放输入文件report的目录',
                        default="/share/data6/PMiD/library/refseq/ascp/GCF/")
    parser.add_argument('-o', '--outdir', help='存放输出文件的路径',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/assembly/outfile")
    parser.add_argument('-y', '--yest', help='是否需要重新下载assembly_summary.txt文件',
                        default="0")
    parser.add_argument('-d', '--debug', help='是否打开调试模式', default=False)
    args = parser.parse_args()
    logging_level = logging.DEBUG if args.debug else logging.INFO
    log_path = os.path.join(args.outdir, "requests.log")
    handler = TimedRotatingFileHandler(log_path,
                                       when="midnight",
                                       interval=1,
                                       backupCount=5)
    logging.basicConfig(level=logging_level,
                        format='[%(asctime)s] p%(process)s {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
                        handlers=[handler])
    start = time.process_time()
    gst = GetSpeciesTaxid(args.inputfile, args.inputdir)
    gst.run()
    end = time.process_time()
    logger.debug("main()开始时间是{}， 结束时间是{} ，总运行时间是{}".format(start, end, end-start))
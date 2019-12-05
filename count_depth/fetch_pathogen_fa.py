import argparse
import base64
import hashlib
import json
import logging
import os
from urllib import parse

import requests

logger = logging.getLogger(__name__)


class ExtractionSequence(object):
    def __init__(self, input_list, inputfna, update_fa, outdir, sign):
        self.input_list = input_list
        self.inputfna = inputfna
        self.update_fa = update_fa
        self.outdir = outdir
        self.sign = int(sign)

    # 根据两个输入文件去获取fa的路径
    def get_fa_path(self):
        panfa_path_list = []
        prokka_path_list = []
        with open(self.input_list, 'r') as files:
            for line in files:
                line_path = line.strip("\n")
                panfa_path_list.append(line_path)
        with open(self.inputfna, 'r') as pf:
            for line in pf:
                line_path = line.strip("\n")
                prokka_path_list.append(line_path)
        return panfa_path_list, prokka_path_list

    # 根据fa路径去获取taxid 和文件中的序列号
    def get_taxid_seq(self, panfa_path_list, prokka_path_list):
        seq2taxid = {}
        for path in panfa_path_list:
            taxid = path.split("/")[-3].split("_")[-1]
            seq_file= os.path.join(self.outdir, 'seq.txt')
            code = os.system("grep '>' {} > {}".format(path, seq_file))
            if code == 0:
                with open(seq_file, 'r') as files:
                    for line in files:
                        seq_name = line.strip("\n").split(" ")[0].strip(">")
                        seq2taxid[seq_name] = taxid
        for patha in prokka_path_list:
            taxid = patha.split("/")[-4].split("_")[-1]
            seq_file = os.path.join(self.outdir, 'prokka_seq.txt')
            code = os.system("grep '>' {} > {}".format(patha, seq_file))
            if code == 0:
                with open(seq_file, 'r') as files:
                    for line in files:
                        seq_name = line.strip("\n").split(" ")[0].strip(">")
                        seq2taxid[seq_name] = taxid
        return seq2taxid

    # 生成一个包含RefSeq-Accn，taxid， species_taxid， assembly_accession四列的文件
    def write_seq2taxid(self, seq2taxid):
        out_file = os.path.join(self.outdir, 'seq_taxid_accession.txt')
        name_list = ["seq_name", "taxid", "species_taxid", "assembly_accession"]
        with open(out_file, 'w') as f:
            f.write("\t".join(name_list) + "\n")
            for key, value in seq2taxid.items():
                f.write(key + "\t" + value + "\t" + value + "\t" + "-" + "\n")

    # 读取生成的seq_taxid_accession.txt过渡文件， 放回一个字典key 为taxid， value为一个包含seq的列表
    def read_seq2taxid(self):
        taxid2seq = {}
        out_file = os.path.join(self.outdir, 'seq_taxid_accession.txt')
        with open(out_file, 'r') as files:
            for line in files:
                line_list = line.strip("\n").split("\t")
                try:
                    seq_list = taxid2seq[line_list[1]]
                    seq_list.append(line_list[0])
                except:
                    taxid2seq[line_list[1]] = [line_list[0]]
        return taxid2seq

    # 发送接口api 获取放回的数据，获取其中的taxid
    def request_api(self):
        api_user = os.environ.get('API_USER')
        api_key = os.environ.get('API_KEY')
        use_url = os.environ.get('API_URL')
        # use_url = "http://dev.ttjbz.com:8888"
        run_url = use_url + "/api/pmid/wikipmid/list/"

        def forSyncAPI(request_data, data_sign, api_user):
            json = {
                'request_data': parse.quote(request_data),
                'data_sign': parse.quote(data_sign),
                'api_user': api_user
            }
            return json

        request_data = json.dumps({})
        md5 = hashlib.md5((request_data + api_key).encode("utf-8"))
        # base64计算结果
        data_sign = base64.b64encode(md5.hexdigest().encode("utf-8"))
        re_data = forSyncAPI(request_data, data_sign, api_user)
        response = requests.post(url=run_url, data=re_data, headers={})
        status_code = response.status_code
        taxid_list = []
        if status_code == 200:
            re_dict = json.loads(response.content.decode())
            print(re_dict)
            logger.info(re_dict)
            re_dict_data = re_dict["data"]
            for dict in re_dict_data:
                taxid_list.append(str(dict["tax_id"]))
        return taxid_list

    def taxid_get_seq(self, taxid_list, taxid2seq):
        seq_list = []
        err_list = []
        logger.info(len(taxid_list))
        for taxid in taxid_list:
            try:
                seq_item = taxid2seq[taxid]
                seq_list.extend(seq_item)
            except:
                err_list.append(taxid)
                logger.error(f"字典中不存在该{taxid}")
        seq_files = os.path.join(self.outdir, 'seq_name.lst')
        logger.info(len(err_list))
        with open(seq_files, 'w') as f:
            for seq in seq_list:
                f.write(seq + "\n")
        return seq_files

    def extract_sequence(self, seq_files):
        out_file = os.path.join(self.outdir, 'out.fa')
        os.system(f"seqtk subseq {self.update_fa} {seq_files} > {out_file}")

    def run(self):
        if self.sign:
            panfa_path_list, prokka_path_list = self.get_fa_path()
            seq2taxid = self.get_taxid_seq(panfa_path_list, prokka_path_list)
            self.write_seq2taxid(seq2taxid)
        taxid2seq = self.read_seq2taxid()
        taxid_list = self.request_api()
        seq_files =self.taxid_get_seq(taxid_list, taxid2seq)
        self.extract_sequence(seq_files)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='从库文件中提取序列')
    parser.add_argument('-i', '--input_list', help='存放输入文件list的路径',
                        default="/share/data5/zhangh/database/Pathogen/pan_genome_reference.list")
    parser.add_argument('-f', '--inputfna', help='存放输入文件ffn的路径',
                        default="/share/data5/zhangh/database/Pathogen/uniq_refseq_taxid_prokka_ffn.list")
    parser.add_argument('-u', '--update_fa', help='存放输入文件update_all.fa的路径',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/count_depth/update_all.fa")
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
    es = ExtractionSequence(args.input_list, args.inputfna, args.update_fa, args.outdir, args.sign)
    es.run()
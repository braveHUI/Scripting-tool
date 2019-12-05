import argparse
import base64
import copy
import hashlib
import json
import logging
import os
from urllib import parse
import time
import numpy as np
import requests

logger = logging.getLogger(__name__)


# 捕获异常的方法
def catch_exception(func):
    def wrapper(self, *args, **kwargs):
        try:
            return func(self, *args, **kwargs)
        except Exception as e:
            logger.exception(e)

    return wrapper


# 写日志类
class write_logger(object):
    def __init__(self):
        pass

    def __call__(self, func):
        def _call(*args, **kwargs):
            logger.debug('%s method is running' % func.__name__)
            return func(*args, **kwargs)

        return _call


class CountGenomeCoverage(object):
    def __init__(self, inputrefseqfile, statdir, sname, tname, name_dmp, outdir):
        self.inputrefseqfile = inputrefseqfile
        self.statdir = statdir
        self.sname = sname
        self.tname = tname
        self.name_dmp = name_dmp
        self.outdir = outdir
        self.api_user = os.environ.get('API_USER')
        self.api_key = os.environ.get('API_KEY')
        # self.use_url = "http://dev.ttjbz.com:8888"
        self.use_url = os.environ.get('API_URL')
        if self.use_url.endswith("/"):
            self.run_post = self.use_url + "api/pmid/species/add/"
        else:
            self.run_post = self.use_url + "/api/pmid/species/add/"

    # 读取reseq_taxid_sepcie_accession.txt文件中的数据生成一个字典refseq2taxid，键值对是：RefSeq-Accn对应的该列数据
    @write_logger()
    def read_reseq_txt(self):
        start_time = time.process_time()
        logger.debug("进入read_reseq_txt方法")
        refseq2speciestaxid = {}
        print(self.inputrefseqfile)
        with open(self.inputrefseqfile, 'r') as pf:
            for line in pf:
                line_list = line.strip("\n").split("\t")
                refseq2speciestaxid[line_list[0]] = line_list
        logger.debug("退出read_reseq_txt方法")
        end_time = time.process_time()
        print(f"read_reseq_txt方法开始运行的时间：{start_time},结束运行的时间: {end_time}, 一共运行的时间： {end_time - start_time}")
        return refseq2speciestaxid

    # 读取name_dmp文件中species_taxid对应的name返回一个字典
    def read_name_dmp(self):
        specise_taxid2name = {}
        with open(self.name_dmp, 'r') as files:
            for line in files:
                line_list = line.strip("\t|\n").split("\t|\t")
                if line_list[3] == "scientific name":
                    specise_taxid2name[line_list[0]] = line_list[1]

        return specise_taxid2name

    # 获取stat文件的路径和包含uniq_read总长度的文件
    def get_stat_path(self):
        start_time = time.process_time()
        stat_name_path_dict = {}
        runid_list = {}
        file_path = os.path.join(self.outdir, "stat.path")
        shell1 = "find {} -name '{}' > {}".format(self.statdir, self.sname, file_path)
        code = os.system(shell1)
        if code == 0:
            files = open(file_path, 'r')
            for line in files:
                line = line.strip("\n")
                file_name = line.split("/")[-5]
                runid = line.split("/")[-6]
                stat_name_path_dict[file_name] = [line]
                runid_list[file_name] = runid
        file1_path = os.path.join(self.outdir, "single_stat.path")
        shell2 = "find {} -name '{}' > {}".format(self.statdir, self.tname, file1_path)
        code = os.system(shell2)
        if code == 0:
            with open(file1_path, 'r') as pf:
                for sline in pf:
                    sline = sline.strip("\n")
                    length_name = sline.split("/")[-3]
                    value = stat_name_path_dict[length_name]
                    value.append(sline)
        end_time = time.process_time()
        print(f"get_stat_path方法开始运行的时间：{start_time},结束运行的时间: {end_time}, 一共运行的时间： {end_time-start_time}")
        return stat_name_path_dict, runid_list

    # 根据包含uniq_read总长度文件的路径获取uniq_read总长度
    def get_length_stat(self, path):
        all_length = 0
        with open(path, 'r')as files:
            for line in files:
                line_list = line.strip("\n").split("\t")
                if line_list[0] == "rmHuman":
                    all_length = int(line_list[1])
        return all_length

    # 读取 single.stat.tsv中的数据
    def read_single_data(self, spath):
        logger.debug("进入read_single_data方法")
        single_list = []
        with open(spath, 'r') as files:
            for line in files:
                line_list = line.strip("\n").split("\t")
                single_list.append(line_list)
        del single_list[0]
        logger.debug("退出read_single_data方法")
        return single_list

    # 根据single_list中每列的Seq_ID从refseq2taxid中查找对应的taxid，RefSeq_assembly_accession并插入和替换之后计算相同基因组的覆盖率
    def analysis_sepcies_taxid(self, single_list, refseq2speciestaxid):
        start_time = time.process_time()
        logger.debug("进入analysis_sepcies_taxid方法")
        sepcies_count_items = {}  # 用于保存最后输出的数据
        logger.debug("进入analysis_sepcies_taxid方法")
        for line in single_list:
            try:
                line_list = refseq2speciestaxid[line[0]]
                del line[0]
                line.insert(0, line_list[1])
                line.insert(0, line_list[2])
            except Exception as e:
                logger.error("{}找不到对应的taxid和species_taxid的数据".format(line[0]))
                logger.warning(e)
        # 计算相同基因组的覆盖率
        logger.debug(single_list)
        for ier in single_list:
            try:
                single_list = sepcies_count_items[ier[0]]
                single_list[2] = str(int(single_list[2]) + int(ier[2]))
                single_list[3] = str(int(single_list[3]) + int(ier[3]))
                single_list[4] = str(int(single_list[3]) / int(single_list[2]))
                single_list[5] = str(int(single_list[5]) + int(ier[5]))
                single_list[6] = str(int(single_list[5]) / int(single_list[2]))
                single_list[7] = str(int(single_list[7]) + int(ier[7]))
                single_list[8] = str(int(single_list[8]) + int(ier[8]))
            except:
                sepcies_count_items[ier[0]] = ier
        logger.debug(sepcies_count_items)
        sepcies_count_list = [line for line in sepcies_count_items.values()]
        if len(sepcies_count_list):
            sepcies_count_array = np.array(sepcies_count_list)
            sepcies_unique_array = np.unique(sepcies_count_array, axis=0)
            sepcies_sort_unique_array = sepcies_unique_array[
                np.lexsort([-1 * sepcies_unique_array[:, 7].astype(int), -1 * sepcies_unique_array[:, 8].astype(int)])]
        else:
            sepcies_sort_unique_array = sepcies_count_list

        logger.debug(sepcies_sort_unique_array)
        logger.debug("退出analysis_sepcies_taxid方法")
        # print(len(sepcies_sort_unique_array))
        end_time = time.process_time()
        print(f"analysis_sepcies_taxid方法开始运行的时间：{start_time},结束运行的时间: {end_time}, 一共运行的时间： {end_time - start_time}")
        return sepcies_sort_unique_array

    # 将计算得到的数据写入输出文件
    @write_logger()
    def write_species_count(self, sepcies_copy_sort_unique, key, specise_taxid2name):
        species_path = os.path.join(self.outdir, key+"_count.tsv")
        f = open(species_path, 'w')
        top = ["species_taxid", "species_name", "taxid", "SeqLen", "CoverLen", "Coverage", "Uniq_CoverLen",
               "Uniq_Coverage", "Reads", "Uniq_reads"]
        if len(sepcies_copy_sort_unique):
            f.write("\t".join(top) + "\n")
            sepcies_copy_sort_unique = sepcies_copy_sort_unique.tolist()
            # sepcies_sort_unique_array = np.insert(sepcies_copy_sort_unique, 0, top, axis=0)
            for item in sepcies_copy_sort_unique:
                try:
                    speceies_name = specise_taxid2name[item[0]]
                except:
                    logger.error(f"{item[0]}在name.dmp文件中找不到")
                    speceies_name = "-"
                item.insert(1, speceies_name)
                strcou = "\t".join(item)
                f.write(strcou + "\n")
        else:
            f.write("\t".join(top) + "\n")
        f.close()
        logger.debug("%s 文件写入成功" % species_path)

    # 计算唯一比上的reads数占所有高质量reads数的比值的百分比
    def calculate_read(self, sepcies_sort_unique_array, all_length):
        taxid2read = []
        if len(sepcies_sort_unique_array):
            sepcies_sort_unique_list = sepcies_sort_unique_array.tolist()
            for line_list in sepcies_sort_unique_list:
                ratio = int(line_list[-1]) / all_length * 100
                ratio = round(ratio, 2)
                taxid2read.append([line_list[1], line_list[-2], line_list[-1], ratio])
        # print(taxid2read)
        return taxid2read

    # bk和zk文件计算唯一比上的reads数占所有高质量reads数的比值的百分比
    def zkbk_calculate_read(self, sepcies_sort_unique_array, all_length):
        taxid2read = {}
        if len(sepcies_sort_unique_array):
            sepcies_sort_unique_list = sepcies_sort_unique_array.tolist()
            for line_lis in sepcies_sort_unique_list:
                ratio = int(line_lis[-1]) / all_length * 100
                ratio = round(ratio, 2)
                taxid2read[line_lis[1]] = [ratio, float(line_lis[-1])]
        return taxid2read

    # 生成发送接口的数据
    def parse_stat_data(self, other_stat_data, bk_stat_data, zk_stat_data, runid_list):
        start_time = time.process_time()
        request_data = []
        for key, value in other_stat_data.items():
            species = []
            # print(value)
            values_array = np.array(value)
            values_array = values_array[np.lexsort([values_array[:, 2].astype(int)])]
            value = values_array.tolist()
            # print(value)
            for line_list in value:
                zk_item = []
                bk_item = []
                items = {}
                try:
                    zk_item = zk_stat_data[line_list[1]]
                except:
                    pass
                try:
                    bk_item = bk_stat_data[line_list[0]]
                except:
                    pass
                items["tax_id"] = int(line_list[0])
                items["reads"] = float(line_list[1])
                items["uniq_reads"] = float(line_list[2])
                items["ratio"] = float(line_list[3])
                if len(zk_item):
                    items["zk_ratio"] = zk_item[0]
                    items["zk_uniq_reads"] = zk_item[1]
                else:
                    items["zk_ratio"] = 0
                    items["zk_uniq_reads"] = 0
                if len(bk_item):
                    items["bk_ratio"] = bk_item[0]
                    items["bk_uniq_reads"] = bk_item[1]
                else:
                    items["bk_ratio"] = 0
                    items["bk_uniq_reads"] = 0
                species.append(items)
            line_dict = {
                "runid": runid_list[key],
                "py_code": key,
                "flow_method": 2,
                "species": species
            }
            logger.info(f'{key}发送了{len(species)}的物种')
            request_data.append(line_dict)
        end_time = time.process_time()
        print(f"parse_stat_data方法开始运行的时间：{start_time},结束运行的时间: {end_time}, 一共运行的时间： {end_time - start_time}")
        return json.dumps(request_data)

    def forSyncAPI(self, request_data, data_sign, api_user):
        json = {
            'request_data': parse.quote(request_data),
            'data_sign': parse.quote(data_sign),
            'api_user': api_user
        }
        return json

    def request_api(self, request_data):
        start_time = time.process_time()
        md5 = hashlib.md5((request_data + self.api_key).encode("utf-8"))
        # base64计算结果
        data_sign = base64.b64encode(md5.hexdigest().encode("utf-8"))
        all_request_data = self.forSyncAPI(request_data, data_sign, self.api_user)
        for i in range(0, 3):
            try:
                response = requests.post(url=self.run_post, data=all_request_data, headers={})
                status_code = response.status_code
                if status_code == 200:
                    re_dict = json.loads(response.content.decode())
                    # print(re_dict)
                    logger.info(re_dict)
                    break
                else:
                    logger.error(status_code)
            except Exception as e:
                logger.error(e)
        end_time = time.process_time()
        print(f"request_api方法开始运行的时间：{start_time},结束运行的时间: {end_time}, 一共运行的时间： {end_time - start_time}")

    def run(self):
        stat_name_path_dict, runid_list = self.get_stat_path()
        refseq2speciestaxid = self.read_reseq_txt()
        specise_taxid2name = self.read_name_dmp()
        zk_stat_data = {}
        bk_stat_data = {}
        other_stat_data = {}
        for key, value in stat_name_path_dict.items():
            all_length = self.get_length_stat(value[1])
            single_list = self.read_single_data(value[0])
            sepcies_sort_unique_array = self.analysis_sepcies_taxid(single_list, refseq2speciestaxid)
            sepcies_copy_sort_unique = copy.deepcopy(sepcies_sort_unique_array)
            self.write_species_count(sepcies_copy_sort_unique, key, specise_taxid2name)
            if key.startswith("BK"):
                bk_stat_data = self.zkbk_calculate_read(sepcies_sort_unique_array, all_length)
            elif key.startswith("ZK"):
                zk_stat_data = self.zkbk_calculate_read(sepcies_sort_unique_array, all_length)
            else:
                taxid2read = self.calculate_read(sepcies_sort_unique_array, all_length)
                if len(taxid2read):
                    other_stat_data[key] = taxid2read
                else:
                    logger.info(f"{key} :物种比对结果为0")
                    # print(key)
        request_data = self.parse_stat_data(other_stat_data, bk_stat_data, zk_stat_data, runid_list)
        self.request_api(request_data)


if __name__ == "__main__":
    '''
        根据single.stat.tsv文件读取相关数据，然后依据parse_assembly脚本生成的reseq_taxid_sepcie_accession.txt文件
        读取seq_id对应的数据，进行计算计算相同基因组的覆盖率

    '''
    parser = argparse.ArgumentParser(description='计算泛基因覆盖度')
    parser.add_argument('-f', '--inputrefseqfile', help='存放输入文件的路径seq_taxid_accession.txt',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/count_depth/outfiles/seq_taxid_accession.txt")
    parser.add_argument('-s', '--statdir', help='存放single.stat.tsv或paired.stat.tsv文件的路径',
                        default="/share/data2/xujm/docker/test_pathogenic_smk/staging/190729_MN00302_0296_A000H2MFWL/")
    parser.add_argument('-n', '--sname', help='single.stat.tsv或paired.stat.tsv文件的名称',
                        default="paired.stat.tsv")
    parser.add_argument('-a', '--tname', help='single_stat.tsv或paired_stat.tsv文件的名称',
                        default="paired_stat.tsv")
    parser.add_argument('-m', '--name_dmp', help='存放names.dmp的文件',
                        default="/share/data6/PMiD/references/taxonomy/20190601/names.dmp")
    parser.add_argument('-o', '--outdir', help='输出文件的路径',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/count_depth/species_taxid_files_10.17")
    parser.add_argument('-d', '--debug', help='是否打开调试模式', default=False)
    args = parser.parse_args()
    logging_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=logging_level, filename='root.log',
                        format='[%(asctime)s] p%(process)s {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s')
    # logger.debug(args.tsvfile, args.txtpath, args.outfna)
    start = time.process_time()
    kresun = CountGenomeCoverage(args.inputrefseqfile, args.statdir, args.sname, args.tname, args.name_dmp, args.outdir)
    kresun.run()
    end =time.process_time()
    print(f"程序运行的时间:{end-start}")

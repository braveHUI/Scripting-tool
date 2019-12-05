import argparse
import logging
import os

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class BypyUpload(object):
    def __init__(self, args):
        self.inputdir = args.inputdir
        self.inputexcel = args.inputexcel
        self.input_list = args.input_list
        self.outdir = args.outdir
        self.bypy = args.bypy

    '''
        读取excel表格中的数据，获取pyid找到对应谱元id fastq文件的路径，生成list文件
    '''
    def excel_tsv(self):
        datas = pd.read_excel(self.inputexcel, header=0, index_col=None, names=None, sheet_name=1)
        datas_array = np.array(datas)
        datas_list = datas_array.tolist()
        pyid_list = [line[0] for line in datas_list]
        print(len(pyid_list))
        pyid_list.remove("CF1900074-RNA")
        pyid_list.remove("F1902194")
        print(len(pyid_list))
        out_path = os.path.join(os.path.dirname(self.outdir), 'path.txt')
        pyid_path_list = {}
        for pyid in pyid_list:
            shell = f"find {self.inputdir} -name *{pyid}* > {out_path}"
            print(shell)
            code = os.system(shell)
            if code == 0:
                with open(out_path, 'r') as file:
                    item = []
                    for path in file:
                        path = path.strip("\n")
                        if path.endswith("fastq.gz") and not "RNA" in path:
                            item.append(path)
                    pyid_path_list[pyid] = item
        os.remove(out_path)
        with open(self.input_list, 'w')as f:
            for key, value in pyid_path_list.items():
                for path in value:
                    f.write(key + "," + path + "\n")

    '''
        读取list文件，获取fastq文件的路径，在输出目录下建立链接
    '''
    def read_list(self):
        os.chdir(self.outdir)
        with open(self.input_list, 'r') as files:
            for line in files:
                line_list = line.strip("\n").split(",")
                shell = f'ln -s {line_list[1]} .'
                os.system(shell)

    '''
        把输出目录上传到百度云
    '''
    def upload(self):
        shell = f"bypy upload {self.outdir} {self.bypy}"
        print(shell)
        os.system(shell)

    def run(self):
        self.excel_tsv()
        self.read_list()
        self.upload()


if __name__ == '__main__':
    '''
    如果已经存在list文件可以注释掉第一个方法self.excel_tsv()
    '''
    parser = argparse.ArgumentParser(description='从库文件中提取序列')
    parser.add_argument('-i', '--inputdir', help='存放fastq文件的目录',
                        default="/share/data6/qinjj/Project_s756g*")
    parser.add_argument('-f', '--inputexcel', help='上传的excel表格',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/test/S19010.xlsx")
    parser.add_argument('-l', '--input_list', help='list文件的路径',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/bypy/PMid_read.list")
    parser.add_argument('-o', '--outdir', help='输出文件的目录',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/bypy/FastQ")
    parser.add_argument('-b', '--bypy', help='上传到百度云的目录',
                        default="/S19010/S19010/FastQ")
    parser.add_argument('-d', '--debug', help='是否打开调试模式', default=False)
    args = parser.parse_args()
    logging_level = logging.DEBUG if args.debug else logging.INFO
    log_path = os.path.join(args.outdir, 'root.log')
    logging.basicConfig(level=logging_level, filename=log_path,
                        format='[%(asctime)s] p%(process)s {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s')
    bu = BypyUpload(args)
    bu.run()
import argparse
import logging
import os
import numpy as np

logger = logging.getLogger(__name__)


class CountDepth(object):
    def __init__(self, inputdir, number, outfile):
        self.outfile = outfile
        self.number = int(number)
        self.inputdir = inputdir

    def get_files(self):
        txt_files = os.listdir(self.inputdir)
        files_list = [os.path.join(self.inputdir, name) for name in txt_files if name.endswith(".txt")]
        file_path_list = files_list
        print(file_path_list)
        return file_path_list

    def get_txt_data(self, file_path_list):
        data_list = []
        for path in file_path_list:
            txt_dict = {}
            with open(path, 'r') as files:
                for line in files:
                    line_list = line.strip("\n").split("\t")
                    key = line_list[0] + "_" + line_list[1]
                    try:
                        row_list = txt_dict[key]
                        print(txt_dict[key][-1])
                        row_list[-1] = line_list[-1] + row_list[-1]
                        print(txt_dict[key][-1])
                    except:
                        line_list[-1] = int(line_list[-1])
                        txt_dict[key] = line_list
            print(f"{path} : {len(txt_dict)}")
            data_list.append(txt_dict)
        return data_list

    def cat_depth_data(self, file_path_list):
        str_path = " ".join(file_path_list)
        all_depth_data = "/share/data5/hegh/project1/5.17/toolkit/Sample/bwa/most_files_bwa/all_detpth_data.txt"
        code = os.system(f"cat {str_path} >> {all_depth_data}")
        if code == 0:
            chromosome_dict = {}
            with open(all_depth_data, 'r') as files:
                for line in files:
                    line_list = line.strip("\n").split("\t")
                    chromosome_dict[line_list[0]] = 1
        print(len(chromosome_dict))


    def parse_txt_data(self, data_list):
        all_txt_dict = {}
        for txt_dict in data_list:
            for key, value in txt_dict.items():
                try:
                    row_list = all_txt_dict[key]
                    row_list[-1] += 1
                    row_list[2] += + value[2]
                except:
                    value.append(1)
                    all_txt_dict[key] = value
        filter_txt_list = []
        for key, value in all_txt_dict.items():
            if value[-1] >= self.number:
                value = [str(value[i]) for i in range(0, 3)]
                filter_txt_list.append(value)
        print(len(filter_txt_list))
        return filter_txt_list

    def write_data(self, filter_txt_list):
        with open(self.outfile, 'w') as f:
            for item in filter_txt_list:
                f.write("\t".join(item) + "\n")

    def run(self):
        file_path_list = self.get_files()
        # self.cat_depth_data(file_path_list)
        data_list = self.get_txt_data(file_path_list)
        filter_txt_list = self.parse_txt_data(data_list)
        self.write_data(filter_txt_list)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Scan Illumina Runs and send status to Promegene BMS')
    parser.add_argument('-f', '--inputdir', help='存放输入文件的目录',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/bwa/most_files_bwa/txt_dirs")
    parser.add_argument('-n', '--number', help='读取输入文件的个数',
                        default="6")
    parser.add_argument('-o', '--outfile', help='输出文件的路径',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/count_depth/filter_all_depth.txt")
    parser.add_argument('-d', '--debug', help='是否打开调试模式', default=False)
    args = parser.parse_args()
    logging_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=logging_level, filename='root.log',
                        format='[%(asctime)s] p%(process)s {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s')
    # logger.debug(args.tsvfile, args.txtpath, args.outfna)
    cd = CountDepth(args.inputdir, args.number, args.outfile)
    cd.run()
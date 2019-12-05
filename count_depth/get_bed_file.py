import argparse
import logging
import os
logger = logging.getLogger(__name__)


class GetBedFile(object):
    def __init__(self, inputfile, outfile):
        self.inputfile = inputfile
        self.outfile = outfile

    def read_depth_data(self):
        name2locus = {}
        with open(self.inputfile, 'r') as files:
            for line in files:
                line_list = line.strip("\n").split("\t")
                try:
                    row_value = name2locus[line_list[0]]
                    row_value.append(int(line_list[1]))
                except:
                    name2locus[line_list[0]] = [int(line_list[1])]
        # print(name2locus)
        return name2locus

    def parse_depth_data(self,name2locus):
        row_list = []
        for key, locus_list in name2locus.items():
            sort_locus_list = sorted(locus_list)
            start = sort_locus_list[0]
            transit = sort_locus_list[0]
            for locus in sort_locus_list:
                if locus != transit:
                    # print(locus, transit)
                    row_list.append([key, str(start), str(transit)])
                    start = locus
                    transit = locus
                transit += 1
        # print(row_list)
        return row_list

    def write_bed(self, row_list):
        with open(self.outfile, 'w') as f:
            for list in row_list:
                f.write("\t".join(list) + "\n")

    def run(self):
        name2locus = self.read_depth_data()
        row_list = self.parse_depth_data(name2locus)
        self.write_bed(row_list)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Scan Illumina Runs and send status to Promegene BMS')
    parser.add_argument('-f', '--inputfile', help='存放输入文件的目录',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/count_depth/filter_all_depth.txt")
    parser.add_argument('-o', '--outfile', help='输出文件的路径',
                        default="/share/data5/hegh/project1/5.17/toolkit/Sample/count_depth/all_depth.bed")
    parser.add_argument('-d', '--debug', help='是否打开调试模式', default=False)
    args = parser.parse_args()
    logging_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=logging_level, filename='root.log',
                        format='[%(asctime)s] p%(process)s {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s')
    # logger.debug(args.tsvfile, args.txtpath, args.outfna)
    gbf = GetBedFile(args.inputfile, args.outfile)
    gbf.run()
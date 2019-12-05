# -*- coding: utf-8 -*-
import argparse
import base64
import hashlib
import json
import logging
import os
import sys
import urllib
from datetime import date

import requests

sys.path.insert(0, '/home/liangzj/.conda/envs/lzjTest/lib/python2.7/site-packages')
from docxtpl import DocxTemplate

defaultencoding = 'utf-8'
if sys.getdefaultencoding() != defaultencoding:
    reload(sys)
    sys.setdefaultencoding(defaultencoding)

logger = logging.getLogger(__name__)
species_dict = {}


class Api(object):
    def __init__(self):
        # API账户，非真实账户
        self.api_user = "prome"
        self.api_user = os.environ.get('API_USER')
        self.api_key = os.environ.get('API_KEY')
        self.use_url = os.environ.get('API_URL')
        self.run_select = self.use_url + "/api/pmid/wikipmid/detail/"
        # self.old_run_user = "https://api.promegene.com/api/pmid/reportjson/"
        self.old_run_user = self.use_url + "/api/pmid/reportjson/"
        self.new_run_user = self.use_url + "/v1/api/pmid/reportjson/"
        self.latin_gram_list = '/share/data2/liangzj/reports/pathogenic/microbial_names_list/bacteria_genus_species_gram_list'
        self.fungi_database = '/share/data2/liangzj/reports/pathogenic/pathogen_detection_develop/fungi_database.txt'
        self.jsonfile = "/share/data5/hegh/project1/7.24/report/ZK190118_report.json"

    def forSyncAPI(self, request_data, data_sign, api_user):
        data_json = {'request_data': urllib.quote(request_data),
                     'data_sign': urllib.quote(data_sign),
                     'api_user': api_user}
        return data_json

    def request_select(self, request_data):
        request_data_json = json.dumps(request_data)
        md5 = hashlib.md5((request_data_json + self.api_key).encode("utf-8"))
        data_sign = base64.b64encode(md5.hexdigest().encode("utf-8"))
        data_json = self.forSyncAPI(request_data_json, data_sign, self.api_user)
        request = requests.post(url=self.run_select, data=data_json, headers={})
        if request.status_code == 200:
            global species_dict
            species_dict = json.loads(request.content.decode('utf-8'))
            return json.loads(request.content.decode('utf-8'))
        else:
            logging.error("can not get json! request.status_code is " + str(request.status_code))
            with open(args.j) as f:
                json_str = f.read()
                data = json.loads(json_str)
                dic1 = data[list(data.keys())[0]]
            return dic1

    # 根据api接口获取用户信息
    def request_userinfo(self, sampleID, flag):
        request_data = [{"py_code": sampleID}]
        request_data_json = json.dumps(request_data)
        md5 = hashlib.md5((request_data_json + self.api_key).encode("utf-8"))
        data_sign = base64.b64encode(md5.hexdigest().encode("utf-8"))
        data_json = self.forSyncAPI(request_data_json, data_sign, self.api_user)
        if flag == 1:
            request = requests.post(url=self.new_run_user, data=data_json, headers={})
        else:
            request = requests.post(url=self.old_run_user, data=data_json, headers={})

        if request.status_code == 200:
            data = json.loads(request.content.decode('utf-8'))
            try:
                # test1 = data['data'][0][sampleID]["sample_infor"]
                dic1 = data['data'][0][sampleID]
                # print(dic1)
                return dic1
            except Exception as e:
                logger.error(e)
        else:
            logger.error("can not get json! request.status_code is " + str(request.status_code))

    # 读取pathogen_database数据
    def read_pathogen_database(self):
        pathogen_database = {}
        file_name = os.path.abspath(__file__)
        update_path = os.path.join(os.path.dirname(file_name), 'pathogen_database.txt')
        with open(update_path) as f:
            i = 0
            for line in f:
                i += 1
                if i > 1:
                    lines = line.strip().split("\t")
                    if len(lines) < 11:
                        continue
                    pathogen_database[lines[0]] = {}
                    pathogen_database[lines[0]]['name'] = lines[1]
                    pathogen_database[lines[0]]['blood'] = lines[3]
                    pathogen_database[lines[0]]['blood_reference'] = lines[4]
                    pathogen_database[lines[0]]['sputum'] = lines[5]
                    pathogen_database[lines[0]]['sputum_reference'] = lines[6]
                    pathogen_database[lines[0]]['CSF'] = lines[7]
                    pathogen_database[lines[0]]['CSF_reference'] = lines[8]
                    pathogen_database[lines[0]]['urine'] = lines[9]
                    pathogen_database[lines[0]]['urine_reference'] = lines[10]
        return pathogen_database

    # 读取 uniqfile文件的数据
    def read_uniqfile_data(self, uniqfile_path):
        name2taxid = {}
        with open(uniqfile_path, 'r') as files:
            for line in files:
                line_list = line.strip("\n").split("\t")
                name2taxid[line_list[1]] = line_list[0]
        return name2taxid


class NewAutoFillWord(Api):
    def __init__(self, new_temple, workdir, sampleID, uniqfile_path, outdir):
        super(NewAutoFillWord, self).__init__()
        self.new_temple = new_temple
        self.content = {'bacteria_concern': [],
                        'fungi_concern': [],
                        'viruses_concern': [],
                        'parasite_concern': [],
                        'microbes': [],
                        'annotate': [],
                        'paper': []}
        self.workdir = workdir
        self.sampleID = sampleID
        self.coverage_files = self.workdir + self.sampleID + '/*/*coverage'
        self.uniqfile_path = uniqfile_path
        self.annotate_list = []
        self.outdir = outdir

    def read_stat_data(self):
        stat_path = os.path.join(self.workdir, self.sampleID+".stat")
        pathogen_database = self.read_pathogen_database()
        lname_list = []
        with open(stat_path, 'r') as f:
            temp_key = ''
            n = 0
            for line in f:
                line = line.strip()
                if len(line) > 1:
                    if line.startswith('Bacteria_species'):
                        temp_key = 'bacteria_concern'
                    elif line.startswith('Fungi_species'):
                        temp_key = 'fungi_concern'
                    elif line.startswith('Viruses_species'):
                        temp_key = 'viruses_concern'
                    elif line.startswith('Fungi_genus') or line.startswith('Viruses_genus'):
                        temp_key = ''
                    else:
                        lines = line.split("\t")
                        if temp_key and len(lines) > 2:
                            lname = lines[0]
                            coverage = lines[1]
                            shell = ['grep', '-w', lname, self.latin_gram_list]
                            shell = " ".join(shell)
                            tem_gram = os.popen(shell).readline()
                            #
                            if temp_key == 'bacteria_concern':
                                gram = 'undetermined' if not tem_gram else tem_gram.strip().split("\t")[1]
                            elif temp_key == 'fungi_concern':
                                shell = ['grep', lname, self.fungi_database]
                                shell = " ".join(shell)
                                tem_gram = os.popen(shell).readline()
                                gram = 'undetermined' if not tem_gram else tem_gram.strip().split("\t")[2]
                            elif temp_key == 'viruses_concern':
                                gram = '病毒'
                            shell = ["cat", self.coverage_files, "|", "grep", lname]
                            shell = " ".join(shell)
                            tem_cov = os.popen(shell).readlines()

                            uniq_reads = 0
                            for line in tem_cov:
                                uniq_reads += int(float(line.strip().split("\t")[-1]))
                            if uniq_reads == 0:
                                continue
                            if lname in pathogen_database.keys():
                                #
                                cname = pathogen_database[lname]['name']
                                ###judge by argv -type (blood, sputum, CSF, urine), default blood
                                blood = pathogen_database[lname]['blood'].strip()
                                blood_reference = pathogen_database[lname]['blood_reference'].strip()
                            else:
                                cname = 'undetermined'
                                ###judge by argv -type (blood, sputum, CSF, urine), default blood
                                blood = 'none'
                            #
                            concern = "*" if blood != 'none' else ""

                            if blood == 'none':
                                tem_dic = {'microbe': gram, 'name': cname, 'latin': lname, 'cov': coverage,
                                           'uniq': uniq_reads, 'con': concern}
                                self.content['microbes'].append(tem_dic)
                            else:
                                n += 1
                                if temp_key == 'bacteria_concern':
                                    tem_dic = {'gram': gram, 'name': cname, 'latin': lname, 'cov': coverage,
                                               'uniq': uniq_reads, 'con': concern}
                                    self.content['bacteria_concern'].append(tem_dic)
                                elif temp_key == 'fungi_concern':
                                    if gram == "真菌":
                                        tem_dic = {'fungi': gram, 'name': cname, 'latin': lname, 'cov': coverage,
                                                   'uniq': uniq_reads, 'con': concern}
                                        self.content['fungi_concern'].append(tem_dic)
                                    if gram == "寄生虫":
                                        tem_dic = {'parasite': gram, 'name': cname, 'latin': lname, 'cov': coverage,
                                                   'uniq': uniq_reads, 'con': concern}
                                        self.content['parasite_concern'].append(tem_dic)
                                elif temp_key == 'viruses_concern':
                                    tem_dic = {'virus': gram, 'name': cname, 'latin': lname, 'cov': coverage,
                                               'uniq': uniq_reads, 'con': concern}
                                    self.content['viruses_concern'].append(tem_dic)
                                lname_list.append(lname)
                                anno_item = {'annid': str(n), 'cname': cname, 'lname': lname, 'blood': blood}
                                self.annotate_list.append(anno_item)
                                self.content['paper'].append("[" + str(n) + "] " + blood_reference)
        content = self.content
        return content, lname_list

    def get_annocnet_content(self, lname_list, content):
        name2taxid = self.read_uniqfile_data(self.uniqfile_path)
        taxid_request_data = []
        for name in lname_list:
            try:
                item = {
                    "tax_id": name2taxid[name]
                }
                taxid_request_data.append(item)
            except Exception as e:
                logger.error(e)
        if len(species_dict):
            data = species_dict
        else:
            data = self.request_select(taxid_request_data)
        try:
            list_item = data['data']
            i = 1
            if len(list_item):
                for list in list_item:
                    if len(list["cites"]):
                        if list["cites"][0]["desc"]:
                            anno_item = {'annid': str(i), 'cname': list["name_ch"], 'lname': list["name_latin"],
                                         'blood': list["cites"][0]["desc"]}
                            content['annotate'].append(anno_item)
                            i += 1
            else:
                logger.info("编号为: {}读取本地解读信息".format(self.sampleID))
                content['annotate'] = self.annotate_list
        except:
            logger.info("编号为: {}读取本地解读信息".format(self.sampleID))
            content['annotate'] = self.annotate_list
        return content

    def parse_userinfo(self, content):
        dic1 = self.request_userinfo(self.sampleID, 1)
        if dic1 != "不存在该谱元编号的病原样品信息":
            for key, value in dic1.items():
                if dic1[key] == "None" or not dic1[key]:
                    content[key] = ""
                else:
                    try:
                        content[key] = str(value)
                    except:
                        content[key] = value.encode('utf-8')

        else:
            logger.info(dic1)
        if "is_a_infectious" in content.keys():
            content['strcheckbox'] = ""
            if content["is_a_infectious"] == "1":
                content['strcheckbox'] = "甲类传染病：☑检出   □未检出"
            elif content["is_a_infectious"] == "-1":
                content['strcheckbox'] = "甲类传染病：□检出   ☑未检出"
            else:
                content['strcheckbox'] = "甲类传染病：□检出   □未检出"

            if content["is_hiv"] == "1":
                content['strcheckbox'] += "     HIV：☑检出    □未检出"
            elif content["is_hiv"] == "-1":
                content['strcheckbox'] += "     HIV：□检出    ☑未检出"
            else:
                content['strcheckbox'] += "     HIV：□检出    □未检出"

            if content["is_syphilis"] == "1":
                content['strcheckbox'] += "      梅毒：☑检出   □未检出"
            elif content["is_syphilis"] == "-1":
                content['strcheckbox'] += "      梅毒：□检出   ☑未检出"
            else:
                content['strcheckbox'] += "      梅毒：□检出   □未检出"

            if content['gender'] == "0":
                content['gender'] = "男"
            else:
                content['gender'] = "女"
            if content['inspection_department'] == "infection":
                content['inspection_department'] = "感染"
            elif content['inspection_department'] == "digestive":
                content['inspection_department'] = "消化"
            elif content['inspection_department'] == "respiration":
                content['inspection_department'] = "呼吸"
            elif content['inspection_department'] == "icu":
                content['inspection_department'] = "ICU"
            elif content['inspection_department'] == "blood":
                content['inspection_department'] = "血液"
            elif content['inspection_department'] == "child":
                content['inspection_department'] = "儿科"
            elif content['inspection_department'] == "surgery":
                content['inspection_department'] = "外科"
            elif content['inspection_department'] == "others":
                content['inspection_department'] = "其他"

            if content['sample_type'] == "fester":
                content['sample_type'] = "脓液"
            elif content['sample_type'] == "hydrothorax_asci":
                content['sample_type'] = "胸腹水"
            elif content['sample_type'] == "microbiome":
                content['sample_type'] = "菌体"
            elif content['sample_type'] == "punctate":
                content['sample_type'] = "穿刺液"

            if content['blood_trans'] == "-1":
                content['blood_trans'] = "输  血：☑无输血   □全血   □浓缩红细胞   □洗涤红细胞   □血浆   □血小板"
            elif content['blood_trans'] == "1":
                content['blood_trans'] = "输  血：□无输血   ☑全血   □浓缩红细胞   □洗涤红细胞   □血浆   □血小板"
            elif content['blood_trans'] == "2":
                content['blood_trans'] = "输  血：□无输血   □全血   ☑浓缩红细胞   □洗涤红细胞   □血浆   □血小板"
            elif content['blood_trans'] == "3":
                content['blood_trans'] = "输  血：□无输血   □全血   □浓缩红细胞   ☑洗涤红细胞   □血浆   □血小板"
            elif content['blood_trans'] == "4":
                content['blood_trans'] = "输  血：□无输血   □全血   □浓缩红细胞   □洗涤红细胞   ☑血浆   □血小板"
            elif content['blood_trans'] == "5":
                content['blood_trans'] = "输  血：□无输血   □全血   □浓缩红细胞   □洗涤红细胞   □血浆   ☑血小板"
            else:
                content['blood_trans'] = "输  血：□无输血   □全血   □浓缩红细胞   □洗涤红细胞   □血浆   □血小板"
            # content['detection_type'] = "fungus"
        if "detection_type" in content.keys():
            if "DNA" in content['detection_type']:
                content['detection_type'] = "☑病原微生物PMiD-DNA版  □病原微生物PMiD-RNA版   □其他："
                content['other_detection_type'] = " "
            elif "RNA" in content['detection_type']:
                content['detection_type'] = "□病原微生物PMiD-DNA版  ☑病原微生物PMiD-RNA版   □其他："
                content['other_detection_type'] = " "
            elif content['detection_type'] == "bacteria":
                content['detection_type'] = "□病原微生物PMiD-DNA版  □病原微生物PMiD-RNA版   ☑其他："
                content['other_detection_type'] = "细菌类病原体检测（16S） "
            elif content['detection_type'] == "fungus":
                content['detection_type'] = "□病原微生物PMiD-DNA版  □病原微生物PMiD-RNA版   ☑其他："
                content['other_detection_type'] = "真菌类病原体检测（ITS）"
        content['datetime'] = date.today()
        content['pycode'] = self.sampleID
        content["undetected"] = ""
        if len(content["bacteria_concern"]):
            flag = False
            for item in content["bacteria_concern"]:
                if item["name"] == "肺炎克雷伯" or item["name"] == "结核分枝杆菌":
                    flag = True
                    break
            if not flag:
                content["undetected"] = "未检出结核分枝杆菌"
        else:
            content["undetected"] = "未检出细菌、结核分枝杆菌"
        if not len(content["fungi_concern"]):
            if len(content["undetected"]):
                content["undetected"] += "、真菌"
            else:
                content["undetected"] += "未检出真菌"
        if not len(content["viruses_concern"]):
            if len(content["undetected"]):
                content["undetected"] += "、病毒"
            else:
                content["undetected"] += "未检出病毒"
        if not len(content["parasite_concern"]):
            if len(content["undetected"]):
                content["undetected"] += "、寄生虫"
            else:
                content["undetected"] += "未检出寄生虫"
        return content

    def run(self):
        content, lname_list = self.read_stat_data()
        content = self.get_annocnet_content(lname_list, content)
        content = self.parse_userinfo(content)
        tpl = DocxTemplate(self.new_temple)
        tpl.render(content)
        outfile = os.path.join(self.outdir, self.sampleID + "_PMiD_V1.0_20190731.docx")
        tpl.save(outfile)


class OldAutoFillWord(Api):
    def __init__(self, old_temple, workdir, sampleID, uniqfile_path, outdir, outfile_name):
        super(OldAutoFillWord, self).__init__()
        self.old_temple = old_temple
        self.content = {'bacteria_concern': [],
                        'fungi_concern': [],
                        'viruses_concern': [],
                        'parasite_concern': [],
                        'others_concern': [],
                        'microbes': [],
                        'annotate': [],
                        'paper': []}
        self.workdir = workdir
        self.sampleID = sampleID
        self.coverage_files = self.workdir + self.sampleID + '/*/*coverage'
        self.uniqfile_path = uniqfile_path
        self.annotate_list = []
        self.outdir = outdir
        self.others_file = "others.tsv"
        self.others_name = {}
        self.outfile_name = outfile_name
        self.paper_list = []

    def read_other_data(self):
        file_name = os.path.abspath(__file__)
        other_path = os.path.join(os.path.dirname(file_name), self.others_file)
        with open(other_path) as files:
            for line in files:
                line_list = line.strip("\n").split("\t")
                self.others_name[line_list[0]] = line_list[1]

    def read_stat_data(self):
        stat_path = os.path.join(self.workdir, self.sampleID+".stat")
        pathogen_database = self.read_pathogen_database()
        lname_list = []
        other_name_list = []
        with open(stat_path, 'r') as f:
            temp_key = ''
            n = 0
            for line in f:
                line = line.strip()
                if len(line) > 1:
                    if line.startswith('Bacteria_species'):
                        temp_key = 'bacteria_concern'
                    elif line.startswith('Fungi_species'):
                        temp_key = 'fungi_concern'
                    elif line.startswith('Viruses_species'):
                        temp_key = 'viruses_concern'
                    elif line.startswith('Fungi_genus') or line.startswith('Viruses_genus'):
                        temp_key = ''
                    else:
                        lines = line.split("\t")
                        if temp_key and len(lines) > 2:
                            #
                            lname = lines[0]
                            #
                            coverage = lines[1]
                            shell = ['grep', '-w', lname, self.latin_gram_list]
                            shell = " ".join(shell)
                            tem_gram = os.popen(shell).readline()
                            #
                            if temp_key == 'bacteria_concern':
                                gram = 'undetermined' if not tem_gram else tem_gram.strip().split("\t")[1]
                            elif temp_key == 'fungi_concern':
                                shell = ['grep', lname, self.fungi_database]
                                shell = " ".join(shell)
                                tem_gram = os.popen(shell).readline()
                                gram = 'undetermined' if not tem_gram else tem_gram.strip().split("\t")[2]
                            elif temp_key == 'viruses_concern':
                                gram = '病毒'
                            shell = ["cat", self.coverage_files, "|", "grep", lname]
                            shell = " ".join(shell)
                            tem_cov = os.popen(shell).readlines()
                            #
                            uniq_reads = 0
                            for line in tem_cov:
                                uniq_reads += int(float(line.strip().split("\t")[-1]))
                            if uniq_reads == 0:
                                continue
                            if lname in pathogen_database.keys():
                                #
                                cname = pathogen_database[lname]['name']
                                ###judge by argv -type (blood, sputum, CSF, urine), default blood
                                blood = pathogen_database[lname]['blood'].strip()
                                blood_reference = pathogen_database[lname]['blood_reference'].strip()
                            else:
                                cname = 'undetermined'
                                ###judge by argv -type (blood, sputum, CSF, urine), default blood
                                blood = 'none'
                            #
                            concern = "*" if blood != 'none' else ""

                            if blood == 'none':
                                if lname in self.others_name.keys():
                                    gram = self.others_name[lname]
                                tem_dic = {'microbe': gram, 'name': cname, 'latin': lname, 'cov': coverage,
                                           'uniq': uniq_reads, 'con': concern}
                                self.content['microbes'].append(tem_dic)
                            else:
                                n += 1
                                if temp_key == 'bacteria_concern':
                                    if self.outfile_name == "PMiD_V0.5_20190924":
                                        try:
                                            other_classify = self.others_name[lname]
                                            tem_dic = {'others': other_classify, 'name': cname, 'latin': lname,
                                                       'cov': coverage,
                                                       'uniq': uniq_reads, 'con': concern}
                                            self.content['others_concern'].append(tem_dic)
                                        except:
                                            tem_dic = {'gram': gram, 'name': cname, 'latin': lname, 'cov': coverage,
                                                       'uniq': uniq_reads, 'con': concern}
                                            self.content['bacteria_concern'].append(tem_dic)
                                    else:
                                        tem_dic = {'gram': gram, 'name': cname, 'latin': lname, 'cov': coverage,
                                                   'uniq': uniq_reads, 'con': concern}
                                        self.content['bacteria_concern'].append(tem_dic)
                                elif temp_key == 'fungi_concern':
                                    if gram == "真菌":
                                        tem_dic = {'fungi': gram, 'name': cname, 'latin': lname, 'cov': coverage,
                                                   'uniq': uniq_reads, 'con': concern}
                                        self.content['fungi_concern'].append(tem_dic)
                                    if gram == "寄生虫":
                                        tem_dic = {'parasite': gram, 'name': cname, 'latin': lname, 'cov': coverage,
                                                   'uniq': uniq_reads, 'con': concern}
                                        self.content['parasite_concern'].append(tem_dic)
                                elif temp_key == 'viruses_concern':
                                    tem_dic = {'virus': gram, 'name': cname, 'latin': lname, 'cov': coverage,
                                               'uniq': uniq_reads, 'con': concern}
                                    self.content['viruses_concern'].append(tem_dic)
                                lname_list.append(lname)
                                anno_item = {'annid': str(n), 'cname': cname, 'lname': lname, 'blood': blood}
                                self.annotate_list.append(anno_item)
                                self.paper_list.append(blood_reference)
                                # self.content['paper'].append("[" + str(n) + "] " + blood_reference)

        if self.outfile_name == "PMiD_V0.5_20190924":
            self.sort_paper()
        self.content['paper'] = self.paper_list
        content = self.content
        return content, lname_list

    # 更改参考文献的顺序
    def sort_paper(self):
        other_paper = []
        original_paper = []
        for i, item in enumerate(self.annotate_list):
            if item["lname"] in self.others_name.keys():
                other_paper.append(self.paper_list[i])
            else:
                original_paper.append(self.paper_list[i])
        del self.paper_list[:]
        j = 1
        if len(original_paper):
            for blood_reference in original_paper:
                self.paper_list.append("[" + str(j) + "] " + blood_reference)
                j += 1
        if len(other_paper):
            for blood_reference in other_paper:
                self.paper_list.append("[" + str(j) + "] " + blood_reference)
                j += 1

    # 更改临床解读的顺序
    def sort_annocenet(self):
        other_items = []
        original_items = []
        for item in self.annotate_list:
            if item["lname"] in self.others_name.keys():
                other_items.append(item)
            else:
                original_items.append(item)
        del self.annotate_list[:]
        i = 1
        if len(original_items):
            for item in original_items:
                anno_item = {'annid': str(i), 'cname': item["cname"], 'lname': item["lname"], 'blood': item["blood"]}
                self.annotate_list.append(anno_item)
                i += 1
        if len(other_items):
            for item in other_items:
                anno_item = {'annid': str(i), 'cname': item["cname"], 'lname': item["lname"], 'blood': item["blood"]}
                self.annotate_list.append(anno_item)
                i += 1

    def get_annocnet_content(self, lname_list, content):
        name2taxid = self.read_uniqfile_data(self.uniqfile_path)
        taxid_request_data = []
        for name in lname_list:
            try:
                item = {
                    "tax_id": name2taxid[name]
                }
                taxid_request_data.append(item)
            except Exception as e:
                logger.error(e)
        if len(species_dict):
            data = species_dict
        else:
            data = self.request_select(taxid_request_data)
            logger.info("第一次获取物种信息")
        try:
            list_item = data['data']
            i = 1
            if len(list_item):
                other_item = []
                for list in list_item:
                    if len(list["cites"]):
                        if list["cites"][0]["desc"]:
                            if self.outfile_name == "PMiD_V0.5_20190924":
                                if list["name_latin"] in self.others_name.keys():
                                    other_item.append(list)
                                else:
                                    anno_item = {'annid': str(i), 'cname': list["name_ch"], 'lname': list["name_latin"],
                                         'blood': list["cites"][0]["desc"]}
                                    content['annotate'].append(anno_item)
                                    i += 1
                            else:
                                anno_item = {'annid': str(i), 'cname': list["name_ch"], 'lname': list["name_latin"],
                                             'blood': list["cites"][0]["desc"]}
                                content['annotate'].append(anno_item)
                            # content['annotate'].append(str(i) + ". " + list["name_ch"] + "(" + list["name_latin"] + "): " + list["cites"][0]["desc"])
                                i += 1
                if len(other_item):
                    for listo in other_item:
                        anno_item = {'annid': str(i), 'cname': listo["name_ch"], 'lname': listo["name_latin"],
                                     'blood': listo["cites"][0]["desc"]}
                        content['annotate'].append(anno_item)
                        i += 1
            else:
                logger.info("编号为: 临床解读{}读取本地解读信息".format(self.sampleID))
                if self.outfile_name == "PMiD_V0.5_20190924":
                    self.sort_annocenet()
                content['annotate'] = self.annotate_list
        except:
            logger.info("编号为: {}读取本地解读信息".format(self.sampleID))
            if self.outfile_name == "PMiD_V0.5_20190924":
                self.sort_annocenet()
            content['annotate'] = self.annotate_list
        return content

    # 发送接口,读取用户信息
    def parse_userinfo(self, content):
        dic1 = self.request_userinfo(self.sampleID, 0)
        if dic1 != "不存在该谱元编号的病原样品信息":
            for key in dic1["clini_infor"]:
                if dic1["clini_infor"][key] == "None" or not dic1["clini_infor"][key]:
                    content[key] = "    "
                else:
                    content[key] = str(dic1["clini_infor"][key])
            for key in dic1["sample_infor"]:
                if dic1["sample_infor"][key] == "None" or not dic1["sample_infor"][key]:
                    content[key] = "    "
                else:
                    content[key] = str(dic1["sample_infor"][key])
        else:
            logger.debug(dic1)
        if "Sample_type" in content.keys():
            if content['Sample_type'] == "EDTA":
                content['Sample_type'] = "EDTA抗凝全血"
            elif content['Sample_type'] == "blood":
                content['Sample_type'] = "血浆"
            elif content['Sample_type'] == "CSF":
                content['Sample_type'] = "脑脊液"
            elif content['Sample_type'] == "BAL":
                content['Sample_type'] = "肺泡灌洗液"
            elif content['Sample_type'] == "Sputum":
                content['Sample_type'] = "痰液"
            elif content['Sample_type'] == "fester":
                content['Sample_type'] = "脓液"
            elif content['Sample_type'] == "hydropsArticuli":
                content['Sample_type'] = "关节积液"
            elif content['Sample_type'] == "pleuroperitonealFluids":
                content['Sample_type'] = "胸腹水"
            elif content['Sample_type'] == "prostaticFluids":
                content['Sample_type'] = "前列腺液"
            elif content['Sample_type'] == "Urine":
                content['Sample_type'] = "尿液"
            elif content['Sample_type'] == "Tumor_puncture_fluid":
                content['Sample_type'] = "肿物穿刺液"
            elif content['Sample_type'] == "faeces":
                content['Sample_type'] = "粪便"
        content['today_date'] = date.today()
        content['Inspection_project'] = "病原微生物PMiD"
        if not len(content["bacteria_concern"]):
            content["nodetected_bacteria"] = "未检出"

        if not len(content["fungi_concern"]):
            content["nodetected_fungi"] = "未检出"

        if not len(content["viruses_concern"]):
            content["nodetected_viruses"] = "未检出"

        if not len(content["parasite_concern"]):
            content["nodetected_parasite"] = "未检出"

        if not len(content["others_concern"]):
            content["nodetected_others"] = "未检出"
        return content

    def run(self):
        self.read_other_data()
        content, lname_list = self.read_stat_data()
        content = self.get_annocnet_content(lname_list, content)
        content = self.parse_userinfo(content)
        tpl = DocxTemplate(self.old_temple)
        tpl.render(content)
        outfile = os.path.join(self.outdir, self.sampleID + "_" + self.outfile_name + ".docx")
        tpl.save(outfile)


def main(workdir, uniqfile, pyid, outdir, temple):
    temple_list = temple.split(",")
    file_name = os.path.abspath(__file__)
    for temple in temple_list:
        name = temple.split("/")[-1]
        if name.startswith("one_v0.5"):
            new_temple = os.path.join(os.path.dirname(file_name), temple + '.docx')
            naf = NewAutoFillWord(new_temple, workdir, pyid, uniqfile, outdir)
            naf.run()
        elif name.startswith("two0.5_20190220"):
            old_temple = os.path.join(os.path.dirname(file_name), temple + '.docx')
            oaf = OldAutoFillWord(old_temple, workdir, pyid, uniqfile, outdir, "PMiD_V0.5_20190220")
            oaf.run()
        elif name.startswith("three0.5_20190924"):
            old_temple = os.path.join(os.path.dirname(file_name), temple + '.docx')
            oaf = OldAutoFillWord(old_temple, workdir, pyid, uniqfile, outdir, "PMiD_V0.5_20190924")
            oaf.run()


if __name__ == '__main__':
    '''

       '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--workdir",  help="workdir, required",
                        default="/share/data2/liangzj/deliver/taoge/2019_Oct_8/6.Analysis/")
    parser.add_argument("-u", "--uniqfile", help="all_species_taxid_latin.txt.uniq",
                        default="/share/data2/xujm/PMiD/tran_16v_taxid/all_species_taxid_latin.txt.uniq")
    parser.add_argument("-j", '--jsonfile', help="json_file", default="/share/data5/hegh/project1/7.24/report/ZK190118_report.json")
    parser.add_argument("-p", '--pyid',  help="PYID", default="SA1901372")
    parser.add_argument("-o", '--outdir',  help="输出目录", default="/share/data5/hegh/project1/5.17/toolkit/Sample/combine_report/notfile")
    parser.add_argument("-t", '--temple', help="word模板名称, 目前支持3个", default="PMiD_V1.0_20190731,PMiD_V0.5_20190220,PMiD_V0.5_20190924")
    parser.add_argument('-d', '--debug', help='是否打开调试模式', default=False)
    args = parser.parse_args()
    logging_level = logging.DEBUG if args.debug else logging.INFO
    log_path = os.path.join(args.outdir, "requests.log")
    logging.basicConfig(filename=log_path, level=logging_level, format='[%(asctime)s] p%(process)s {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s')

    main(args.workdir, args.uniqfile, args.pyid, args.outdir, args.temple)



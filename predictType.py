#!/usr/bin/env python
# -*- coding:utf-8-*-
# auther:Kewei Liu time:2019/8/26 20:11  QQ:422209303  e-mail:Liukeweiaway@hotmail.com
# --------------------------------------------------------------------------
import pandas as pd
import re
from joblib import load
import Threshold
from featureEtraction import iMRMFeatureEtraction
# ------------------------------------ HG --------------------------------------------


def human_XG_iRNA_pse(sequence, species, modification, M_threshold):
    # -------------------------------------------------------------
    hg_importance_dict_sorted_list_num_Pse = [16, 22, 1441, 1389, 1692, 300, 661, 406, 1666, 1660, 1432, 1698, 1695, 338, 738, 10, 1680, 808, 1414, 204, 77, 940, 69, 1408, 74, 980, 31, 234, 8, 292, 178, 1382, 1458, 1, 1689, 13, 1373, 39, 237, 1369, 1444, 245, 66, 54, 1434, 1417, 1396]
    XG_iRNA_hg_Pse = load('model/hg/hg_xgb_Pse_10fold')
    # --------------------------------------------------------------
    line = sequence
    hg_n_base_21 = 9
    hg_n_base_list_pse_proba = []
    hg_pse_list_num = []
    M_threshold = Threshold.choiceThreshold(species, modification, M_threshold)
    # -----------------------------PSE ---------------------------------------------------
    for base in line.strip()[10:-10]:
        hg_n_base_21 += 1
        file_inter_pse = open('inter_use/hg_inter_use_pse.txt', 'w')
        if base == 'U':
            new_line = line.strip()[hg_n_base_21 - 10:hg_n_base_21 + 11]
            # ---f nc dnc tnc tenc pnc nd dnd ncp dpcp onehot dbe-----------------
            print(*iMRMFeatureEtraction.nc(new_line, 21), *iMRMFeatureEtraction.dnc(new_line, 21), *iMRMFeatureEtraction.tnc(new_line, 21), *iMRMFeatureEtraction.tenc(new_line, 21), *iMRMFeatureEtraction.pnc(new_line, 21), *iMRMFeatureEtraction.nd(new_line, 21), *iMRMFeatureEtraction.dnd(new_line, 21), *iMRMFeatureEtraction.ncp(new_line, 21), *iMRMFeatureEtraction.dpcp(new_line, 21), *iMRMFeatureEtraction.onehot(new_line, 21), *iMRMFeatureEtraction.dbe(new_line, 21), file=file_inter_pse)
            file_inter_pse.close()
            pddtest = pd.read_csv("inter_use/hg_inter_use_pse.txt", sep=' ', header=None)
            X_test = pddtest.values[:, [x - 1 for x in hg_importance_dict_sorted_list_num_Pse]]
            prediction = XG_iRNA_hg_Pse.predict_proba(X_test)
            if prediction[0][1] > M_threshold:
                hg_n_base_list_pse_proba.append(prediction[0][1])
                hg_pse_list_num.append(hg_n_base_21 + 1)
    if len(hg_pse_list_num) > 0:
        hg_pse_list_num.reverse()
        hg_n_base_list_pse_proba.reverse()
        n = 0
        for i in [i-1 for i in hg_pse_list_num]:
            line = list(line)
            line[i] = '<font color=#FF3333>U</font>(' + str(hg_n_base_list_pse_proba[n]) + ')'
            n += 1
    line = str(line).replace(',', '')
    line = line.replace('[', '')
    line = line.replace(']', '')
    line = line.replace(' ', '')
    line = line.replace('\'', '')
    line = line.replace('\\r', '')
    line_pro = line.replace('fontcolor', 'font color')
    line = line_pro
    pro_list = re.findall('\(.*?\)', line_pro)
    for i in pro_list:
        line = line.replace(i, '')
    return line, line_pro, hg_pse_list_num


def human_XG_iRNA_AI(sequence, species, modification, M_threshold):
    hg_importance_dict_sorted_list_num_AI = [1899, 1895, 5, 1, 69, 79, 252, 1892, 19, 1540, 1390, 1441, 1391, 52, 299, 1931, 1406, 25, 572, 1547, 83, 1983, 303, 1367, 72, 1414, 23, 250, 61, 110, 1584, 24, 401, 1402, 1449, 1903, 1437, 341, 214, 1894, 1833, 36, 1574, 739, 1382, 44, 1404, 1526]
    XG_iRNA_hg_AI = load('model/hg/hg_xgb_AI_10fold')

    # --------------------------------------------------------------
    line = sequence
    hg_n_base_51 = 24
    hg_n_base_list_AI_proba = []
    hg_AI_list_num = []
    M_threshold = Threshold.choiceThreshold(species, modification, M_threshold)
    # -----------------------------AI ---------------------------------------------------
    for base in line.strip()[25:-25]:
        hg_n_base_51 += 1
        file_inter_AI = open('inter_use/hg_inter_use_AI.txt', 'w')
        if base == 'A':
            new_line = line.strip()[hg_n_base_51 - 25:hg_n_base_51 + 26]
            # ---f nc dnc tnc tenc pnc nd dnd ncp dpcp onehot dbe-----------------
            print(*iMRMFeatureEtraction.nc(new_line, 51), *iMRMFeatureEtraction.dnc(new_line, 51), *iMRMFeatureEtraction.tnc(new_line, 51), *iMRMFeatureEtraction.tenc(new_line, 51), *iMRMFeatureEtraction.pnc(new_line, 51), *iMRMFeatureEtraction.nd(new_line, 51), *iMRMFeatureEtraction.dnd(new_line, 51), *iMRMFeatureEtraction.ncp(new_line, 51), *iMRMFeatureEtraction.dpcp(new_line, 51), *iMRMFeatureEtraction.onehot(new_line, 51), *iMRMFeatureEtraction.dbe(new_line, 51), file=file_inter_AI)
            file_inter_AI.close()
            pddtest = pd.read_csv("inter_use/hg_inter_use_AI.txt", sep=' ', header=None)
            X_test = pddtest.values[:, [x - 1 for x in hg_importance_dict_sorted_list_num_AI]]
            prediction = XG_iRNA_hg_AI.predict_proba(X_test)
            if prediction[0][1] > M_threshold:
                hg_n_base_list_AI_proba.append(prediction[0][1])
                hg_AI_list_num.append(hg_n_base_51 + 1)
    if len(hg_AI_list_num) > 0:
        hg_AI_list_num.reverse()
        hg_n_base_list_AI_proba.reverse()
        n = 0
        for i in [i-1 for i in hg_AI_list_num]:
            line = list(line)
            line[i] = '<font color=#9400D3>A</font>(' + str(hg_n_base_list_AI_proba[n]) + ')'
            n += 1
    line = str(line).replace(',', '')
    line = line.replace('[', '')
    line = line.replace(']', '')
    line = line.replace(' ', '')
    line = line.replace('\'', '')
    line = line.replace('\\r', '')
    line_pro = line.replace('fontcolor', 'font color')
    line = line_pro
    pro_list = re.findall('\(.*?\)', line_pro)
    for i in pro_list:
        line = line.replace(i, '')
    return line, line_pro, hg_AI_list_num


def human_XG_iRNA_M1A(sequence, species, modification, M_threshold):
    # -------------------------------------------------------------
    importance_dict_sorted_list_num_M1A = [1834, 1831, 16, 1500, 1503, 3, 9, 1, 19, 4, 20, 6, 2, 1378, 1386, 1872, 30, 1520, 1529, 1484, 243, 1002, 1372, 259, 1873, 200, 1817, 312, 1400, 584, 572, 1404, 943, 239, 41, 190, 1415, 34, 1455, 64, 1856]
    XG_iRNA_hg_M1A = load('model/hg/hg_xgb_M1A_10fold')
    # --------------------------------------------------------------
    M_threshold = Threshold.choiceThreshold(species, modification, M_threshold)
    line = sequence  # bug
    hg_n_base_41_M1A = 19
    hg_n_base_list_M1A_proba = []
    hg_M1A_list_num = []
    for base in line.strip()[20:-20]:
        hg_n_base_41_M1A += 1
        file_inter_M1A = open('inter_use/hg_inter_use_M1A.txt', 'w')
        if base == 'A':
            new_line = line.strip()[hg_n_base_41_M1A - 20:hg_n_base_41_M1A + 21]
            # ---f nc dnc tnc tenc pnc nd dnd ncp dpcp onehot dbe-----------------
            print(*iMRMFeatureEtraction.nc(new_line, 41), *iMRMFeatureEtraction.dnc(new_line, 41), *iMRMFeatureEtraction.tnc(new_line, 41), *iMRMFeatureEtraction.tenc(new_line, 41), *iMRMFeatureEtraction.pnc(new_line, 41), *iMRMFeatureEtraction.nd(new_line, 41), *iMRMFeatureEtraction.dnd(new_line, 41), *iMRMFeatureEtraction.ncp(new_line, 41), *iMRMFeatureEtraction.dpcp(new_line, 41), *iMRMFeatureEtraction.onehot(new_line, 41), *iMRMFeatureEtraction.dbe(new_line, 41), file=file_inter_M1A)
            file_inter_M1A.close()
            pddtest = pd.read_csv("inter_use/hg_inter_use_M1A.txt", sep=' ', header=None)
            X_test = pddtest.values[:, [x - 1 for x in importance_dict_sorted_list_num_M1A]]
            prediction = XG_iRNA_hg_M1A.predict_proba(X_test)
            if prediction[0][1] > M_threshold:
                hg_n_base_list_M1A_proba.append(prediction[0][1])
                hg_M1A_list_num.append(hg_n_base_41_M1A + 1)
    if len(hg_M1A_list_num) > 0:
        hg_M1A_list_num.reverse()
        hg_n_base_list_M1A_proba.reverse()
        n = 0
        for i in [i-1 for i in hg_M1A_list_num]:
            line = list(line)
            line[i] = '<font color=#1E90FF>A</font>(' + str(hg_n_base_list_M1A_proba[n]) + ')'
            n += 1
    line = str(line).replace(',', '')
    line = line.replace('[', '')
    line = line.replace(']', '')
    line = line.replace(' ', '')
    line = line.replace('\'', '')
    line = line.replace('\\r', '')
    line_pro = line.replace('fontcolor', 'font color')
    line = line_pro
    pro_list = re.findall('\(.*?\)', line_pro)
    for i in pro_list:
        line = line.replace(i, '')
    return line, line_pro, hg_M1A_list_num


def human_XG_iRNA_M6A(sequence, species, modification, M_threshold):
    importance_dict_sorted_list_num_M6A = [1834, 1500, 1831, 1503, 16, 1446, 1489, 1882, 30, 113, 997, 1545, 1857, 1892, 87, 4, 147, 1866, 209, 369, 1395, 1366, 1520, 1401, 1422, 1397, 1398, 279, 1454, 40, 243, 1417, 1498, 42, 81, 377, 1441, 51, 12, 1557, 17, 404, 96, 724]
    XG_iRNA_hg_M6A = load('model/hg/hg_xgb_M6A_10fold')

    M_threshold = Threshold.choiceThreshold(species, modification, M_threshold)
    # -------------------------------------------------------------
    line = sequence
    line = line.upper()  # bug
    if re.search(r'[^AUCG]', line.strip()):
        nota = "Please input the RNA sequence 'AUCG'."
        return nota
    length = len(line.strip())
    if length < 41:
        nota = "Please input the RNA sequence at least 51 bases! "
        return nota
    hg_n_base_41_M6A = 19
    hg_n_base_list_M6A_proba = []
    hg_M6A_list_num = []
    for base in line.strip()[20:-20]:
        hg_n_base_41_M6A += 1
        file_inter_M6A = open('inter_use/hg_inter_use_M6A.txt', 'w')
        if base == 'A':
            new_line = line.strip()[hg_n_base_41_M6A - 20:hg_n_base_41_M6A + 21]
            # ---f nc dnc tnc tenc pnc nd dnd ncp dpcp onehot dbe-----------------
            print(*iMRMFeatureEtraction.nc(new_line, 41), *iMRMFeatureEtraction.dnc(new_line, 41), *iMRMFeatureEtraction.tnc(new_line, 41), *iMRMFeatureEtraction.tenc(new_line, 41), *iMRMFeatureEtraction.pnc(new_line, 41), *iMRMFeatureEtraction.nd(new_line, 41), *iMRMFeatureEtraction.dnd(new_line, 41), *iMRMFeatureEtraction.ncp(new_line, 41), *iMRMFeatureEtraction.dpcp(new_line, 41), *iMRMFeatureEtraction.onehot(new_line, 41), *iMRMFeatureEtraction.dbe(new_line, 41), file=file_inter_M6A)
            file_inter_M6A.close()
            pddtest = pd.read_csv("inter_use/hg_inter_use_M6A.txt", sep=' ', header=None)
            X_test = pddtest.values[:, [x - 1 for x in importance_dict_sorted_list_num_M6A]]
            prediction = XG_iRNA_hg_M6A.predict_proba(X_test)
            if prediction[0][1] > M_threshold:
                hg_n_base_list_M6A_proba.append(prediction[0][1])
                hg_M6A_list_num.append(hg_n_base_41_M6A + 1)
    if len(hg_M6A_list_num) > 0:
        hg_M6A_list_num.reverse()
        hg_n_base_list_M6A_proba.reverse()
        n = 0
        for i in [i-1 for i in hg_M6A_list_num]:
            line = list(line)
            line[i] = '<font color=#00FF00>A</font>(' + str(hg_n_base_list_M6A_proba[n]) + ')'
            n += 1
    line = str(line).replace(',', '')
    line = line.replace('[', '')
    line = line.replace(']', '')
    line = line.replace(' ', '')
    line = line.replace('\'', '')
    line = line.replace('\\r', '')
    line_pro = line.replace('fontcolor', 'font color')
    line = line_pro
    pro_list = re.findall('\(.*?\)', line_pro)
    for i in pro_list:
        line = line.replace(i, '')
    return line, line_pro, hg_M6A_list_num


def human_XG_iRNA_M5C(sequence, species, modification, M_threshold):
    importance_dict_sorted_list_num_M5C = [16, 144, 1442, 82, 18, 1225, 210, 1401, 1430, 80, 57, 68, 1562, 74, 950]
    XG_iRNA_hg_M5C = load('model/hg/hg_xgb_M5C_10fold')
    M_threshold = Threshold.choiceThreshold(species, modification, M_threshold)
    # --------------------------------------------------------------
    line = sequence
    line = line.upper()  # bug
    if re.search(r'[^AUCG]', line.strip()):
        nota = "Please input the RNA sequence 'AUCG'."
        return nota
    length = len(line.strip())
    if length < 41:
        nota = "Please input the RNA sequence at least 51 bases! "
        return nota
    hg_n_base_41_M5C = 19
    hg_n_base_list_M5C_proba = []
    hg_M5C_list_num = []
    for base in line.strip()[20:-20]:
        hg_n_base_41_M5C += 1
        file_inter_M5C = open('inter_use/hg_inter_use_M5C.txt', 'w')
        if base == 'C':
            new_line = line.strip()[hg_n_base_41_M5C - 20:hg_n_base_41_M5C + 21]
            # ---f nc dnc tnc tenc pnc nd dnd ncp dpcp onehot dbe-----------------
            print(*iMRMFeatureEtraction.nc(new_line, 41), *iMRMFeatureEtraction.dnc(new_line, 41), *iMRMFeatureEtraction.tnc(new_line, 41), *iMRMFeatureEtraction.tenc(new_line, 41), *iMRMFeatureEtraction.pnc(new_line, 41), *iMRMFeatureEtraction.nd(new_line, 41), *iMRMFeatureEtraction.dnd(new_line, 41), *iMRMFeatureEtraction.ncp(new_line, 41), *iMRMFeatureEtraction.dpcp(new_line, 41), *iMRMFeatureEtraction.onehot(new_line, 41), *iMRMFeatureEtraction.dbe(new_line, 41), file=file_inter_M5C)
            file_inter_M5C.close()
            pddtest = pd.read_csv("inter_use/hg_inter_use_M5C.txt", sep=' ', header=None)
            X_test = pddtest.values[:, [x - 1 for x in importance_dict_sorted_list_num_M5C]]
            prediction = XG_iRNA_hg_M5C.predict_proba(X_test)
            if prediction[0][1] > M_threshold:
                hg_n_base_list_M5C_proba.append(prediction[0][1])
                hg_M5C_list_num.append(hg_n_base_41_M5C + 1)
    if len(hg_M5C_list_num) > 0:
        hg_M5C_list_num.reverse()
        hg_n_base_list_M5C_proba.reverse()
        n = 0
        for i in [i-1 for i in hg_M5C_list_num]:
            line = list(line)
            line[i] = '<font color=#FFA500>C</font>(' + str(hg_n_base_list_M5C_proba[n]) + ')'
            n += 1
    line = str(line).replace(',', '')
    line = line.replace('[', '')
    line = line.replace(']', '')
    line = line.replace(' ', '')
    line = line.replace('\'', '')
    line = line.replace('\\r', '')
    line_pro = line.replace('fontcolor', 'font color')
    line = line_pro
    pro_list = re.findall('\(.*?\)', line_pro)
    for i in pro_list:
        line = line.replace(i, '')
    return line, line_pro, hg_M5C_list_num


# ------------------------------------ SC --------------------------------------------
def sc_XG_iRNA_pse(sequence, species, modification, M_threshold):
    M_threshold = Threshold.choiceThreshold(species, modification, M_threshold)
    # -------------------------------------------------------------
    sc_importance_dict_sorted_list_num_Pse = [1477, 201, 1748, 1759, 1380, 1784, 10, 1121, 693, 1492, 135, 1377, 1792, 808, 1383, 251, 1800, 1766, 1746, 1, 302, 146, 78, 1488, 34, 242, 76, 230, 12, 1747, 13, 6, 269, 1753, 714, 38, 182, 77, 324, 97, 165, 1788, 664, 82, 893, 194]
    XG_iRNA_sc_Pse = load('model/sc/sc_xgb_Pse_10fold')
    # --------------------------------------------------------------
    line = sequence
    line = line.upper()  # bug
    if re.search(r'[^AUCG]', line.strip()):
        nota = "Please input the RNA sequence 'AUCG'."
        return nota
    length = len(line.strip())
    if length < 31:
        nota = "Please input the RNA sequence at least 51 bases! "
        return nota
    sc_n_base_31 = 14
    sc_n_base_list_pse_proba = []
    sc_pse_list_num = []
    # -----------------------------PSE ---------------------------------------------------
    for base in line.strip()[15:-15]:
        sc_n_base_31 += 1
        file_inter_pse = open('inter_use/sc_inter_use_pse.txt', 'w')
        if base == 'U':
            new_line = line.strip()[sc_n_base_31 - 15:sc_n_base_31 + 16]
            # ---f nc dnc tnc tenc pnc nd dnd ncp dpcp onehot dbe-----------------
            print(*iMRMFeatureEtraction.nc(new_line, 31), *iMRMFeatureEtraction.dnc(new_line, 31), *iMRMFeatureEtraction.tnc(new_line, 31), *iMRMFeatureEtraction.tenc(new_line, 31), *iMRMFeatureEtraction.pnc(new_line, 31), *iMRMFeatureEtraction.nd(new_line, 31), *iMRMFeatureEtraction.dnd(new_line, 31), *iMRMFeatureEtraction.ncp(new_line, 31), *iMRMFeatureEtraction.dpcp(new_line, 31), *iMRMFeatureEtraction.onehot(new_line, 31), *iMRMFeatureEtraction.dbe(new_line, 31), file=file_inter_pse)
            file_inter_pse.close()
            pddtest = pd.read_csv("inter_use/sc_inter_use_pse.txt", sep=' ', header=None)
            X_test = pddtest.values[:, [x - 1 for x in sc_importance_dict_sorted_list_num_Pse]]
            prediction = XG_iRNA_sc_Pse.predict_proba(X_test)
            if prediction[0][1] > M_threshold:
                sc_n_base_list_pse_proba.append(prediction[0][1])
                sc_pse_list_num.append(sc_n_base_31 + 1)
    if len(sc_pse_list_num) > 0:
        sc_pse_list_num.reverse()
        sc_n_base_list_pse_proba.reverse()
        n = 0
        for i in [i-1 for i in sc_pse_list_num]:
            line = list(line)
            line[i] = '<font color=#FF3333>U</font>(' + str(sc_n_base_list_pse_proba[n]) + ')'
            n += 1
    line = str(line).replace(',', '')
    line = line.replace('[', '')
    line = line.replace(']', '')
    line = line.replace(' ', '')
    line = line.replace('\'', '')
    line = line.replace('\\r', '')
    line_pro = line.replace('fontcolor', 'font color')
    line = line_pro
    pro_list = re.findall('\(.*?\)', line_pro)
    for i in pro_list:
        line = line.replace(i, '')
    return line, line_pro, sc_pse_list_num


def sc_XG_iRNA_M1A(sequence, species, modification, M_threshold):
    # -------------------------------------------------------------
    M_threshold = Threshold.choiceThreshold(species, modification, M_threshold)
    importance_dict_sorted_list_num_M1A = [1824, 1831, 1500]
    XG_iRNA_sc_M1A = load('model/sc/sc_xgb_M1A_10fold')
    # --------------------------------------------------------------
    line = sequence
  # bug
    if re.search(r'[^AUCG]', line.strip()):
        nota = "Please input the RNA sequence 'AUCG'."
        return nota
    length = len(line.strip())
    if length < 41:
        nota = "Please input the RNA sequence at least 51 bases! "
        return nota
    sc_n_base_41_M1A = 19
    sc_n_base_list_M1A_proba = []
    sc_M1A_list_num = []
    for base in line.strip()[20:-20]:
        sc_n_base_41_M1A += 1
        file_inter_M1A = open('inter_use/sc_inter_use_M1A.txt', 'w')
        if base == 'A':
            new_line = line.strip()[sc_n_base_41_M1A - 20:sc_n_base_41_M1A + 21]
            # ---f nc dnc tnc tenc pnc nd dnd ncp dpcp onehot dbe-----------------
            print(*iMRMFeatureEtraction.nc(new_line, 41), *iMRMFeatureEtraction.dnc(new_line, 41), *iMRMFeatureEtraction.tnc(new_line, 41), *iMRMFeatureEtraction.tenc(new_line, 41), *iMRMFeatureEtraction.pnc(new_line, 41), *iMRMFeatureEtraction.nd(new_line, 41), *iMRMFeatureEtraction.dnd(new_line, 41), *iMRMFeatureEtraction.ncp(new_line, 41), *iMRMFeatureEtraction.dpcp(new_line, 41), *iMRMFeatureEtraction.onehot(new_line, 41), *iMRMFeatureEtraction.dbe(new_line, 41), file=file_inter_M1A)
            file_inter_M1A.close()
            pddtest = pd.read_csv("inter_use/sc_inter_use_M1A.txt", sep=' ', header=None)
            X_test = pddtest.values[:, [x - 1 for x in importance_dict_sorted_list_num_M1A]]
            prediction = XG_iRNA_sc_M1A.predict_proba(X_test)
            if prediction[0][1] > M_threshold:
                sc_n_base_list_M1A_proba.append(prediction[0][1])
                sc_M1A_list_num.append(sc_n_base_41_M1A + 1)
    if len(sc_M1A_list_num) > 0:
        sc_M1A_list_num.reverse()
        sc_n_base_list_M1A_proba.reverse()
        n = 0
        for i in [i-1 for i in sc_M1A_list_num]:
            line = list(line)
            line[i] = '<font color=#1E90FF>A</font>(' + str(sc_n_base_list_M1A_proba[n]) + ')'
            n += 1
    line = str(line).replace(',', '')
    line = line.replace('[', '')
    line = line.replace(']', '')
    line = line.replace(' ', '')
    line = line.replace('\'', '')
    line = line.replace('\\r', '')
    line_pro = line.replace('fontcolor', 'font color')
    line = line_pro
    pro_list = re.findall('\(.*?\)', line_pro)
    for i in pro_list:
        line = line.replace(i, '')
    return line, line_pro, sc_M1A_list_num


def sc_XG_iRNA_M6A(sequence, species, modification, M_threshold):
    M_threshold = Threshold.choiceThreshold(species, modification, M_threshold)
    importance_dict_sorted_list_num_M6A = [1912, 1879, 1905, 1310, 327, 11, 1, 1536, 1309, 53, 18, 1311, 762, 280, 474, 4, 34, 71, 1908, 1555, 1549, 1142, 728, 186, 29, 7, 1554, 43, 46, 1117, 75, 1943, 100, 1559, 872, 126, 135, 650, 58, 1411, 981, 542, 1577, 119, 568, 216, 1441, 552]
    XG_iRNA_sc_M6A = load('model/sc/sc_xgb_M6A_10fold')
    # -------------------------------------------------------------
    line = sequence
    if re.search(r'[^AUCG]', line.strip()):
        nota = "Please input the RNA sequence 'AUCG'."
        return nota
    length = len(line.strip())
    if length < 41:
        nota = "Please input the RNA sequence at least 51 bases! "
        return nota
    sc_n_base_41_M6A = 19
    sc_n_base_list_M6A_proba = []
    sc_M6A_list_num = []
    for base in line.strip()[20:-20]:
        sc_n_base_41_M6A += 1
        file_inter_M6A = open('inter_use/sc_inter_use_M6A.txt', 'w')
        if base == 'A':
            new_line = line.strip()[sc_n_base_41_M6A - 20:sc_n_base_41_M6A + 21]
            # ---f nc dnc tnc tenc pnc nd dnd ncp dpcp onehot dbe-----------------
            print(*iMRMFeatureEtraction.nc(new_line, 41), *iMRMFeatureEtraction.dnc(new_line, 41), *iMRMFeatureEtraction.tnc(new_line, 41), *iMRMFeatureEtraction.tenc(new_line, 41), *iMRMFeatureEtraction.pnc(new_line, 41), *iMRMFeatureEtraction.nd(new_line, 41), *iMRMFeatureEtraction.dnd(new_line, 41), *iMRMFeatureEtraction.ncp(new_line, 41), *iMRMFeatureEtraction.dpcp(new_line, 41), *iMRMFeatureEtraction.onehot(new_line, 41), *iMRMFeatureEtraction.dbe(new_line, 41), file=file_inter_M6A)
            file_inter_M6A.close()
            pddtest = pd.read_csv("inter_use/sc_inter_use_M6A.txt", sep=' ', header=None)
            X_test = pddtest.values[:, [x - 1 for x in importance_dict_sorted_list_num_M6A]]
            prediction = XG_iRNA_sc_M6A.predict_proba(X_test)
            if prediction[0][1] > M_threshold:
                sc_n_base_list_M6A_proba.append(prediction[0][1])
                sc_M6A_list_num.append(sc_n_base_41_M6A + 1)
    if len(sc_M6A_list_num) > 0:
        sc_M6A_list_num.reverse()
        sc_n_base_list_M6A_proba.reverse()
        n = 0
        for i in [i-1 for i in sc_M6A_list_num]:
            line = list(line)
            line[i] = '<font color=#00FF00>A</font>(' + str(sc_n_base_list_M6A_proba[n]) + ')'
            n += 1
    line = str(line).replace(',', '')
    line = line.replace('[', '')
    line = line.replace(']', '')
    line = line.replace(' ', '')
    line = line.replace('\'', '')
    line = line.replace('\\r', '')
    line_pro = line.replace('fontcolor', 'font color')
    line = line_pro
    pro_list = re.findall('\(.*?\)', line_pro)
    for i in pro_list:
        line = line.replace(i, '')
    return line, line_pro, sc_M6A_list_num


def sc_XG_iRNA_M5C(sequence, species, modification, M_threshold):
    M_threshold = Threshold.choiceThreshold(species, modification, M_threshold)
    importance_dict_sorted_list_num_M5C = [299, 1009, 705, 1542, 1522, 1008, 1764, 4, 1862, 1200, 142, 1385, 1358, 1214, 1842, 27, 693, 1086, 193, 1527, 1904, 1533, 1554, 920, 323, 1879, 16]
    XG_iRNA_sc_M5C = load('model/sc/sc_xgb_M5C_10fold')
    # --------------------------------------------------------------
    line = sequence  # bug
    if re.search(r'[^AUCG]', line.strip()):
        nota = "Please input the RNA sequence 'AUCG'."
        return nota
    length = len(line.strip())
    if length < 41:
        nota = "Please input the RNA sequence at least 51 bases! "
        return nota
    sc_n_base_41_M5C = 19
    sc_n_base_list_M5C_proba = []
    sc_M5C_list_num = []
    for base in line.strip()[20:-20]:
        sc_n_base_41_M5C += 1
        file_inter_M5C = open('inter_use/sc_inter_use_M5C.txt', 'w')
        if base == 'C':
            new_line = line.strip()[sc_n_base_41_M5C - 20:sc_n_base_41_M5C + 21]
            # ---f nc dnc tnc tenc pnc nd dnd ncp dpcp onehot dbe-----------------
            print(*iMRMFeatureEtraction.nc(new_line, 41), *iMRMFeatureEtraction.dnc(new_line, 41), *iMRMFeatureEtraction.tnc(new_line, 41), *iMRMFeatureEtraction.tenc(new_line, 41), *iMRMFeatureEtraction.pnc(new_line, 41), *iMRMFeatureEtraction.nd(new_line, 41), *iMRMFeatureEtraction.dnd(new_line, 41), *iMRMFeatureEtraction.ncp(new_line, 41), *iMRMFeatureEtraction.dpcp(new_line, 41), *iMRMFeatureEtraction.onehot(new_line, 41), *iMRMFeatureEtraction.dbe(new_line, 41), file=file_inter_M5C)
            file_inter_M5C.close()
            pddtest = pd.read_csv("inter_use/sc_inter_use_M5C.txt", sep=' ', header=None)
            X_test = pddtest.values[:, [x - 1 for x in importance_dict_sorted_list_num_M5C]]
            prediction = XG_iRNA_sc_M5C.predict_proba(X_test)
            if prediction[0][1] > M_threshold:
                sc_n_base_list_M5C_proba.append(prediction[0][1])
                sc_M5C_list_num.append(sc_n_base_41_M5C + 1)
    if len(sc_M5C_list_num) > 0:
        sc_M5C_list_num.reverse()
        sc_n_base_list_M5C_proba.reverse()
        n = 0
        for i in [i-1 for i in sc_M5C_list_num]:
            line = list(line)
            line[i] = '<font color=#FFA500>C</font>(' + str(sc_n_base_list_M5C_proba[n]) + ')'
            n += 1
    line = str(line).replace(',', '')
    line = line.replace('[', '')
    line = line.replace(']', '')
    line = line.replace(' ', '')
    line = line.replace('\'', '')
    line = line.replace('\\r', '')
    line_pro = line.replace('fontcolor', 'font color')
    line = line_pro
    pro_list = re.findall('\(.*?\)', line_pro)
    for i in pro_list:
        line = line.replace(i, '')
    return line, line_pro, sc_M5C_list_num


# ------------------------------------ MM --------------------------------------------
def mm_XG_iRNA_pse(sequence, species, modification, M_threshold):
    M_threshold = Threshold.choiceThreshold(species, modification, M_threshold)
    # -------------------------------------------------------------
    mm_importance_dict_sorted_list_num_Pse = [1441, 1442, 85, 1692, 55, 190, 1678, 101, 25, 221, 122, 128, 60, 692, 93, 16, 255, 173, 1696, 12, 35, 701, 1715, 67, 183, 150, 808, 147, 1431, 43, 1463, 5, 1378, 71, 44, 1377, 1385, 1440, 1446, 1373, 278, 1649, 1697, 1368, 49, 662, 1457, 327, 1384, 767]
    XG_iRNA_mm_Pse = load('model/mm/mm_xgb_Pse_10fold')
    # --------------------------------------------------------------
    line = sequence  # bug
    if re.search(r'[^AUCG]', line.strip()):
        nota = "Please input the RNA sequence 'AUCG'."
        return nota
    length = len(line.strip())
    if length < 21:
        nota = "Please input the RNA sequence at least 51 bases! "
        return nota
    mm_n_base_21 = 9
    mm_n_base_list_pse_proba = []
    mm_pse_list_num = []
    # -----------------------------PSE ---------------------------------------------------
    for base in line.strip()[10:-10]:
        mm_n_base_21 += 1
        file_inter_pse = open('inter_use/mm_inter_use_pse.txt', 'w')
        if base == 'U':
            new_line = line.strip()[mm_n_base_21 - 10:mm_n_base_21 + 11]
            # ---f nc dnc tnc tenc pnc nd dnd ncp dpcp onehot dbe-----------------
            print(*iMRMFeatureEtraction.nc(new_line, 21), *iMRMFeatureEtraction.dnc(new_line, 21), *iMRMFeatureEtraction.tnc(new_line, 21), *iMRMFeatureEtraction.tenc(new_line, 21), *iMRMFeatureEtraction.pnc(new_line, 21), *iMRMFeatureEtraction.nd(new_line, 21), *iMRMFeatureEtraction.dnd(new_line, 21), *iMRMFeatureEtraction.ncp(new_line, 21), *iMRMFeatureEtraction.dpcp(new_line, 21), *iMRMFeatureEtraction.onehot(new_line, 21), *iMRMFeatureEtraction.dbe(new_line, 21), file=file_inter_pse)
            file_inter_pse.close()
            pddtest = pd.read_csv("inter_use/mm_inter_use_pse.txt", sep=' ', header=None)
            X_test = pddtest.values[:, [x - 1 for x in mm_importance_dict_sorted_list_num_Pse]]
            prediction = XG_iRNA_mm_Pse.predict_proba(X_test)
            if prediction[0][1] > M_threshold:
                mm_n_base_list_pse_proba.append(prediction[0][1])
                mm_pse_list_num.append(mm_n_base_21 + 1)
    if len(mm_pse_list_num) > 0:
        mm_pse_list_num.reverse()
        mm_n_base_list_pse_proba.reverse()
        n = 0
        for i in [i-1 for i in mm_pse_list_num]:
            line = list(line)
            line[i] = '<font color=#FF3333>U</font>(' + str(mm_n_base_list_pse_proba[n]) + ')'
            n += 1
    line = str(line).replace(',', '')
    line = line.replace('[', '')
    line = line.replace(']', '')
    line = line.replace(' ', '')
    line = line.replace('\'', '')
    line = line.replace('\\r', '')
    line_pro = line.replace('fontcolor', 'font color')
    line = line_pro
    pro_list = re.findall('\(.*?\)', line_pro)
    for i in pro_list:
        line = line.replace(i, '')
    return line, line_pro, mm_pse_list_num


def mm_XG_iRNA_M1A(sequence, species, modification, M_threshold):
    M_threshold = Threshold.choiceThreshold(species, modification, M_threshold)
    # -------------------------------------------------------------
    mm_importance_dict_sorted_list_num_M1A = [1834, 1831, 1503, 1500, 16, 87, 5, 9, 276, 1385, 1817, 19, 350, 1442, 69, 1528, 197, 1003, 1542, 4, 1530, 286, 81, 1389, 1386, 1853, 1384, 885]
    XG_iRNA_mm_M1A = load('model/mm/mm_xgb_M1A_10fold')
    # --------------------------------------------------------------
    line = sequence  # bug
    if re.search(r'[^AUCG]', line.strip()):
        nota = "Please input the RNA sequence 'AUCG'."
        return nota
    length = len(line.strip())
    if length < 41:
        nota = "Please input the RNA sequence at least 51 bases! "
        return nota
    mm_n_base_41_M1A = 19
    mm_n_base_list_M1A_proba = []
    mm_M1A_list_num = []
    for base in line.strip()[20:-20]:
        mm_n_base_41_M1A += 1
        file_inter_M1A = open('inter_use/mm_inter_use_M1A.txt', 'w')
        if base == 'A':
            new_line = line.strip()[mm_n_base_41_M1A - 20:mm_n_base_41_M1A + 21]
            # ---f nc dnc tnc tenc pnc nd dnd ncp dpcp onehot dbe-----------------
            print(*iMRMFeatureEtraction.nc(new_line, 41), *iMRMFeatureEtraction.dnc(new_line, 41), *iMRMFeatureEtraction.tnc(new_line, 41), *iMRMFeatureEtraction.tenc(new_line, 41), *iMRMFeatureEtraction.pnc(new_line, 41), *iMRMFeatureEtraction.nd(new_line, 41), *iMRMFeatureEtraction.dnd(new_line, 41), *iMRMFeatureEtraction.ncp(new_line, 41), *iMRMFeatureEtraction.dpcp(new_line, 41), *iMRMFeatureEtraction.onehot(new_line, 41), *iMRMFeatureEtraction.dbe(new_line, 41), file=file_inter_M1A)
            file_inter_M1A.close()
            pddtest = pd.read_csv("inter_use/mm_inter_use_M1A.txt", sep=' ', header=None)
            X_test = pddtest.values[:, [x - 1 for x in mm_importance_dict_sorted_list_num_M1A]]
            prediction = XG_iRNA_mm_M1A.predict_proba(X_test)
            if prediction[0][1] > M_threshold:
                mm_n_base_list_M1A_proba.append(prediction[0][1])
                mm_M1A_list_num.append(mm_n_base_41_M1A + 1)
    if len(mm_M1A_list_num) > 0:
        mm_M1A_list_num.reverse()
        mm_n_base_list_M1A_proba.reverse()
        n = 0
        for i in [i-1 for i in mm_M1A_list_num]:
            line = list(line)
            line[i] = '<font color=#1E90FF>A</font>(' + str(mm_n_base_list_M1A_proba[n]) + ')'
            n += 1
    line = str(line).replace(',', '')
    line = line.replace('[', '')
    line = line.replace(']', '')
    line = line.replace(' ', '')
    line = line.replace('\'', '')
    line = line.replace('\\r', '')
    line_pro = line.replace('fontcolor', 'font color')
    line = line_pro
    pro_list = re.findall('\(.*?\)', line_pro)
    for i in pro_list:
        line = line.replace(i, '')
    return line, line_pro, mm_M1A_list_num


def mm_XG_iRNA_M5C(sequence, species, modification, M_threshold):
    M_threshold = Threshold.choiceThreshold(species, modification, M_threshold)
    mm_importance_dict_sorted_list_num_M5C = [16, 48, 1458, 1197, 299, 859, 1359, 1879, 383, 65, 381, 82, 1519]
    XG_iRNA_mm_M5C = load('model/mm/mm_xgb_M5C_10fold')
    # --------------------------------------------------------------
    line = sequence  # bug
    if re.search(r'[^AUCG]', line.strip()):
        nota = "Please input the RNA sequence 'AUCG'."
        return nota
    length = len(line.strip())
    if length < 41:
        nota = "Please input the RNA sequence at least 51 bases! "
        return nota
    mm_n_base_41_M5C = 19
    mm_n_base_list_M5C_proba = []
    mm_M5C_list_num = []
    for base in line.strip()[20:-20]:
        mm_n_base_41_M5C += 1
        file_inter_M5C = open('inter_use/mm_inter_use_M5C.txt', 'w')
        if base == 'C':
            new_line = line.strip()[mm_n_base_41_M5C - 20:mm_n_base_41_M5C + 21]
            # ---f nc dnc tnc tenc pnc nd dnd ncp dpcp onehot dbe-----------------
            print(*iMRMFeatureEtraction.nc(new_line, 41), *iMRMFeatureEtraction.dnc(new_line, 41), *iMRMFeatureEtraction.tnc(new_line, 41), *iMRMFeatureEtraction.tenc(new_line, 41), *iMRMFeatureEtraction.pnc(new_line, 41), *iMRMFeatureEtraction.nd(new_line, 41), *iMRMFeatureEtraction.dnd(new_line, 41), *iMRMFeatureEtraction.ncp(new_line, 41), *iMRMFeatureEtraction.dpcp(new_line, 41), *iMRMFeatureEtraction.onehot(new_line, 41), *iMRMFeatureEtraction.dbe(new_line, 41), file=file_inter_M5C)
            file_inter_M5C.close()
            pddtest = pd.read_csv("inter_use/mm_inter_use_M5C.txt", sep=' ', header=None)
            X_test = pddtest.values[:, [x - 1 for x in mm_importance_dict_sorted_list_num_M5C]]
            prediction = XG_iRNA_mm_M5C.predict_proba(X_test)
            if prediction[0][1] > M_threshold:
                mm_n_base_list_M5C_proba.append(prediction[0][1])
                mm_M5C_list_num.append(mm_n_base_41_M5C + 1)
    if len(mm_M5C_list_num) > 0:
        mm_M5C_list_num.reverse()
        mm_n_base_list_M5C_proba.reverse()
        n = 0
        for i in [i-1 for i in mm_M5C_list_num]:
            line = list(line)
            line[i] = '<font color=#FFA500>C</font>(' + str(mm_n_base_list_M5C_proba[n]) + ')'
            n += 1
    line = str(line).replace(',', '')
    line = line.replace('[', '')
    line = line.replace(']', '')
    line = line.replace(' ', '')
    line = line.replace('\'', '')
    line = line.replace('\\r', '')
    line_pro = line.replace('fontcolor', 'font color')
    line = line_pro
    pro_list = re.findall('\(.*?\)', line_pro)
    for i in pro_list:
        line = line.replace(i, '')
    return line, line_pro, mm_M5C_list_num


def mm_XG_iRNA_M6A(sequence, species, modification, M_threshold):
    M_threshold = Threshold.choiceThreshold(species, modification, M_threshold)
    mm_importance_dict_sorted_list_num_M6A = [1834, 1500, 1549, 1367, 1310, 1459, 5, 1854, 1877, 58, 1843, 498, 42, 71, 1389, 700, 1775, 213, 906, 942, 81, 16, 1429, 888, 1373, 1817, 1501, 218, 124, 1380, 1436, 60, 12, 1536]
    XG_iRNA_mm_M6A = load('model/mm/mm_xgb_M6A_10fold')
    # --------------------------------------------------------------
    line = sequence  # bug
    if re.search(r'[^AUCG]', line.strip()):
        nota = "Please input the RNA sequence 'AUCG'."
        return nota
    length = len(line.strip())
    if length < 41:
        nota = "Please input the RNA sequence at least 51 bases! "
        return nota
    mm_n_base_41_M6A = 19
    mm_n_base_list_M6A_proba = []
    mm_M6A_list_num = []
    for base in line.strip()[20:-20]:
        mm_n_base_41_M6A += 1
        file_inter_M6A = open('inter_use/mm_inter_use_M6A.txt', 'w')
        if base == 'A':
            new_line = line.strip()[mm_n_base_41_M6A - 20:mm_n_base_41_M6A + 21]
            # ---f nc dnc tnc tenc pnc nd dnd ncp dpcp onehot dbe-----------------
            print(*iMRMFeatureEtraction.nc(new_line, 41), *iMRMFeatureEtraction.dnc(new_line, 41), *iMRMFeatureEtraction.tnc(new_line, 41), *iMRMFeatureEtraction.tenc(new_line, 41), *iMRMFeatureEtraction.pnc(new_line, 41), *iMRMFeatureEtraction.nd(new_line, 41), *iMRMFeatureEtraction.dnd(new_line, 41), *iMRMFeatureEtraction.ncp(new_line, 41), *iMRMFeatureEtraction.dpcp(new_line, 41), *iMRMFeatureEtraction.onehot(new_line, 41), *iMRMFeatureEtraction.dbe(new_line, 41), file=file_inter_M6A)
            file_inter_M6A.close()
            pddtest = pd.read_csv("inter_use/mm_inter_use_M6A.txt", sep=' ', header=None)
            X_test = pddtest.values[:, [x - 1 for x in mm_importance_dict_sorted_list_num_M6A]]
            prediction = XG_iRNA_mm_M6A.predict_proba(X_test)
            if prediction[0][1] > M_threshold:
                mm_n_base_list_M6A_proba.append(prediction[0][1])
                mm_M6A_list_num.append(mm_n_base_41_M6A + 1)
    if len(mm_M6A_list_num) > 0:
        mm_M6A_list_num.reverse()
        mm_n_base_list_M6A_proba.reverse()
        n = 0
        for i in [i-1 for i in mm_M6A_list_num]:
            line = list(line)
            line[i] = '<font color=#00FF00>A</font>(' + str(mm_n_base_list_M6A_proba[n]) + ')'
            n += 1
    line = str(line).replace(',', '')
    line = line.replace('[', '')
    line = line.replace(']', '')
    line = line.replace(' ', '')
    line = line.replace('\'', '')
    line = line.replace('\\r', '')
    line_pro = line.replace('fontcolor', 'font color')
    line = line_pro
    pro_list = re.findall('\(.*?\)', line_pro)
    for i in pro_list:
        line = line.replace(i, '')
    return line, line_pro, mm_M6A_list_num

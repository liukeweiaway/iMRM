#!/usr/bin/env python
# -*- coding:utf-8-*-
# author: Liukewei time:2019/1/21 QQ:422209303 e-mail:Liukeweiaway@hotmail.com
# ----------------------------------------------------------------------------
from itertools import product


def nc(seq, seq_length):
    A = round(seq.count('A') / seq_length, 3)
    U = round(seq.count('U') / seq_length, 3)
    C = round(seq.count('C') / seq_length, 3)
    G = round(seq.count('G') / seq_length, 3)
    return A, U, C, G


def dnc(seq, seq_length):
    AA = round(seq.count('AA') / seq_length, 3)
    AU = round(seq.count('AU') / seq_length, 3)
    AC = round(seq.count('AC') / seq_length, 3)
    AG = round(seq.count('AG') / seq_length, 3)
    UA = round(seq.count('UA') / seq_length, 3)
    UU = round(seq.count('UU') / seq_length, 3)
    UC = round(seq.count('UC') / seq_length, 3)
    UG = round(seq.count('UG') / seq_length, 3)
    CA = round(seq.count('CA') / seq_length, 3)
    CU = round(seq.count('CU') / seq_length, 3)
    CC = round(seq.count('CC') / seq_length, 3)
    CG = round(seq.count('CG') / seq_length, 3)
    GA = round(seq.count('GA') / seq_length, 3)
    GU = round(seq.count('GU') / seq_length, 3)
    GC = round(seq.count('GC') / seq_length, 3)
    GG = round(seq.count('GG') / seq_length, 3)
    return AA, AU, AC, AG, UA, UU, UC, UG, CA, CU, CC, CG, GA, GU, GC, GG


def tnc(seq, seq_length):
    RNA = ['A', 'U', 'C', 'G']
    tri_nucleotide_values = []
    tri_nucleotide_dict = {"".join(i): 0 for i in product(RNA, repeat=3)}
    for tri_nucleotide in tri_nucleotide_dict.keys():
        tri_nucleotide_dict[tri_nucleotide] = round(seq.count(tri_nucleotide) / seq_length, 3)
    for dict_value in tri_nucleotide_dict.values():
        tri_nucleotide_values.append(dict_value)
    return tri_nucleotide_values


def tenc(seq, seq_length):
    RNA = ['A', 'U', 'C', 'G']
    fourth_nucleotide_values = []
    fourth_nucleotide_dict = {"".join(i): 0 for i in product(RNA, repeat=4)}
    for fourth_nucleotide in fourth_nucleotide_dict.keys():
        fourth_nucleotide_dict[fourth_nucleotide] = round(seq.count(fourth_nucleotide) / seq_length, 3)
    for dict_value in fourth_nucleotide_dict.values():
        fourth_nucleotide_values.append(dict_value)
    return fourth_nucleotide_values


def pnc(seq, seq_length):
    RNA = ['A', 'U', 'C', 'G']
    fifth_nucleotide_values = []
    fifth_nucleotide_dict = {"".join(i): 0 for i in product(RNA, repeat=5)}
    for fifth_nucleotide in fifth_nucleotide_dict.keys():
        fifth_nucleotide_dict[fifth_nucleotide] = round(seq.count(fifth_nucleotide) / seq_length, 3)
    for dict_value in fifth_nucleotide_dict.values():
        fifth_nucleotide_values.append(dict_value)
    return fifth_nucleotide_values


def ncp(seq, seq_length):
    ncp_lsit = [None] * seq_length * 3
    for j in range(seq_length):
        if seq[j] == 'A':
            ncp_lsit[j * 3] = 1
            ncp_lsit[j * 3 + 1] = 1
            ncp_lsit[j * 3 + 2] = 1
        elif seq[j] == 'U':
            ncp_lsit[j * 3] = 0
            ncp_lsit[j * 3 + 1] = 0
            ncp_lsit[j * 3 + 2] = 1
        elif seq[j] == 'C':
            ncp_lsit[j * 3] = 0
            ncp_lsit[j * 3 + 1] = 1
            ncp_lsit[j * 3 + 2] = 0
        elif seq[j] == 'G':
            ncp_lsit[j * 3] = 1
            ncp_lsit[j * 3 + 1] = 0
            ncp_lsit[j * 3 + 2] = 0
    return ncp_lsit


def nd(seq, seq_length):
    nd_list = [None] * seq_length
    for j in range(seq_length):
        if seq[j] == 'A':
            nd_list[j] = round(seq[0:j].count('A') / (j + 1), 3)
        elif seq[j] == 'U':
            nd_list[j] = round(seq[0:j].count('U') / (j + 1), 3)
        elif seq[j] == 'C':
            nd_list[j] = round(seq[0:j].count('C') / (j + 1), 3)
        elif seq[j] == 'G':
            nd_list[j] = round(seq[0:j].count('G') / (j + 1), 3)
    return nd_list


def dnd(seq, seq_length):
    nnd_list = [None] * (seq_length-1)
    for j in range(seq_length-1):
        if seq[j:j + 2] == 'AA':
            nnd_list[j] = round(seq[0:j+2].count('AA') / (j+1), 3)
        elif seq[j:j + 2] == 'AU':
            nnd_list[j] = round(seq[0:j+2].count('AU') / (j+1), 3)
        elif seq[j:j + 2] == 'AC':
            nnd_list[j] = round(seq[0:j+2].count('AC') / (j+1), 3)
        elif seq[j:j + 2] == 'AG':
            nnd_list[j] = round(seq[0:j+2].count('AG') / (j+1), 3)
        elif seq[j:j + 2] == 'UA':
            nnd_list[j] = round(seq[0:j+2].count('UA') / (j+1), 3)
        elif seq[j:j + 2] == 'UU':
            nnd_list[j] = round(seq[0:j+2].count('UU') / (j+1), 3)
        elif seq[j:j + 2] == 'UC':
            nnd_list[j] = round(seq[0:j+2].count('UC') / (j+1), 3)
        elif seq[j:j + 2] == 'UG':
            nnd_list[j] = round(seq[0:j+2].count('UG') / (j+1), 3)
        elif seq[j:j + 2] == 'CA':
            nnd_list[j] = round(seq[0:j+2].count('CA') / (j+1), 3)
        elif seq[j:j + 2] == 'CU':
            nnd_list[j] = round(seq[0:j+2].count('CU') / (j+1), 3)
        elif seq[j:j + 2] == 'CC':
            nnd_list[j] = round(seq[0:j+2].count('CC') / (j+1), 3)
        elif seq[j:j + 2] == 'CG':
            nnd_list[j] = round(seq[0:j+2].count('CG') / (j+1), 3)
        elif seq[j:j + 2] == 'GA':
            nnd_list[j] = round(seq[0:j+2].count('GA') / (j+1), 3)
        elif seq[j:j + 2] == 'GU':
            nnd_list[j] = round(seq[0:j+2].count('GU') / (j+1), 3)
        elif seq[j:j + 2] == 'GC':
            nnd_list[j] = round(seq[0:j+2].count('GC') / (j+1), 3)
        elif seq[j:j + 2] == 'GG':
            nnd_list[j] = round(seq[0:j+2].count('GG') / (j+1), 3)
    return nnd_list


def onehot(seq, seq_length):
    global ot
    ONE_HOT = []
    for i in seq:
        if i == 'A':
            ot = [1, 0, 0, 0]
        elif i == 'U':
            ot = [0, 1, 0, 0]
        elif i == 'C':
            ot = [0, 0, 1, 0]
        elif i == 'G':
            ot = [0, 0, 0, 1]
        ONE_HOT = ONE_HOT + ot
    return ONE_HOT


def dbe(seq, seq_length):
    seq_length = len(seq)
    global dot
    ONE_HOT = []
    for j in range(seq_length - 1):
        if seq[j:j + 2] == 'AA':
            dot = [0, 0, 0, 0]
        elif seq[j:j + 2] == 'AU':
            dot = [0, 0, 0, 1]
        elif seq[j:j + 2] == 'AC':
            dot = [0, 0, 1, 0]
        elif seq[j:j + 2] == 'AG':
            dot = [0, 0, 1, 1]
        elif seq[j:j + 2] == 'UA':
            dot = [0, 1, 0, 0]
        elif seq[j:j + 2] == 'UU':
            dot = [0, 1, 0, 1]
        elif seq[j:j + 2] == 'UC':
            dot = [0, 1, 1, 0]
        elif seq[j:j + 2] == 'UG':
            dot = [0, 1, 1, 1]
        elif seq[j:j + 2] == 'CA':
            dot = [1, 0, 0, 0]
        elif seq[j:j + 2] == 'CU':
            dot = [1, 0, 0, 1]
        elif seq[j:j + 2] == 'CC':
            dot = [1, 0, 1, 0]
        elif seq[j:j + 2] == 'CG':
            dot = [1, 0, 1, 1]
        elif seq[j:j + 2] == 'GA':
            dot = [1, 1, 0, 0]
        elif seq[j:j + 2] == 'GU':
            dot = [1, 1, 0, 1]
        elif seq[j:j + 2] == 'GC':
            dot = [1, 1, 1, 0]
        elif seq[j:j + 2] == 'GG':
            dot = [1, 1, 1, 1]
        ONE_HOT = ONE_HOT + dot
    return ONE_HOT


def dpcp(seq, seq_length):
    PC15DNC16 = []
    Shift = [-0.08, -0.06, 0.23, -0.04, 0.07, 0.23, 0.07, -0.01, 0.11, -0.04, -0.01, 0.3, -0.02, -0.08, 0.07, 0.11]
    Slide = [-1.27, -1.36, -1.43, -1.5, -1.7, -1.43, -1.39, -1.78, -1.46, -1.5, -1.78, -1.89, -1.45, -1.27, -1.7, -1.46]
    Rise = [3.18, 3.24, 3.24, 3.3, 3.38, 3.24, 3.22, 3.32, 3.09, 3.3, 3.32, 3.3, 3.26, 3.18, 3.38, 3.09]
    Tilt = [-0.8, 1.1, 0.8, 0.5, 1.3, 0.8, 0, 0.3, 1, 0.5, 0.3, -0.1, -0.2, -0.8, 1.3, 1]
    Roll = [7, 7.1, 4.8, 8.5, 9.4, 4.8, 6.1, 12.1, 9.9, 8.5, 8.7, 12.1, 10.7, 7, 9.4, 9.9]
    Twist = [31, 33, 32, 30, 32, 32, 35, 32, 31, 30, 32, 27, 32, 31, 32, 31]
    Stacking_energy = [-13.7, -15.4, -13.8, -14, -14.2, -13.8, -16.9, -11.1, -14.4, -14, -11.1, -15.6, -16, -13.7, -14.2, -14.4]
    Enthalpy = [-6.6, -5.7, -10.2, -7.6, -13.3, -10.2, -14.2, -12.2, -10.5, -7.6, -12.2, -8, -8.1, -6.6, -10.2, -7.6]
    Entropy = [-18.4, -15.5, -26.2, -19.2, -35.5, -26.2, -34.9, -29.7, -27.8, -19.2, -29.7, -19.4, -22.6, -18.4, -26.2, -19.2]
    Free_energy = [-0.93, -1.1, -2.24, -2.08, -2.35, -2.24, -3.42, -3.26, -2.11, -2.08, -3.26, -2.36, -1.33, -0.93, -2.35,-2.11]
    Hydrophilicity = [0.04, 0.14, 0.14, 0.08, 0.1, 0.27, 0.26, 0.17, 0.21, 0.52, 0.49, 0.35, 0.21, 0.44, 0.48, 0.34]
    name_list = [Shift, Slide, Rise, Tilt, Roll, Twist, Stacking_energy, Enthalpy, Entropy, Free_energy, Hydrophilicity]
    for i in name_list:
        AA = [round(seq.count('AA')*i[0] / seq_length, 3)]
        AU = [round(seq.count('AU')*i[1] / seq_length, 3)]
        AC = [round(seq.count('AC')*i[2] / seq_length, 3)]
        AG = [round(seq.count('AG')*i[3] / seq_length, 3)]
        UA = [round(seq.count('UA')*i[4] / seq_length, 3)]
        UU = [round(seq.count('UU')*i[5] / seq_length, 3)]
        UC = [round(seq.count('UC')*i[6] / seq_length, 3)]
        UG = [round(seq.count('UG')*i[7] / seq_length, 3)]
        CA = [round(seq.count('CA')*i[8] / seq_length, 3)]
        CU = [round(seq.count('CU')*i[9] / seq_length, 3)]
        CC = [round(seq.count('CC')*i[10] / seq_length, 3)]
        CG = [round(seq.count('CG')*i[11] / seq_length, 3)]
        GA = [round(seq.count('GA')*i[12] / seq_length, 3)]
        GU = [round(seq.count('GU')*i[13] / seq_length, 3)]
        GC = [round(seq.count('GC')*i[14] / seq_length, 3)]
        GG = [round(seq.count('GG')*i[15] / seq_length, 3)]
        PC15DNC16 = PC15DNC16 + AA + AU + AC + AG + UA + UU + UC + UG + CA + CU + CC + CG + GA + GU + GC +GG
    return PC15DNC16

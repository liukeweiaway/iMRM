#!/usr/bin/env python
# -*- coding:utf-8-*-
# author: Liukewei time:2020/1/21 QQ:422209303 e-mail:Liukeweiaway@hotmail.com
# ----------------------------------------------------------------------------
import predictType


def cuttingText(seq, num):
    seq_split = []  ## 空列表
    while (seq != ''):
        seq_split.append(seq[0:num])
        seq = seq[num:]
    return seq_split


def readFile(fileName):
    with open(fileName, 'r') as f:
        file_data = f.readlines()
    return file_data[0], file_data[1]


def preprocesser(inputFile, outputFile, species, modification, threshold):
    Fasta_name, sequence = readFile(inputFile)
    sequence = sequence.upper()
    htmlf = open('iRMRModel.html', 'r', encoding="utf-8")
    outputFileR = 'results/' + outputFile
    outputFileP = 'results/' + 'Possibility_' + outputFile
    result_real = open(outputFileR, 'w', encoding="utf-8")
    result_poosi = open(outputFileP, 'w', encoding="utf-8")
    if species == 'Human':
        if modification == 'm1A':
            hg_line_M1A, hg_line_M1A_pro, hg_M1A_list = predictType.human_XG_iRNA_M1A(sequence, species, modification,
                                                                                  threshold)
            seq_split = cuttingText(sequence, 60)
            M1A_seq_ = ['-'] * len(sequence)
            for i in hg_M1A_list:
                M1A_seq_[i - 1] = 'A'  # <font color=#1E90FF>A</font>
            M1A_seq__split = cuttingText(''.join(M1A_seq_), 60)
            for line in htmlf.readlines():
                if '<li>>' not in line:
                    result_real.write(line.strip())
                    result_real.write('\n')
                elif '<li>>N1 :GUGAUAUAACUCAGUGGCAGA</li>' in line:
                    result_real.write(Fasta_name)
                    result_real.write('<br>')
                    for i in range(len(seq_split)):
                        M1A_seq = M1A_seq__split[i].replace('A', '<font color=#1E90FF>A</font>')
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        result_real.write(''.join(seq_split[i]))
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        if i < len(seq_split) - 1:
                            result_real.write(str((i + 1) * 60))
                        result_real.write('<br>')
                        result_real.write('m1A:')
                        result_real.write(''.join(M1A_seq))
                        result_real.write('<br>')
            result_poosi.write('m1A:')
            result_poosi.write('<br>')
            result_poosi.write(hg_line_M1A_pro)

        elif modification == 'm6A':
            hg_line_M6A, hg_line_M6A_pro, hg_M6A_list = predictType.human_XG_iRNA_M6A(sequence, species, modification,
                                                                                  threshold)
            seq_split = cuttingText(sequence, 60)
            M6A_seq_ = ['-'] * len(sequence)
            for i in hg_M6A_list:
                M6A_seq_[i - 1] = 'A'  # <font color=#00FF00>A</font>
            M6A_seq__split = cuttingText(''.join(M6A_seq_), 60)
            for line in htmlf.readlines():
                if '<li>>' not in line:
                    result_real.write(line.strip())
                    result_real.write('\n')
                elif '<li>>N1 :GUGAUAUAACUCAGUGGCAGA</li>' in line:
                    result_real.write(Fasta_name)
                    result_real.write('<br>')
                    for i in range(len(seq_split)):
                        M6A_seq = M6A_seq__split[i].replace('A', '<font color=#00FF00>A</font>')
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        result_real.write(''.join(seq_split[i]))
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        if i < len(seq_split) - 1:
                            result_real.write(str((i + 1) * 60))
                        result_real.write('<br>')
                        result_real.write('m6A:')
                        result_real.write(''.join(M6A_seq))
                        result_real.write('<br>')

            result_poosi.write('m6A:')
            result_poosi.write('<br>')
            result_poosi.write(hg_line_M6A_pro)
        elif modification == 'm5C':
            hg_line_M5C, hg_line_M5C_pro, hg_M5C_list = predictType.human_XG_iRNA_M5C(sequence, species, modification,
                                                                                  threshold)
            seq_split = cuttingText(sequence, 60)
            M5C_seq_ = ['-'] * len(sequence)
            for i in hg_M5C_list:
                M5C_seq_[i - 1] = 'C'  # <font color=#FF7F50>C</font>
            M5C_seq__split = cuttingText(''.join(M5C_seq_), 60)
            for line in htmlf.readlines():
                if '<li>>' not in line:
                    result_real.write(line.strip())
                    result_real.write('\n')
                elif '<li>>N1 :GUGAUAUAACUCAGUGGCAGA</li>' in line:
                    result_real.write(Fasta_name)
                    result_real.write('<br>')
                    for i in range(len(seq_split)):
                        M5C_seq = M5C_seq__split[i].replace('C', '<font color=#FF7F50>C</font>')
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        result_real.write(''.join(seq_split[i]))
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        if i < len(seq_split) - 1:
                            result_real.write(str((i + 1) * 60))
                        result_real.write('<br>')
                        result_real.write('m5C:')
                        result_real.write(''.join(M5C_seq))
                        result_real.write('<br>')

            result_poosi.write('<br>')
            result_poosi.write('m5C:')
            result_poosi.write(hg_line_M5C_pro)
        elif modification == 'pseudouridine':
            hg_line_pse, hg_line_pse_pro, hg_pse_list = predictType.human_XG_iRNA_pse(sequence, species, modification,
                                                                                  threshold)
            seq_split = cuttingText(sequence, 60)
            pse_seq_ = ['-'] * len(sequence)
            for i in hg_pse_list:
                pse_seq_[i - 1] = 'U'  # <font color=#FF7F50>C</font>
            pse_seq__split = cuttingText(''.join(pse_seq_), 60)
            for line in htmlf.readlines():
                if '<li>>' not in line:
                    result_real.write(line.strip())
                    result_real.write('\n')
                elif '<li>>N1 :GUGAUAUAACUCAGUGGCAGA</li>' in line:
                    result_real.write(Fasta_name)
                    result_real.write('<br>')
                    for i in range(len(seq_split)):
                        pse_seq = pse_seq__split[i].replace('U', '<font color=#DC143C>U</font>')
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        result_real.write(''.join(seq_split[i]))
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        if i < len(seq_split) - 1:
                            result_real.write(str((i + 1) * 60))
                        result_real.write('<br>')
                        result_real.write('Pse:')
                        result_real.write(''.join(pse_seq))
                        result_real.write('<br>')

            result_poosi.write('<br>')
            result_poosi.write('Pse:')
            result_poosi.write(hg_line_pse_pro)
        elif modification == 'A-to-I':
            hg_line_AI, hg_line_AI_pro, hg_AI_list = predictType.human_XG_iRNA_AI(sequence, species, modification,
                                                                              threshold)
            seq_split = cuttingText(sequence, 60)
            AI_seq_ = ['-'] * len(sequence)
            for i in hg_AI_list:
                AI_seq_[i - 1] = 'A'  # <font color=#FF7F50>C</font>
            AI_seq__split = cuttingText(''.join(AI_seq_), 60)
            for line in htmlf.readlines():
                if '<li>>' not in line:
                    result_real.write(line.strip())
                    result_real.write('\n')
                elif '<li>>N1 :GUGAUAUAACUCAGUGGCAGA</li>' in line:
                    result_real.write(Fasta_name)
                    result_real.write('<br>')
                    for i in range(len(seq_split)):
                        AI_seq = AI_seq__split[i].replace('A', '<font color=#9932CC>A</font>')
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        result_real.write(''.join(seq_split[i]))
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        if i < len(seq_split) - 1:
                            result_real.write(str((i + 1) * 60))
                        result_real.write('<br>')
                        result_real.write('A-I:')
                        result_real.write(''.join(AI_seq))
                        result_real.write('<br>')

            result_poosi.write('<br>')
            result_poosi.write('A-I:')
            result_poosi.write(hg_line_AI_pro)
        elif modification == 'all':
            modification = 'pseudouridine'
            hg_line_Pse, hg_line_Pse_pro, hg_Pse_list = predictType.human_XG_iRNA_pse(sequence, species, modification,
                                                                                  threshold)
            modification = 'A-to-I'
            hg_line_AI, hg_line_AI_pro, hg_AI_list = predictType.human_XG_iRNA_AI(sequence, species, modification,
                                                                              threshold)
            modification = 'm5C'
            hg_line_M5C, hg_line_M5C_pro, hg_M5C_list = predictType.human_XG_iRNA_M5C(sequence, species, modification,
                                                                                  threshold)
            modification = 'm6A'
            hg_line_M6A, hg_line_M6A_pro, hg_M6A_list = predictType.human_XG_iRNA_M6A(sequence, species, modification,
                                                                                  threshold)
            modification = 'm1A'
            hg_line_M1A, hg_line_M1A_pro, hg_M1A_list = predictType.human_XG_iRNA_M1A(sequence, species, modification,
                                                                                  threshold)
            seq_split = cuttingText(sequence, 60)
            pse_seq_ = ['-'] * len(sequence)
            AI_seq_ = ['-'] * len(sequence)
            M5C_seq_ = ['-'] * len(sequence)
            M6A_seq_ = ['-'] * len(sequence)
            M1A_seq_ = ['-'] * len(sequence)
            for i in hg_Pse_list:
                pse_seq_[i - 1] = 'U'  # <font color=#DC143C>U</font>
            pse_seq__split = cuttingText(''.join(pse_seq_), 60)

            for i in hg_AI_list:
                AI_seq_[i - 1] = 'A'  # <font color=#9932CC>A</font>
            AI_seq__split = cuttingText(''.join(AI_seq_), 60)

            for i in hg_M5C_list:
                M5C_seq_[i - 1] = 'C'  # <font color=#FF7F50>C</font>
            M5C_seq__split = cuttingText(''.join(M5C_seq_), 60)

            for i in hg_M6A_list:
                M6A_seq_[i - 1] = 'A'  # <font color=#00FF00>A</font>
            M6A_seq__split = cuttingText(''.join(M6A_seq_), 60)

            for i in hg_M1A_list:
                M1A_seq_[i - 1] = 'A'  # <font color=#1E90FF>A</font>
            M1A_seq__split = cuttingText(''.join(M1A_seq_), 60)

            for line in htmlf.readlines():
                if '<li>>' not in line:
                    result_real.write(line.strip())
                    result_real.write('\n')
                elif '<li>>N1 :GUGAUAUAACUCAGUGGCAGA</li>' in line:
                    result_real.write(Fasta_name)
                    result_real.write('<br>')
                    for i in range(len(seq_split)):
                        pse_seq = pse_seq__split[i].replace('U', '<font color=#DC143C>U</font>')
                        AI_seq = AI_seq__split[i].replace('A', '<font color=#9932CC>A</font>')
                        M5C_seq = M5C_seq__split[i].replace('C', '<font color=#FF7F50>C</font>')
                        M6A_seq = M6A_seq__split[i].replace('A', '<font color=#00FF00>A</font>')
                        M1A_seq = M1A_seq__split[i].replace('A', '<font color=#1E90FF>A</font>')
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        result_real.write(''.join(seq_split[i]))
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        if i < len(seq_split) - 1:
                            result_real.write(str((i + 1) * 60))
                        result_real.write('<br>')
                        result_real.write('PSE:')
                        result_real.write(''.join(pse_seq))
                        result_real.write('<br>')
                        result_real.write('A-I:')
                        result_real.write(''.join(AI_seq))
                        result_real.write('<br>')
                        result_real.write('m6A:')
                        result_real.write(''.join(M6A_seq))
                        result_real.write('<br>')
                        result_real.write('m1A:')
                        result_real.write(''.join(M1A_seq))
                        result_real.write('<br>')
                        result_real.write('m5C:')
                        result_real.write(''.join(M5C_seq))
                        result_real.write('<br>')

            result_poosi.write('PSE:')
            result_poosi.write(hg_line_Pse_pro)
            result_poosi.write('<br>')
            result_poosi.write('AI:')
            result_poosi.write(hg_line_AI_pro)
            result_poosi.write('<br>')
            result_poosi.write('m6A:')
            result_poosi.write(hg_line_M6A_pro)
            result_poosi.write('<br>')
            result_poosi.write('m1A:')
            result_poosi.write(hg_line_M1A_pro)
            result_poosi.write('<br>')
            result_poosi.write('m5C:')
            result_poosi.write(hg_line_M5C_pro)

    if species == 'Yeast':  # --- Yeast -----------------------------------------------------
        if modification == 'm1A':
            sc_line_M1A, sc_line_M1A_pro, sc_M1A_list = predictType.sc_XG_iRNA_M1A(sequence, species, modification,
                                                                               threshold)
            seq_split = cuttingText(sequence, 60)
            M1A_seq_ = ['-'] * len(sequence)
            for i in sc_M1A_list:
                M1A_seq_[i - 1] = 'A'  # <font color=#1E90FF>A</font>
            M1A_seq__split = cuttingText(''.join(M1A_seq_), 60)
            for line in htmlf.readlines():
                if '<li>>' not in line:
                    result_real.write(line.strip())
                    result_real.write('\n')
                elif '<li>>N1 :GUGAUAUAACUCAGUGGCAGA</li>' in line:
                    result_real.write(Fasta_name)
                    result_real.write('<br>')
                    for i in range(len(seq_split)):
                        M1A_seq = M1A_seq__split[i].replace('A', '<font color=#1E90FF>A</font>')
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        result_real.write(''.join(seq_split[i]))
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        if i < len(seq_split) - 1:
                            result_real.write(str((i + 1) * 60))
                        result_real.write('<br>')
                        result_real.write('m1A:')
                        result_real.write(''.join(M1A_seq))
                        result_real.write('<br>')
            result_poosi.write('m1A:')
            result_poosi.write(sc_line_M1A_pro)


        elif modification == 'm6A':
            sc_line_M6A, sc_line_M6A_pro, sc_M6A_list = predictType.sc_XG_iRNA_M6A(sequence, species, modification,
                                                                               threshold)
            seq_split = cuttingText(sequence, 60)
            M6A_seq_ = ['-'] * len(sequence)
            for i in sc_M6A_list:
                M6A_seq_[i - 1] = 'A'  # <font color=#00FF00>A</font>
            M6A_seq__split = cuttingText(''.join(M6A_seq_), 60)
            for line in htmlf.readlines():
                if '<li>>' not in line:
                    result_real.write(line.strip())
                    result_real.write('\n')
                elif '<li>>N1 :GUGAUAUAACUCAGUGGCAGA</li>' in line:
                    result_real.write(Fasta_name)
                    result_real.write('<br>')
                    for i in range(len(seq_split)):
                        M6A_seq = M6A_seq__split[i].replace('A', '<font color=#00FF00>A</font>')
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        result_real.write(''.join(seq_split[i]))
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        if i < len(seq_split) - 1:
                            result_real.write(str((i + 1) * 60))
                        result_real.write('<br>')
                        result_real.write('m6A:')
                        result_real.write(''.join(M6A_seq))
                        result_real.write('<br>')
            result_poosi.write('m6A:')
            result_poosi.write(sc_line_M6A_pro)

        elif modification == 'm5C':
            sc_line_M5C, sc_line_M5C_pro, sc_M5C_list = predictType.sc_XG_iRNA_M5C(sequence, species, modification,
                                                                               threshold)
            seq_split = cuttingText(sequence, 60)
            M5C_seq_ = ['-'] * len(sequence)
            for i in sc_M5C_list:
                M5C_seq_[i - 1] = 'C'  # <font color=#FF7F50>C</font>
            M5C_seq__split = cuttingText(''.join(M5C_seq_), 60)
            for line in htmlf.readlines():
                if '<li>>' not in line:
                    result_real.write(line.strip())
                    result_real.write('\n')
                elif '<li>>N1 :GUGAUAUAACUCAGUGGCAGA</li>' in line:
                    result_real.write(Fasta_name)
                    result_real.write('<br>')
                    for i in range(len(seq_split)):
                        M5C_seq = M5C_seq__split[i].replace('C', '<font color=#FF7F50>C</font>')
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        result_real.write(''.join(seq_split[i]))
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        if i < len(seq_split) - 1:
                            result_real.write(str((i + 1) * 60))
                        result_real.write('<br>')
                        result_real.write('m5C:')
                        result_real.write(''.join(M5C_seq))
                        result_real.write('<br>')

            result_poosi.write('m5C:')
            result_poosi.write(sc_line_M5C_pro)

        elif modification == 'pseudouridine':
            sc_line_pse, sc_line_pse_pro, sc_pse_list = predictType.sc_XG_iRNA_pse(sequence, species, modification,
                                                                               threshold)
            seq_split = cuttingText(sequence, 60)
            pse_seq_ = ['-'] * len(sequence)
            for i in sc_pse_list:
                pse_seq_[i - 1] = 'U'  # <font color=#FF7F50>C</font>
            pse_seq__split = cuttingText(''.join(pse_seq_), 60)
            for line in htmlf.readlines():
                if '<li>>' not in line:
                    result_real.write(line.strip())
                    result_real.write('\n')
                elif '<li>>N1 :GUGAUAUAACUCAGUGGCAGA</li>' in line:
                    result_real.write(Fasta_name)
                    result_real.write('<br>')
                    for i in range(len(seq_split)):
                        pse_seq = pse_seq__split[i].replace('U', '<font color=#DC143C>U</font>')
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        result_real.write(''.join(seq_split[i]))
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        if i < len(seq_split) - 1:
                            result_real.write(str((i + 1) * 60))
                        result_real.write('<br>')
                        result_real.write('PSE:')
                        result_real.write(''.join(pse_seq))
                        result_real.write('<br>')

            result_poosi.write('PSE:')
            result_poosi.write(sc_line_pse_pro)

        elif modification == 'all':
            modification = 'pseudouridine'
            sc_line_Pse, sc_line_Pse_pro, sc_Pse_list = predictType.sc_XG_iRNA_pse(sequence, species, modification,
                                                                               threshold)
            modification = 'm5C'
            sc_line_M5C, sc_line_M5C_pro, sc_M5C_list = predictType.sc_XG_iRNA_M5C(sequence, species, modification,
                                                                               threshold)
            modification = 'm6A'
            sc_line_M6A, sc_line_M6A_pro, sc_M6A_list = predictType.sc_XG_iRNA_M6A(sequence, species, modification,
                                                                               threshold)
            modification = 'm1A'
            sc_line_M1A, sc_line_M1A_pro, sc_M1A_list = predictType.sc_XG_iRNA_M1A(sequence, species, modification,
                                                                               threshold)
            seq_split = cuttingText(sequence, 60)
            pse_seq_ = ['-'] * len(sequence)
            M5C_seq_ = ['-'] * len(sequence)
            M6A_seq_ = ['-'] * len(sequence)
            M1A_seq_ = ['-'] * len(sequence)
            for i in sc_Pse_list:
                pse_seq_[i - 1] = 'U'  # <font color=#DC143C>U</font>
            pse_seq__split = cuttingText(''.join(pse_seq_), 60)

            for i in sc_M5C_list:
                M5C_seq_[i - 1] = 'C'  # <font color=#FF7F50>C</font>
            M5C_seq__split = cuttingText(''.join(M5C_seq_), 60)

            for i in sc_M6A_list:
                M6A_seq_[i - 1] = 'A'  # <font color=#00FF00>A</font>
            M6A_seq__split = cuttingText(''.join(M6A_seq_), 60)

            for i in sc_M1A_list:
                M1A_seq_[i - 1] = 'A'  # <font color=#1E90FF>A</font>
            M1A_seq__split = cuttingText(''.join(M1A_seq_), 60)

            for line in htmlf.readlines():
                if '<li>>' not in line:
                    result_real.write(line.strip())
                    result_real.write('\n')
                elif '<li>>N1 :GUGAUAUAACUCAGUGGCAGA</li>' in line:
                    result_real.write(Fasta_name)
                    result_real.write('<br>')
                    for i in range(len(seq_split)):
                        pse_seq = pse_seq__split[i].replace('U', '<font color=#DC143C>U</font>')
                        M5C_seq = M5C_seq__split[i].replace('C', '<font color=#FF7F50>C</font>')
                        M6A_seq = M6A_seq__split[i].replace('A', '<font color=#00FF00>A</font>')
                        M1A_seq = M1A_seq__split[i].replace('A', '<font color=#1E90FF>A</font>')
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        result_real.write(''.join(seq_split[i]))
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        if i < len(seq_split) - 1:
                            result_real.write(str((i + 1) * 60))
                        result_real.write('<br>')
                        result_real.write('PSE:')
                        result_real.write(''.join(pse_seq))
                        result_real.write('<br>')
                        result_real.write('m6A:')
                        result_real.write(''.join(M6A_seq))
                        result_real.write('<br>')
                        result_real.write('m1A:')
                        result_real.write(''.join(M1A_seq))
                        result_real.write('<br>')
                        result_real.write('m5C:')
                        result_real.write(''.join(M5C_seq))
                        result_real.write('<br>')

            result_poosi.write('PSE:')
            result_poosi.write(sc_line_Pse_pro)
            result_poosi.write('<br>')
            result_poosi.write('m6A:')
            result_poosi.write(sc_line_M6A_pro)
            result_poosi.write('<br>')
            result_poosi.write('m1A:')
            result_poosi.write(sc_line_M1A_pro)
            result_poosi.write('<br>')
            result_poosi.write('m5C:')
            result_poosi.write(sc_line_M5C_pro)

    if species == 'Mouse':  # Mouse
        if modification == 'm1A':
            mm_line_M1A, mm_line_M1A_pro, mm_M1A_list = predictType.mm_XG_iRNA_M1A(sequence, species, modification,
                                                                               threshold)
            seq_split = cuttingText(sequence, 60)
            M1A_seq_ = ['-'] * len(sequence)
            for i in mm_M1A_list:
                M1A_seq_[i - 1] = 'A'  # <font color=#1E90FF>A</font>
            M1A_seq__split = cuttingText(''.join(M1A_seq_), 60)
            for line in htmlf.readlines():
                if '<li>>' not in line:
                    result_real.write(line.strip())
                    result_real.write('\n')
                elif '<li>>N1 :GUGAUAUAACUCAGUGGCAGA</li>' in line:
                    result_real.write(Fasta_name)
                    result_real.write('<br>')
                    for i in range(len(seq_split)):
                        M1A_seq = M1A_seq__split[i].replace('A', '<font color=#1E90FF>A</font>')
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        result_real.write(''.join(seq_split[i]))
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        if i < len(seq_split) - 1:
                            result_real.write(str((i + 1) * 60))
                        result_real.write('<br>')
                        result_real.write('m1A:')
                        result_real.write(''.join(M1A_seq))
                        result_real.write('<br>')

            result_poosi.write('m1A:')
            result_poosi.write(mm_line_M1A_pro)

        elif modification == 'm5C':
            mm_line_M5C, mm_line_M5C_pro, mm_M5C_list = predictType.mm_XG_iRNA_M5C(sequence, species, modification,
                                                                               threshold)
            seq_split = cuttingText(sequence, 60)
            M5C_seq_ = ['-'] * len(sequence)
            for i in mm_M5C_list:
                M5C_seq_[i - 1] = 'C'  # <font color=#FF7F50>C</font>
            M5C_seq__split = cuttingText(''.join(M5C_seq_), 60)
            for line in htmlf.readlines():
                if '<li>>' not in line:
                    result_real.write(line.strip())
                    result_real.write('\n')
                elif '<li>>N1 :GUGAUAUAACUCAGUGGCAGA</li>' in line:
                    result_real.write(Fasta_name)
                    result_real.write('<br>')
                    for i in range(len(seq_split)):
                        M5C_seq = M5C_seq__split[i].replace('C', '<font color=#FF7F50>C</font>')
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        result_real.write(''.join(seq_split[i]))
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        if i < len(seq_split) - 1:
                            result_real.write(str((i + 1) * 60))
                        result_real.write('<br>')
                        result_real.write('m5C:')
                        result_real.write(''.join(M5C_seq))
                        result_real.write('<br>')
            result_poosi.write('m5C:')
            result_poosi.write(mm_line_M5C_pro)

        elif modification == 'pseudouridine':
            mm_line_pse, mm_line_pse_pro, mm_pse_list = predictType.mm_XG_iRNA_pse(sequence, species, modification,
                                                                               threshold)
            seq_split = cuttingText(sequence, 60)
            pse_seq_ = ['-'] * len(sequence)
            for i in mm_pse_list:
                pse_seq_[i - 1] = 'U'  # <font color=#FF7F50>C</font>
            pse_seq__split = cuttingText(''.join(pse_seq_), 60)
            for line in htmlf.readlines():
                if '<li>>' not in line:
                    result_real.write(line.strip())
                    result_real.write('\n')
                elif '<li>>N1 :GUGAUAUAACUCAGUGGCAGA</li>' in line:
                    result_real.write(Fasta_name)
                    result_real.write('<br>')
                    for i in range(len(seq_split)):
                        pse_seq = pse_seq__split[i].replace('U', '<font color=#DC143C>U</font>')
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        result_real.write(''.join(seq_split[i]))
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        if i < len(seq_split) - 1:
                            result_real.write(str((i + 1) * 60))
                        result_real.write('<br>')
                        result_real.write('PSE:')
                        result_real.write(''.join(pse_seq))
                        result_real.write('<br>')

            result_poosi.write('PSE:')
            result_poosi.write(mm_line_pse_pro)

        elif modification == 'm6A':
            mm_line_M6A, mm_line_M6A_pro, mm_M6A_list = predictType.mm_XG_iRNA_M6A(sequence, species, modification,
                                                                               threshold)
            seq_split = cuttingText(sequence, 60)
            M6A_seq_ = ['-'] * len(sequence)
            for i in mm_M6A_list:
                M6A_seq_[i - 1] = 'A'  # <font color=#00FF00>A</font>
            M6A_seq__split = cuttingText(''.join(M6A_seq_), 60)
            for line in htmlf.readlines():
                if '<li>>' not in line:
                    result_real.write(line.strip())
                    result_real.write('\n')
                elif '<li>>N1 :GUGAUAUAACUCAGUGGCAGA</li>' in line:
                    result_real.write(Fasta_name)
                    result_real.write('<br>')
                    for i in range(len(seq_split)):
                        M6A_seq = M6A_seq__split[i].replace('A', '<font color=#00FF00>A</font>')
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        result_real.write(''.join(seq_split[i]))
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        if i < len(seq_split) - 1:
                            result_real.write(str((i + 1) * 60))
                        result_real.write('<br>')
                        result_real.write('m6A:')
                        result_real.write(''.join(M6A_seq))
                        result_real.write('<br>')

            result_poosi.write('m6A:')
            result_poosi.write(mm_line_M6A_pro)

        elif modification == 'all':
            modification = 'pseudouridine'
            mm_line_Pse, mm_line_Pse_pro, mm_Pse_list = predictType.mm_XG_iRNA_pse(sequence, species, modification,
                                                                               threshold)
            modification = 'm5C'
            mm_line_M5C, mm_line_M5C_pro, mm_M5C_list = predictType.mm_XG_iRNA_M5C(sequence, species, modification,
                                                                               threshold)
            modification = 'm1A'
            mm_line_M1A, mm_line_M1A_pro, mm_M1A_list = predictType.mm_XG_iRNA_M1A(sequence, species, modification,
                                                                               threshold)
            modification = 'm6A'
            mm_line_M6A, mm_line_M6A_pro, mm_M6A_list = predictType.mm_XG_iRNA_M6A(sequence, species, modification,
                                                                               threshold)
            seq_split = cuttingText(sequence, 60)
            pse_seq_ = ['-'] * len(sequence)
            M5C_seq_ = ['-'] * len(sequence)
            M1A_seq_ = ['-'] * len(sequence)
            M6A_seq_ = ['-'] * len(sequence)
            for i in mm_Pse_list:
                pse_seq_[i - 1] = 'U'  # <font color=#DC143C>U</font>
            pse_seq__split = cuttingText(''.join(pse_seq_), 60)

            for i in mm_M5C_list:
                M5C_seq_[i - 1] = 'C'  # <font color=#FF7F50>C</font>
            M5C_seq__split = cuttingText(''.join(M5C_seq_), 60)

            for i in mm_M1A_list:
                M1A_seq_[i - 1] = 'A'  # <font color=#1E90FF>A</font>
            M1A_seq__split = cuttingText(''.join(M1A_seq_), 60)

            for i in mm_M6A_list:
                M6A_seq_[i - 1] = 'A'  # <font color=#1E90FF>A</font>
            M6A_seq__split = cuttingText(''.join(M6A_seq_), 60)

            for line in htmlf.readlines():
                if '<li>>' not in line:
                    result_real.write(line.strip())
                    result_real.write('\n')
                elif '<li>>N1 :GUGAUAUAACUCAGUGGCAGA</li>' in line:
                    result_real.write(Fasta_name)
                    result_real.write('<br>')
                    for i in range(len(seq_split)):
                        pse_seq = pse_seq__split[i].replace('U', '<font color=#DC143C>U</font>')
                        M5C_seq = M5C_seq__split[i].replace('C', '<font color=#FF7F50>C</font>')
                        M1A_seq = M1A_seq__split[i].replace('A', '<font color=#1E90FF>A</font>')
                        M6A_seq = M6A_seq__split[i].replace('A', '<font color=#00FF00>A</font>')
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        result_real.write(''.join(seq_split[i]))
                        result_real.write('&nbsp&nbsp&nbsp&nbsp')
                        if i < len(seq_split) - 1:
                            result_real.write(str((i + 1) * 60))
                        result_real.write('<br>')
                        result_real.write('PSE:')
                        result_real.write(''.join(pse_seq))
                        result_real.write('<br>')
                        result_real.write('m6A:')
                        result_real.write(''.join(M6A_seq))
                        result_real.write('<br>')
                        result_real.write('m1A:')
                        result_real.write(''.join(M1A_seq))
                        result_real.write('<br>')
                        result_real.write('m5C:')
                        result_real.write(''.join(M5C_seq))
                        result_real.write('<br>')

            result_poosi.write('PSE:')
            result_poosi.write(mm_line_Pse_pro)
            result_poosi.write('<br>')
            result_poosi.write('m6A:')
            result_poosi.write(mm_line_M6A_pro)
            result_poosi.write('<br>')
            result_poosi.write('m1A:')
            result_poosi.write(mm_line_M1A_pro)
            result_poosi.write('<br>')
            result_poosi.write('m5C:')
            result_poosi.write(mm_line_M5C_pro)

    htmlf.close()
    result_real.close()
    result_poosi.close()

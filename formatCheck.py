#!/usr/bin/env python
# -*- coding:utf-8-*-
# author: Liukewei time:2020/1/21 QQ:422209303 e-mail:Liukeweiaway@hotmail.com
# ----------------------------------------------------------------------------
from preproces import readFile
import re


def formatCheck(inputFile, outputFile, species, modification, threshold):
    Fasta_name, sequence = readFile(inputFile)
    speciesList = ['Human', 'Mouse', 'Yeast']
    modificationList = ['m1A', 'm6A', 'm5C', 'A-to-I', 'pseudouridine', 'all']
    thresholdList = ['low', 'normal', 'high']
    line = sequence.upper()  # bug
    if re.search(r'[^AUCG]', line.strip()):
        note = "Please input the RNA sequence 'AUCG' or 'aucg'."
        return note
    length = len(line.strip())
    if length < 51:
        note = "Please input the RNA sequence at least 51 bases! "
        return note
    if species not in speciesList:
        note = "eror: please input: -s Human/Yeast/Mouse"
        return note
    if modification not in modificationList:
        note = "eror: please input: -m m1A/m6A/m5C/A-to-I/pseudouridine/all"
        return note
    if threshold not in thresholdList:
        note = "eror: please input: -t low/normal/high"
        return note
    note = 'checking over'
    return note


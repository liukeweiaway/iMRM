#!/usr/bin/env python
# -*- coding:utf-8-*-
# author: Liukewei time:2020/1/21 QQ:422209303 e-mail:Liukeweiaway@hotmail.com
# ----------------------------------------------------------------------------
import argparse
import preproces
import warnings
from formatCheck import formatCheck
warnings.filterwarnings("ignore")


def parse_args():
    """
    :return:进行参数的解析
    """
    description = "iMRMF is able to simultaneously identify m6A, m5C, m1A, ψ and A-to-I modifications in Homo sapiens, Mus musculus and Saccharomyces cerevisiae. Thank you very much for submitting the error message to Liukeweiaway@hotmail.com.  Example: python iMRM.py -i sequence.txt -o ccc.html -s Human -m all -t normal "
    parser = argparse.ArgumentParser(description=description)
    help = "The correct path of address"
    parser.add_argument('--addresses', help=help)
    parser.add_argument('-i', '--inputFile', help='-i input.txt (The input file is a complete Fasta format sequence.)')
    parser.add_argument('-o', '--outputFile', help='-o output.html (Results are saved under results folder.)')
    parser.add_argument('-s', '--species', help='-m Human/Mouse/Yeast (Choose one from three species to use.)')
    parser.add_argument('-m', '--modification', help='-m m6A/m5C/A-toI/pseudouridine/m1A/all (We can choose the predicted modification site, or select all to predict multiple modifications together. It should be noted that A-to-I modification is limited by the amount of data and can only be predicted in human.)')
    parser.add_argument('-t', '--threshold', help='-m low/normal/high (We offer 3 options based on the difference in specificity, which are low, normal and high.)')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    formatCheck = formatCheck(args.inputFile, args.outputFile, args.species, args.modification, args.threshold)
    print(formatCheck)
    preproces.preprocesser(args.inputFile, args.outputFile, args.species, args.modification, args.threshold)
    # previous_deal(args.inputFile, args.outputFile, args.featureList)  # feature extraction

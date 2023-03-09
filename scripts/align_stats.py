#
# ============================================================================
#
#        Authors: Christina SHI <Christina.hshi@gmail.com>
#
# ============================================================================
#

import math
import random
import scipy.stats
import datetime
import argparse
import numpy as np
import sys
import csv
import re
import os
from os import listdir
from os import path
import pandas as pd

# from os.path import isfile, join

def getArguments():
  parser = argparse.ArgumentParser(description="Compute statistics of OM aligns.")
  parser.add_argument('-r', '--ref', required=True, help="Reference map file in CMAP format")
  parser.add_argument('-q', '--query', required=True, help='Query map file in CMAP format')
  parser.add_argument('-a', '--align', required=True, help="Aligment file in [XMAP,OMD] format")
  # parser.add_argument('-p', '--per_align', type=float, default=30, help="For filtering: minimum percentage aligned.")
  parser.add_argument("-p", '--per_align', nargs=3, metavar=('min', 'max', 'step'), type=int, default=(0, 50, 10), help="For filtering: minimum percentage aligned, multi-settings [min, max, step]")
  parser.add_argument('-c', '--confidence', nargs=3, metavar=('min', 'max', 'step'), type=int, default=(0, 20, 5), help="For filtering: minimum confidence score of alignment, multi-settings [min, max, step].")
  # parser.add_argument('-c', '--confidence', type=float, default=13, help="For filtering: minimum confidence score of alignment.")
  # parser.add_argument('-s', '--resolution', type=int, default=20, help="resolution (default:40) of OM red labels.")
  parser.add_argument('-o', '--output', default="align_stats.csv", help='output file')

  args = parser.parse_args()
  return args

args = getArguments()

print("[Stage] Load reference info..")
with open(args.ref, 'r') as fin:
    line2skip = 0
    line = fin.readline().strip(" \r\n")
    while line and line[0] == '#':
        line2skip = line2skip + 1
        if len(line) > 1 and line[1] == 'h':
            colnames = line[3:].split("\t")
        line = fin.readline().strip(" \r\n")
    # print(line2skip, "lines in header to skip.")
    # print("#Colnames:", ",".join(colnames), flush=True)

ref_data = pd.read_csv(args.ref, sep='\t', skiprows=line2skip, names=colnames)
cols=[i for i in ref_data.columns if i not in ["CMapId"]]
for col in cols:
    ref_data[col]=pd.to_numeric(ref_data[col])
# ref_data_simple = ref_data.apply(pd.to_numeric)
ref_data_simple = ref_data.drop_duplicates("CMapId")

# for index, row in ref_data_simple.iterrows():
#     ref_len = ref_len + row['RefEndPos'] - row['RefStartPos']
reflentotal = sum(ref_data_simple['ContigLength'])
print(reflentotal, "bp in reference sequences.")

print("[Stage] Load query info..")
with open(args.query, 'r') as fin:
    line2skip = 0
    line = fin.readline().strip(" \r\n")
    while line and line[0] == '#':
        line2skip = line2skip + 1
        if len(line) > 1 and line[1] == 'h': # maybe no h
            colnames = line[3:].split("\t")
        line = fin.readline().strip(" \r\n")
    # print(line2skip, "lines in header to skip.")
    # print("#Colnames:", ",".join(colnames))

query_data = pd.read_csv(args.query, sep='\t', skiprows=line2skip, names=colnames)
cols=[i for i in query_data.columns if i not in ["CMapId"]]
for col in cols:
    query_data[col]=pd.to_numeric(query_data[col])
# query_data_simple = query_data.apply(pd.to_numeric).drop_duplicates("CMapId")
query_data_simple = query_data.drop_duplicates("CMapId")
print(query_data_simple.shape[0], "molecules to align.")

# determine the format of the alignment file
format={"XMAP":{},
        "OMD": {"RefID":"RefContigID", "RefStartCoord":"RefStartPos", "RefStopCoord":"RefEndPos",
                # "QryStartPos", "QryEndPos",
                "QuerySize":"QryLen", "Confidence":"Confidence",
                "AlignedSegRatio":"PerAlign"}}

print("[Stage] Load alignment info..")
alignFileFormat = args.align.split(".")[-1].upper()
with open(args.align, 'r') as fin:
    line2skip = 0
    line = fin.readline().strip(" \r\n")
    while line and line[0] == '#':
        line2skip = line2skip + 1
        if alignFileFormat == "XMAP" and len(line) > 1 and line[1] == 'h': # the header line of XMAP starting with "#h "
            colnames = line[3:].split("\t")
        if alignFileFormat == "OMD": # the last line starting with # shows the column names
            colnames = line[1:].split("\t")
        line = fin.readline().strip(" \r\n")
    # print(line2skip, "lines in header to skip.")
    # print("#Colnames:", ",".join(colnames))

align_data = pd.read_csv(args.align, sep='\t', skiprows=line2skip, names=colnames)
if alignFileFormat == "XMAP":
    # numericCols = ["RefContigID",
    numericCols = ["RefStartPos", "RefEndPos", "QryStartPos", "QryEndPos", "QryLen", "Confidence"]
    align_data[numericCols] = align_data[numericCols].apply(pd.to_numeric)
    align_data['PerAlign'] = abs(align_data['QryEndPos'] - align_data['QryStartPos'])*100/align_data['QryLen']
elif alignFileFormat == "OMD":
    align_data.rename(columns = format[alignFileFormat], inplace=True)
    numericCols = list(format[alignFileFormat].values())
    align_data[numericCols] = align_data[numericCols].apply(pd.to_numeric)
    align_data['PerAlign'] = align_data['PerAlign'] * 100
    align_data['Confidence'] = align_data['Confidence'] * 100
else:
    print("Unrecognized align file format.")
print("PerAlign: min. ", min(align_data['PerAlign']), ", max. ", max(align_data['PerAlign']))
print("Confidence: min. ", min(align_data['Confidence']), ", max. ", max(align_data['Confidence']))

result = pd.DataFrame(columns=["Percent aligned", "Confidence", "Coverage"])
for per_align in range(args.per_align[0], args.per_align[1]+args.per_align[2], args.per_align[2]):
    per_align = min(per_align, args.per_align[1])

    for confidence in range(args.confidence[0], args.confidence[1]+args.confidence[2], args.confidence[2]):
        confidence = min(confidence, args.confidence[1])

        align_data_filter = align_data[(align_data["Confidence"] >= confidence) & (align_data['PerAlign'] >= per_align)]
        align_data_filter.reset_index(inplace=True)
        align_data_sorted = align_data_filter.sort_values(by=["RefContigID", "RefStartPos", "RefEndPos"])

        notoverlap_per=0
        if align_data_sorted.shape[0] > 0:
            notoverlap_aligns = pd.DataFrame(columns=["RefContigID", "RefStartPos", "RefEndPos"])
            refid = align_data_sorted['RefContigID'][0]
            refstart = align_data_sorted['RefStartPos'][0]
            refend = align_data_sorted['RefEndPos'][0]
            idx = 1
            # print(align_data_sorted[["RefContigID", "RefStartPos", "RefEndPos"]].head(10))
            while idx < align_data_sorted.shape[0]:
                if align_data_sorted['RefContigID'][idx] == refid and refend >= align_data_sorted['RefStartPos'][idx]:
                    refend = max(refend, align_data_sorted['RefEndPos'][idx])
                else:
                    notoverlap_aligns.loc[len(notoverlap_aligns.index)] = [refid, refstart, refend]
                    refid = align_data_sorted['RefContigID'][idx]
                    refstart = align_data_sorted['RefStartPos'][idx]
                    refend = align_data_sorted['RefEndPos'][idx]
                idx = idx + 1
            notoverlap_aligns.loc[len(notoverlap_aligns.index)] = [refid, refstart, refend]
            notoverlap_aligns["Len"] = notoverlap_aligns["RefEndPos"] - notoverlap_aligns["RefStartPos"]
            notoverlap_lentotal = sum(notoverlap_aligns["Len"])
            notoverlap_per = round(notoverlap_lentotal*100/reflentotal, 2)
        result.loc[len(result.index)] = [per_align, confidence, notoverlap_per]
        # print(notoverlap_lentotal, "bp covered by reads.")
        # print(notoverlap_per, "% reference covered by alignment.")

result.to_csv(args.output, index=False, sep=',')
# print(notoverlap_aligns[["RefContigID", "Len"]].groupby("RefContigID").sum())
# notoverlap_aligns.to_csv("test.tsv", sep='\t')
# print(notoverlap_aligns.head(10))
# print(notoverlap_aligns.tail(10))

# print(align_data_sorted.loc[align_data_sorted['RefContigID'] == 2][["RefContigID", "RefStartPos", "RefEndPos"]].head(10))

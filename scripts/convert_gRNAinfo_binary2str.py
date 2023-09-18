import sys
import argparse

def getArguments():
  parser = argparse.ArgumentParser(description="Convert gRNAinfo from binary to string.")
  parser.add_argument('-i', '--input', required=True, help='The gRNA info file in tsv format')
  parser.add_argument('-o', '--output', required=True, help='Output file')

  args = parser.parse_args()
  return args

args = getArguments()

DNA_BASES = ['A', 'C', 'G', 'T']
DNA_BASES_to_Binary = {'A':0b00, 'C':0b01, 'G':0b10, 'T':0b11}

# not support N
def gRNA_bin2str(kmer_b, gRNA_len=20):
  kmer = ""
  for x in range(0, gRNA_len):
    idx = (kmer_b & 0b11)
    kmer = DNA_BASES[idx] + kmer
    kmer_b = (kmer_b >> 2)
  return kmer

def DNA2Binary(seq):
  kmer = 0
  for x in range(len(seq)):
    kmer = (kmer | DNA_BASES_to_Binary[seq[x].upper()])
    kmer = (kmer << 2)
  kmer = (kmer >> 2)
  return kmer

with open(args.input, 'r') as fin, open(args.output, 'w') as fout:
  # the first line is default to the header line
  line = fin.readline()
  fout.write(line)
  # Loop through all gRNA lines, the first element should be the gRNA sequence
  line = fin.readline().strip(" \r\n")
  while line:
    eles = line.split("\t")
    eles[0] = gRNA_bin2str(int(eles[0]))
    fout.write("\t".join(eles)+"\n")
    line = fin.readline().strip(" \r\n")

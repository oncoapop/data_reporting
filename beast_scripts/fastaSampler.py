#!/usr/bin/env python
import random
import sys

#name of the input file (fasta or fastq)
#assumes input file is standard fasta/fastq format
fileName = sys.argv[1]
#number of sequences to subsample
numSeq = int(sys.argv[2])
increment = 0

#if it's a fasta file
if (fileName.find(".fasta") != -1):
  increment = 2
#else if it's a fastq file
elif (fileName.find(".fastq") != -1):
  increment = 4
#quit if neither
else:
  sys.stdout.write("not a fasta/fastq file\n")
  sys.exit()

FILE = open(fileName, 'r')
reads = list()
index = 0
#put the entire file in a list
for line in FILE:
  if(index % increment == 0):
    reads.append(line)
  else:
    reads[(index/increment)] += line
  index += 1
FILE.close()

#figure out total number of reads, error if applicable
totalReads = index/increment
if(totalReads < numSeq):
  sys.stdout.write("You only have "+str(totalReads)+" reads!\n")
  sys.exit()

#shuffle the reads!
random.shuffle(reads)
#output the reads to stdout
for i in range(0, numSeq):
  sys.stdout.write(reads[i])

#!/usr/bin/env python

import random
import glob
inFiles = glob.glob('*.fasta') + glob.glob('*.fastq')
outFiles = []
num = int(raw_input("Enter number of random sequences to select:\n"))

for i in range(len(inFiles)):

        for k in range(3):
                fNames = []
                fSeqs = []
                same = 0
                outFiles.append(file(str(inFiles[i])+'_'+'Rand_'+str(num)+'-'+str(k+1)+'.fa', 'wt'))
                for line in open(inFiles[i]):
                        if ((line[0] == '>') or (line[0] == '@')):
                                fNames.append(line)
                                same = 0
                        else:
                                if (same != 1):
                                        fSeqs.append(line)
                                        same = 1
                                else:
                                        fSeqs[(len(fSeqs)-1)] += line
                curr = (len(outFiles)-1)
                for j in range(num):
                        a = random.randint(0, (len(fNames)-1))
                        outFiles[curr].write(fNames.pop(a))
                        outFiles[curr].write(fSeqs.pop(a))
raw_input("Done.")


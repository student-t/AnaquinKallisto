#!/usr/bin/python

#
# Simulate on sequins given a FASTA file and mixture file.
#

import os
import sys
import math
import subprocess
from random import randint

# Split a file of sequin into individual sequins
def split(file, seq_path):
    os.system('mkdir -p ' + seq_path)

    with open(file) as f:
        while True:
            l1 = f.readline()
            l2 = f.readline()            
            if (not l2):
                break

            file = l1.replace(">", "")
            file = file.replace("?", "")
            file = file.replace("\n", "")
    
            w = open(seq_path + os.sep + file + '.fa', 'w')
            w.write(l1)
            w.write(l2)

def readMix(file):
    r = {}

    with open(file) as f:
        l = f.readline()
        split = ',' if ',' in l else '\t'

        while True:
            l = f.readline()
            if (not l):
                break

            toks = l.strip().split(split)

            if (toks[0] == 'id'):
                continue

            r[toks[0]] = { 'id': toks[0], 'A': float(toks[1]) }                           
                           
    return r

# Generate simulated reads for each sequin for a given mixture
def simulate(file, scale, basePath, mix, min_=0, max_=sys.maxsize, isTooLowError=True):
    mixFile = readMix(file)
    covs = []

    for f in os.listdir(basePath):
        key = f.split('.')[0]

        if key in mixFile:            
            # Length of the sequin
            length = 1000

            # The reads depend on the concentration level and sequin length
            reads = scale * mixFile[key][mix]

            #print '\n------------------ ' + key + ' ------------------'
            
            # Path for the sequin
            path  = basePath

            # This is the number of reads that we'll need
            reads = int(reads)

            #
            # The number of reads need to be adjusted for the sequin length.
            # The implementation borrows from Wendy's script. For example,
            #
            #   - length is 1689
            #   - reads is 468750
            #
            # We would calculate: 468750 * (1689 / 1000) to get reads per KB.
            #
            
            reads = max(min_, reads)
            reads = min(max_, reads)
            reads = reads * (length / 1000)

            # Don't bother if the abundance is too low or too high
            if (reads > 1):
                #print 'Generating: ' + str(reads)

                # Simulate reads from a given sequin
                i  = path + '/' + key + '.fa'
                o1 = path + '/' + key + '.R1.fq'
                o2 = path + '/' + key + '.R2.fq'

                #cmd = 'wgsim -1150 -2150 -d 130 -S ' + str(randint(1,100)) + ' -N ' + str(reads) + ' ' + i + ' ' + o1 + ' ' + o2 + ' > /dev/null'
                cmd = 'wgsim -r 0 -R 0 -X 0 -e 0 -1150 -2150 -d 130 -S ' + str(randint(1,100)) + ' -N ' + str(reads) + ' ' + i + ' ' + o1 + ' ' + o2 + ' > /dev/null'
                                
                # What's the estimated coverage?
                cov = (150 * 2 * reads) / length
                
                covs.append(cov)
                
                #print(cmd)
                os.system(cmd)
                
                if (os.stat(o1).st_size == 0):
                    raise Exception(key + ' not generated!')
                if (os.stat(o2).st_size == 0):
                    raise Exception(key + ' not generated!')                    
            else:
                if isTooLowError:
                    raise Exception('Error: ' + key + ' not generated! Reads: ' + str(reads))
                else:
                    print('Warning: ' + key + ' not generated!')
        else:
            print('-------- Warning --------: ' + key + ' not found in the mixture!')

    os.system('cat ' + basePath + '*R1.fq > ' + basePath + 'simulated_1.fastq')
    os.system('cat ' + basePath + '*R2.fq > ' + basePath + 'simulated_2.fastq')

def print_usage():
    print('Usage: python3 simulate.py sequins.fa sequins_allele_freqs.tsv 10000')

if __name__ == '__main__':
    if (len(sys.argv) != 4):
        print_usage()
    else:
        split(sys.argv[1], 'Simulation/')
        simulate(sys.argv[2], int(sys.argv[3]), 'Simulation/', 'A')

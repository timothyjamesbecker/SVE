#!/usr/bin/env python
import sys
import numpy as np
import argparse
des = """Random SAM/BAM/CRAM sampler for piped data processing v0.0.1
Timothy James Becker, PhD candidate, UCONN 04/27/2017
------------------------------------------------------------
Given a BAM input, pass the header while randomly sampling with bernoulli probabilty -s
[USAGE]
samtools view -Sh sample.bam | random_sample.py -s 1E-3"""
parser = argparse.ArgumentParser(description=des,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-s', '--sample_prob',type=float,help='sampling probabilty\t[1E-6]')
args = parser.parse_args()
if args.sample_prob is not None:
    s = args.sample_prob
    if s <= 0.0: s = float(1E-6) #clip
    elif s > 1.0: s *= (1.0/s)   #scale
else: s = float(1E-6)
for line in sys.stdin:
    if line.startswith('@'):
        print(line)
    elif np.random.choice([True,False],p=[s,1.0-s]):
        print(line)
sys.exit()

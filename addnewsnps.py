#! /usr/bin/env python
# Input:
#       <frequency file, 4 columns: (chrom, pos, allele1:freq1, allel2:freq2) e.g. (chr1 23000000 A:0.1 B:0.9)>
#           This could most likely come from either PLINK or VCFTOOLS.
#       <new positions (chr pos)> i.e. output from 'makenewpos.py'
#       <output>
#       <searchwindow [min,]max> e.g. 100,1000
# Script will search within a window for the highest MAF SNP, or for the first SNP outside of that window

import sys

freqfile = sys.argv[1]
newposfile = sys.argv[2]
outfile = sys.argv[3]
w = sys.argv[4].split(',')

if len(w) == 1:
    minwindow,maxwindow = 0,int(w[0])
else:
    minwindow,maxwindow = int(w[0]),int(w[1])

# Read the estimated positions we want to fill to have a uniform SNP coverage of the genome.
fill = []
with open(newposfile,'r') as fin:
    for line in fin:
        l = line.strip().split()
        if len(l) < 1:
            continue
        fill.append([l[0],l[1]])

# Read SNP information for the SNPs we can use as potential fillers.
freqall = {}
with open(freqfile,'r') as fin:
    header = next(fin)
    for line in fin:
        l = line.strip().split()
        if len(l) < 1:
            continue
        chrom,pos,a1,a2 = l[0],int(l[1]),float(l[4].split(':')[1]),float(l[5].split(':')[1])
        maf = min(a1,a2)
        if chrom not in freqall:
            freqall[chrom] = []
        freqall[chrom].append([chrom,pos,maf])

# Deciding which SNP to select for each of the "empty holes"
freq = []
with open(outfile,'w') as fout:
    while len(fill) > 0:
        # Pop out the first position to fill
        chrom1,pos1 = fill.pop(0)
        pos1 = int(pos1)
        if chrom1 not in freqall:
            fout.write('{}\t{}\tNAN\n'.format(chrom1,pos1))
            continue
        freq = freqall[chrom1]
        # Locate the two positions closest to pos1
        freqind = 0
        while freqind < len(freq) and freq[freqind][0] == chrom1 and freq[freqind][1] < pos1:
            freqind += 1
        up = freqind - 1
        if freqind < len(freq): down = freqind
        else: down = freqind - 1
        if (pos1 - freq[up][1]) < (freq[down][1] - pos1): # Up is closer
            collect = freq[up]+[up]
        else:
            collect = freq[down]+[down]
        # Search up (backwards) until outside window
        while up >= 0 and freq[up][1] > (pos1 - maxwindow):
            if freq[up][2] > collect[2] and (freq[up][1] < (pos1-minwindow)):
                collect = freq[up]+[up]
            up -= 1
        # Search down (forwards) until outside window
        while down < len(freq) and freq[down][0] == chrom1 and freq[down][1] < (pos1 + maxwindow):
            if freq[down][2] > collect[2] and (freq[down][1] > (pos1+minwindow)):
                collect = freq[down]+[down]
            down += 1
        if collect[0] == '0':
            # If down is still within current chromosome, compare
            if down < len(freq) and freq[down][0] == chrom1 and freq[down][2] > freq[up][2]:
                collect = freq[down]+[down]
            else:
                collect = freq[up]+[up]
        if collect[0] == chrom1:
            fout.write('{}\t{}\t{:.3f}\n'.format(collect[0],collect[1],collect[2]))
            freq.pop(collect[3])
        else: fout.write('{}\t{}\tNAN\n'.format(chrom1,pos1))

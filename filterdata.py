#! /usr/bin/env python
# Removes SNP from a map-file so that no SNP is within a given distance from each other.
# All SNPs must be sorted according to chromosome and position
# Input: One file containing SNP information.
#        First line is a header line that will be copied directly to the output
#        One SNP pr. line and the second column must contain the physical position for the SNP
# Output: One file with same format as input file, containing only the included lines.
import sys

snpdist = int(sys.argv[3])
inputfile = sys.argv[1]
outputfile = sys.argv[2]

line1ist = []
pos1ist = []
with open(inputfile,'r') as fin, open(outputfile,'w') as fout:
    fout.write(next(fin))
    # Read all inputlines and the position of each SNP.
    for line in fin:
        l = line.strip().split()
        if len(l) < 1:
            continue
        try:
            pos = int(l[1])
        except:
            continue
        line1ist.append(line)
        pos1ist.append(pos)
    # Remove all lines closer than the given snpdist to both neighbors
    # Skips first and last SNP pr. chromosome, assumes that the bp-value of last SNP on one chromosome
    # is higher than bp-value of first SNP on next chromosome.
    i = 1
    while i < len(line1ist)-1:
        p1,p2,p3 = pos1ist[i-1],pos1ist[i],pos1ist[i+1]
        if p3 < p2 or p1 > p2:
            i += 1
            continue
        if abs(p2-p1) < snpdist and abs(p3-p2) < snpdist:
            line1ist.pop(i)
            pos1ist.pop(i)
            continue
        i += 1
    # Remove all lines closer than dist to the previous line
    i = 1
    while i < len(line1ist)-1:
        p1,p2 = pos1ist[i-1],pos1ist[i]
        if p1 > p2:
            i += 1
            continue
        if abs(p2-p1) < snpdist:
            line1ist.pop(i)
            pos1ist.pop(i)
            continue
        i += 1
    for el in line1ist:
        fout.write(el)

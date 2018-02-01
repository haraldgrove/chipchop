#! /usr/bin/env python
# The program requires 4 inputs on the command line with the order given below.
# Input: <input file: SNP list with at least two columns: chromosome and bp-position, sorted>
#        <output file>
#        <position of the 2 columns in input file with: chromosome and basepair position> e.g. 1,2
#        <inter-SNP distance)
# Output: File with chromosome and bp-position where a SNP should be inserted
#         The file contain one special line at the end of each chromosome, where the word 'END' is added to the
#         bp-position of the last SNP. (e.g. chr1 END123000000)
import sys



info = sys.argv[3].split(',')
bpcol = int(info[1])-1
chromcol = int(info[0])-1
snpdist = int(sys.argv[4])
inputfile = sys.argv[1]
outputfile = sys.argv[2]

prevpos = 0
prevchrom = '0'
with open(inputfile,'r') as fin, open(outputfile,'w') as fout:
    header = next(fin)
    for line in fin:
        l = line.strip().split()
        if len(l) < 1:
            continue
        pos = int(l[bpcol])
        if pos < prevpos: # First position in new chromosome
            fout.write('{}\t{}{}\n'.format(prevchrom,'END',prevpos))
            prevpos = 0
            prevchrom = l[chromcol]
        # Calculate number of SNPs needed to fill the gap between the previous SNP and current,
        # given the wanted inter-SNP distance.
        newsnps = ((pos-prevpos)//snpdist)
        if newsnps > 1:
            locdist = (pos-prevpos)//(newsnps)
            nextpos = prevpos+locdist
            for i in range(0,newsnps-1):
                fout.write('{}\t{}\n'.format(l[chromcol],nextpos))
                nextpos += locdist
        prevpos = pos
        prevchrom = l[chromcol]

# analysis of G4 related sequence
# python=3.9
# Bin Zhang
# Data: May 29, 2022
# version: 1.0
# ENV: 

import re
import sys

# calcaulte [ACGT] content of the giving sequence file: column 1 is id and column2 is the sequence 
def calNContent (in_txt,out_file,nuc):
    with open (out_file,'w') as w:
        n = 0
        with open (in_txt,'r') as r:
            for _line in r:
                n += 1
                if n % 10000 == 0:
                    print('### prcocessed {} sequence'.format(n))
                eles = re.split('[\t\ ]',_line.rstrip())
                seq = eles[1].upper()
                subs = seq.split(nuc)
                prop = float(len(subs) - 1)/len(seq)
                #print(seq,len(subs))
                w.write(eles[0] + '\t' + str(prop) + '\t' + str(len(seq)) + '\n')
    print('### calculating content of {} Done'.format(nuc))
    
def G4predict (in_txt,out_file):
    with open (out_file,'w') as w:
        n = 0
        with open (in_txt,'r') as r:
            for _line in r:
                n += 1
                if n % 10000 == 0:
                    print('### prcocessed {} sequence'.format(n))
                eles = re.split('[\t\ ]',_line.rstrip())
                seq = eles[1].upper()
                p1 = 'NA'
                p2 = 'NA'
                p3 = 'NA'
                p4 = 'NA'
                # pqs1: 4G
                pqsf1 = re.search(r'G{3,}(.{1,7}?G{3,}){3,}',seq,re.I)
                if pqsf1:
                    p1 = pqsf1.span()
                # pqs2: 4GL15
                for n in range(0,3):
                    pattern = '(G{3,}.{1,7}?){' + str(n) + ',}G{3,}.{1,15}?G{3,}(.{1,7}?G{3,}){' + str(2-n) + ',}'
                    pqsf2 = re.search(pattern,seq,re.I)
                    if pqsf2:
                        p2 = pqsf2.span()
                        break
                # pqs3: Bugle
                for n in range(0,4):
                    pattern = '(G{3,}.{1,7}?){' + str(n) + ',}(GG[^G]GG?|GG?[^G]GG)(.{1,7}?G{3,}){' + str(3-n) + ',}'
                    pqsf3 = re.search(pattern,seq,re.I)
                    if pqsf3:
                        p3 = pqsf3.span()
                        break
                # pqs4: GVBQ
                for n in range(0,4):
                    pattern = '(G{3,}.{1,7}?){' + str(n) + ',}G{2}(.{1,7}? G{3,}){' + str(3-n) + ',}'
                    pqsf4 = re.search(pattern,seq,re.I)
                    if pqsf4:
                        p4 = pqsf4.span()
                        break
                w.write('\t'.join(str(e) for e in [eles[0],p1,p2,p3,p4]) + '\n')
    print('### Predicting G4 Done')
                                

if len(sys.argv) < 3:
    print("Usage:python ana.g4.seq.py method paramters\nmethods should be nucleotide [content] or [G4] prediction")
    print('python ana.g4.seq.py content input_txt output AorGorCorT')
    print('python ana.g4.seq.py G4 input_txt output')
    sys.exit(1)
else:
    if 'content' in sys.argv[1]:
        print('### Calculating nulceotide content') 
        calNContent(sys.argv[2],sys.argv[3],sys.argv[4])
    elif 'G4' in sys.argv[1]:
        print('### Predicting G4 based on motif')
        G4predict(sys.argv[2],sys.argv[3])
                
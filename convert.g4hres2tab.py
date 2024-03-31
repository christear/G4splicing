# convert the results from G4Hunter to table format 
# python=3.9
# Bin Zhang
# Data: Aug 4, 2022
# version: 1.0
# ENV: 

import re
import sys

# the results of G4Hunter in a format: >id\n head \n 5_columns_tab \n 
# extract the id, position from the 5_columns_tab and fill missing values with NA
def g4hres2tab (g4h_s,out_file):
    loci_win = {}
    score_dict = {}
    _tmp_id = ''
    _i = 0
    with open(g4h_s,'r') as r:
        for _line in r:
            _line = _line.rstrip()
            if '>' in _line:
                _tmp_id = _line[1:]
                score_dict[_tmp_id] = {}
            else:
                eles = re.split('[\t\ ]+',_line)
                #print(eles)
                if len(eles) == 5:
                    #break
                    loci_win[eles[0]] = 1
                    score_dict[_tmp_id][eles[0]] = float(eles[4])
                else:
                    #_i += 1
                    continue
                    print(_tmp_id)
                    #if _i > 10:
                    #    break
    #print(score_dict)
    #print(loci_win)
    with open (out_file,'w') as w:
        _head = 'id\t' + '\t'.join(loci_win.keys())
        w.write(_head + '\n')
        #print(_head)
        for _id in score_dict:
            out_col = []
            for _loci in loci_win:
                if _loci in score_dict[_id].keys():
                    out_col.append(score_dict[_id][_loci])
                else:
                    out_col.append('NA')
            _out_line = _id + '\t' + '\t'.join(str(e) for e in out_col)
            #print(_out_line)
            w.write(_out_line + '\n')
                

if len(sys.argv) < 2:
    print("Usage:python covert.g4hres2tab.py input_from_g4h output_file")
    sys.exit(1)
else:
    g4hres2tab(sys.argv[1],sys.argv[2])


                
        

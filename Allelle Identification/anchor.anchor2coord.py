#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#@File : anchor.anchor2coord.py
#@Author: Qionghou Li

import sys
import numpy as np

with open(sys.argv[1],'r') as op1,open(sys.argv[2],'w') as op2:
    for i in op1:
        if i.startswith(('#','refChr')):
            continue
        else:
            line=i.strip().split('\t')
            len1=str(int(line[2])-int(line[1])) if str(int(line[2])-int(line[1])) != "-1" else "0"
            len2=str(int(line[5])-int(line[4])) if str(int(line[5])-int(line[4])) != "-1" else "0"
            if len1==len2=="-1":
                continue
            else:
                if line[6] == "+":
                    wrt=[str(int(line[1])),line[2],str(int(line[4])-1),line[5],len1,len2,"0","1","1",line[0],line[3]]
                elif line[6] == "-":
                    wrt=[str(int(line[1])),line[2],str(int(line[4])-1),line[5],len1,len2,"0","1","-1",line[0],line[3]]
                else:
                    raise TypeError
                op2.write('\t'.join(wrt)+'\n')

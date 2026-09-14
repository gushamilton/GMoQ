#!/bin/bash
# list every completed raw dir
B=/user/work/fh6520/kde_lumpy
ls -d $B/run[ABCD]/*/ 2>/dev/null | sed "s:/$::" > $B/joblistS2.txt
wc -l < $B/joblistS2.txt

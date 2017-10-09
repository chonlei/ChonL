import numpy as np
import glob
import re

PATH_TO_DIR = "../out_PaciSweep_failed_re2/"
FILE_NAMES = "outputPaciSimulation*.txt"

all_existing_files = glob.glob(PATH_TO_DIR+FILE_NAMES)
list_idx = []

for existing_file in all_existing_files:
    idx = int(re.findall("\.\.\/out\_PaciSweep\_failed\_re2\/outputPaciSimulation(\d+).+",existing_file)[0])
    list_idx.append(idx)

if len(list_idx) != len(all_existing_files):
    print "WARNING: Not extracting all existing files"

list_idx.sort()


readfile = open('conductanceDataTest_out_PaciSweep_failed_re2.txt', 'r')
writefile = open('conductanceDataTest.txt', 'w')

for line in readfile:
    sline = line.split()
    if int(sline[0]) not in list_idx:
        writefile.write("%d"%int(sline[0]))
        for a in sline[1:]:
            writefile.write("\t%s"%a)
        writefile.write("\n")

readfile.close()
writefile.close()

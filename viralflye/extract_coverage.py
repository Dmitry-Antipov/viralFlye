############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2015 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from __future__ import with_statement
import os
import shlex
import shutil
import stat
import subprocess
import sys
import platform
import re
import gzip
import time
import random
from collections import defaultdict
from genericpath import isdir, exists

from os.path import join

#Do not use sra tools from chihua! it's outdated and bugged
sra_tools = "/home/dantipov/other_tools/sratoolkit.2.10.7-ubuntu64/bin/"



def main():
    file = sys.argv[1]
    stats = sys.argv[2]
    covs = {}
    alts = {}
    for line in open (stats,'r'):
        arr = line.split()
        if arr[2] == 'cov.':
            continue
        covs[arr[0]] = float(arr[2])
        alts[arr[0]] = arr[-2]
    count = 0
    alt_group = 0
    total_cov = 0
    for line in open (file,'r'):
        arr = line.split(',')
        if arr[3] == 'Partial':
            count +=1
            if alts[arr[0]] != '*':
                alt_group +=1
            total_cov += covs[arr[0]]
    print (f'{count} {total_cov/count} {alt_group}') 

if __name__ == '__main__':
    main()




# Kill a list of processes specified in a file ###

import sys
import subprocess as sp

with open(sys.argv[1]) as f:
    for line in f:
        fields = str.strip(line).split()
        killcmd = 'kill ' + fields[1]
        sp.call(killcmd, shell = True)

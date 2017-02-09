import sys
import subprocess as sp

sourceDirectory = sys.argv[1]
targetDirectory = sys.argv[2]

cmd1 = 'cp $(find ' + sourceDirectory  + ' -type f -name "*.py") ' + targetDirectory
cmd2 = 'cp $(find ' + sourceDirectory  + ' -type f -name "*.sh") ' + targetDirectory

sp.call(cmd1,shell=True)
sp.call(cmd2,shell=True)

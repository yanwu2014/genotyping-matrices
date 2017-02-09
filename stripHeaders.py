import sys
import os

filesDir = sys.argv[1]
fileType = sys.argv[2]

files = os.listdir(filesDir)
files = [x for x in files if x.endswith(fileType)]

print 'Stripping Headers from following files: '
for fi in files:
    print fi

os.chdir(filesDir)

for fi in files:
    lines = []
    with open(fi) as f:
        header = f.readline()
        for line in f:
            lines.append(line)
    with open(fi,'w') as f:
        for line in lines:
            f.write(line)

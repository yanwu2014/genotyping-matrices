### Script to combine fastq files from multiple lanes into a single fastq file ### 

import sys
import subprocess as sp

def main():
    sampleName = sys.argv[1]
    sample_ids = []
    for i in range(2,len(sys.argv)):
        sample_ids.append(sys.argv[i])
    
    combined_R1_file = sampleName + '_R1'
    combined_R2_file = sampleName + '_R2'

    # generate read 1
    cat_files(sample_ids, '/*_R1_*', combined_R1_file)
    
    # generate read 2
    cat_files(sample_ids, '/*_R2_*', combined_R2_file)

def cat_files(sample_ids, wildcard, out_prefix):
    in_files = []
    for sample_dir in sample_ids:
        cmd = 'ls ' + sample_dir + wildcard
        in_paths = sp.check_output(cmd, shell = True)
        in_paths = str.strip(in_paths).split('\n')

        in_files.extend(in_paths)

    cat_str = ''
    for fi in in_files:
        cat_str = fi + ' ' + cat_str
    cat_str = 'cat ' + cat_str + '> ' + out_prefix + '.fastq.gz'
    
    sp.call(cat_str, shell = True)

        
if __name__ == '__main__': main()

"""Index HapMap LD files for fast lookups

Output (chromosome, position, file offset) pairs for the forward and reverse
directions.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import os

def index(filename, key, chromosome, outfile):
    pos = 0
    last = ''
    with open(filename, 'rb') as f:
        line = f.readline()
        while line:
            k = key(line)
            if k != last:
                last = k
                pos = f.tell() - len(line)
                print(chromosome, k, pos, sep=',', file=outfile)
            line = f.readline()

if __name__ == '__main__':
    with open('forward_index.csv', 'w') as f:
        for i in range(1, 23):
            index(os.path.expanduser('~/hp/hapmap/ld_chr{}_CEU.txt').format(i),
                  lambda s: int(s.split()[0]), i, f)
    with open('reverse_index.csv', 'w') as f:
        for i in range(1, 23):
            index(os.path.expanduser('~/hp/hapmap/ld_chr{}_rev.txt').format(i),
                  lambda s: int(s.split()[0]), i, f)

"""Filter PLINK dosage files

Usage: python filter.py RSIDS

Expects dosage on stdin. Writes dosage on stdout.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import sys
import argparse

p = argparse.ArgumentParser()
p.add_argument('-v', '--invert-match', help='Invert match', action='store_true', dest='invert')
p.add_argument('rsids', nargs=1)
args = p.parse_args()

with open(args.rsids[0]) as f:
    rsids = set(line.strip() for line in f)

if args.invert:
    for line in sys.stdin:
        if line.split()[0] not in rsids:
            print(line, end='')
else:
    for line in sys.stdin:
        if line.split()[0] in rsids:
            print(line, end='')

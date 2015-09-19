import signal
import sys

def parse(row):
    for link in row[11].split(';'):
        if link:
            eid, _ = link.split('_')
            yield row[1:4] + [row[5], eid]

signal.signal(signal.SIGPIPE, signal.SIG_DFL)
data = (line.split() for line in sys.stdin)
for row in data:
    for link in parse(row):
        print(*link, sep='\t')

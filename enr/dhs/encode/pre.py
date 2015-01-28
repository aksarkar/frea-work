import re
import sys

field = re.compile(r'([^= ]*)=([^;]*);')
for line in sys.stdin:
    filename, rest = line.split('\t')
    fields = dict(field.findall(rest))
    if 'treatment' in fields and fields['treatment'] != 'None':
        print(filename, '{cell}-{treatment}'.format(**fields))
    else:
        print(filename, '{cell}'.format(**fields))

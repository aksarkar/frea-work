import sys

print('ID_1 ID_2 missing gender cohort supplier well region ethnicity age_recruitment age_onset case')
print('0 0 0 D D D D D D D D B')

data = (line.split() for line in sys.stdin)
next(data)
for row in data:
    print(row[0], row[0], 0, *row[1:])

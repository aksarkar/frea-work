import sys

with open(sys.argv[1]) as f:
    data = (line.split() for line in f)
    test = next(data)
    n_eigenvectors = len(test) - 2
    eigenvectors = {(row[0], row[1]): [float(x) for x in row[2:]] for row in data}
    eigenvectors[(test[0], test[1])] = [float(x) for x in test[2:]]

with open(sys.argv[2]) as f:
    header1 = next(f).split()
    header1[-1] = 't1d'
    print(' '.join(header1), *['pc{}'.format(i) for i in range(n_eigenvectors)])
    print(next(f).strip(), *['C' for _ in range(n_eigenvectors)])
    data = (line.split() for line in f)
    for row in data:
        row[6] = str(int(row[6]) - 1)  # fix case-control status
        print(' '.join(row), *eigenvectors[(row[0], row[1])])
        

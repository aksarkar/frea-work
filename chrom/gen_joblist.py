if __name__ == '__main__':
    candidates = [[4, 5, 11, 12], [11, 12], [4, 5], [11], [12], [4], [5]]
    cell_types = ['GM12878','H1','HMEC','HSMM','HepG2','Huvec','K562','NHEK','NHLF']
    for t in cell_types:
        for c in candidates:
            print('{},"{}"'.format(t, ','.join(str(x) for x in c)))

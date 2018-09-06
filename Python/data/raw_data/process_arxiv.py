mapping = {}
ind = 1
with open('arxiv_compressed.txt', 'wb') as o:
    with open('arxiv.txt', 'rb') as f:
        for l in f:
            l = l.split()
            n1, n2 = l[0], l[1]
            for a in [n1, n2]:
                if not a in mapping:
                    mapping[a] = str(ind)
                    ind += 1
            o.write(mapping[n1]+' '+mapping[n2]+'\n')


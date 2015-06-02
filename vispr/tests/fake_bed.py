from __future__ import print_function, division, absolute_import

import sys, csv, random


reader = csv.reader(sys.stdin)
gene = None
next(reader)
for l in reader:
    name = l[0]
    if l[1] != gene:
        gene = l[1]
        gene_start = random.randrange(1, 240000000)
    start = gene_start + random.randrange(0, 10000)
    score = random.uniform(0, 1)
    print("1", start, start + 20, name, score, sep="\t")

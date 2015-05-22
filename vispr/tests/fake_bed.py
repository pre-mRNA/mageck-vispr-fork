import sys, csv, random

reader = csv.reader(sys.stdin)
next(reader)
for l in reader:
    name = l[0]
    start = random.randrange(10000, 16000)
    score = random.uniform(0, 1)
    print("1", start, start + 20, name, score, sep="\t")

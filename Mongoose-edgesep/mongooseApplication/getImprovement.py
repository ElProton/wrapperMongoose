import csv
import sys

reader = csv.reader(open(sys.argv[1], 'r'))
basedic = {}
head = next(reader)
head = head[1:]
for row in reader:
    k = row[0]
    v = row[1:]
    basedic[k] = v


reader = csv.reader(open(sys.argv[2], 'r'))
moduledic = {}
next(reader)
for row in reader:
    k = row[0]
    v = row[1:]
    moduledic[k] = v

improvementindex = dict()

for key in basedic.keys():
    if key in moduledic:
        baseval = basedic[key]
        moduleval = moduledic[key]

        for i in range(len(baseval)):
            if baseval[i] != moduleval[i]:
                if key in improvementindex:
                    improvementindex[key] += [i]
                else:
                    improvementindex[key] = [i]
    else:
        print("le fichier " + key + " n'a pas de version compacte")

for key,val in improvementindex.items():
    print("Amélioration de la matrice " + key + " :")
    for i in val:
        print("    "+head[i] + " : base : " + basedic[key][i] + ", compacte : " + moduledic[key][i])
    print()

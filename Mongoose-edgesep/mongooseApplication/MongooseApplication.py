#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import csv
import sys

# All the file that it will be partitionning by Mongoose
filenames = [os.path.join(root,filename) for root, directories, filenames in os.walk('../Matrix/') for filename in filenames if filename]


with open(sys.argv[1], 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    # Each column of the csv file
    header = ["Total Edge Separator Time", "Matching", "Coarsening", "Refinement", "FM", "QP", "IO", "Cut Size", "Cut Cost", "Imbalance"]
    spamwriter.writerow(["matrixName"]+header)

    for testfile in filenames:
        try:

            cmd = './../Mongoose/Mongoose-edgesep/build/bin/mongoose '+testfile
            # Get back the result of Mongoose execution
            mongres = os.popen(cmd).read()

            # Delete first part of the output
            laststar = mongres.rindex("*")
            mongres = mongres[laststar+2:]

            # Split each line
            tab = mongres.split("\n")
            tab = tab[:len(tab)-1]
            resdic = dict()

            for line in tab:
                # Remove additional quote and split key:value
                line = line.replace("'","")
                key, value = line.split(":")
                if value != "":
                    resdic[key.strip()] = value.strip()

            reslist = [testfile.replace("../Matrix/","")]
            for col in header:
                reslist += [resdic[col]]
            spamwriter.writerow(reslist)
        except :
            # If Mongoose didn't return a good result (cause by a non sparse matrix or other). We keep a trace of matrix name but we didn't stop the script
            print(testfile + " fichier error")
            pass






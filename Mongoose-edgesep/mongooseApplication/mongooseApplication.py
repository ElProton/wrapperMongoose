#!/usr/bin/py
# -*- coding: utf-8 -*-
import subprocess
import os
import csv
import sys
import io

# All the file that it will be partitionning by Mongoose
filenames = [os.path.join(root,filename) for root, directories, filenames in os.walk('../../../Matrix/') for filename in filenames if filename]


with open(sys.argv[1], 'w', newline='') as csvfile:
	spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|',quoting=csv.QUOTE_MINIMAL)
	# Each column of the csv file
	
	header = ["Edge Cut",  "cut size", "cut cost", "imbalance"]
	
	spamwriter.writerow(["matrixName"]+header)
	

	for testfile in filenames:
		try:
			resdic = dict()

			cmd = 'cat '+testfile+' | ./connectorTest ' if sys.argv[2] != "B" else 'cat '+testfile+' | ./connectorTest --B'
			# Get back the result of Mongoose execution
			mongres = subprocess.Popen([cmd],stdout=subprocess.PIPE, shell=True)
			for line in mongres.stdout.readlines():
				line = line.decode("utf-8")
				# Remove additional quote and split key:value
				#line = line.replace("'","")
				print(line)
				key, value = line.split(":")
				if value != "":
					resdic[key.strip()] = value.strip()

			reslist = [testfile.replace("../../../Matrix/","")]
			for col in header:
				reslist += [resdic[col]]

			spamwriter.writerow(reslist)
				
		except :
			print(sys.exc_info()[0])
			pass


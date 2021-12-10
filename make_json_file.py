#!/usr/bin/env python3
import sys,re,io
from collections import defaultdict
import numpy as np
import pandas as pd
import resource

def main():
	meta, inFiles = getArgs()
	meta = storeMeta(meta = meta)
	getNames(meta = meta, inFiles = inFiles)

def getArgs():
	import argparse
	parser = argparse.ArgumentParser(description = "make json input for genome paint from ccle data and input file list")
	parser.add_argument('-meta', action = 'store', type = str, help = 'metadata table with ccle info')	
	parser.add_argument('-inFiles', action = 'store', type = str, help = 'list of bigwig files')	
	args = parser.parse_args()
	meta, inFiles =  args.meta, args.inFiles
	return(meta, inFiles)

def getNames(meta, inFiles):
	namesDict = defaultdict(dict)
	levelTracker = defaultdict(dict)
	levelCounter = 0
	### define max number of samples per level			
	maxSamps = 10
	for line in open(inFiles, 'r'):
		line = line.strip()
		srr = line.split('_')[2].split('.')[0]
		ind = meta[meta['rail_id'] == srr].index[0]
		
	
		cellLine = meta['cell_line'][ind]
		primary_site = meta['primary_site'][ind]
		if primary_site in levelTracker:
			levelTracker[primary_site]["samples"].append(line)
			if(len(levelTracker[primary_site]["samples"]) <= maxSamps):
				if not (meta.iloc[[ind]].isnull().values.any()):
					namesDict[line]["srr"] = srr
					namesDict[line]["cell_line"] = cellLine
					namesDict[line]["primary_site"] = primary_site
				
		else:
			levelTracker[primary_site]["level"] = levelCounter
			levelTracker[primary_site]["samples"] = [line]
			levelCounter += 1
			if not (meta.iloc[[ind]].isnull().values.any()):
				namesDict[line]["srr"] = srr
				namesDict[line]["cell_line"] = cellLine
				namesDict[line]["primary_site"] = primary_site
		
	### print the json file

	print("[\n" + 
		"\t{\n" + 
		"\t\t\"isfacet\": true,\n" + 
		"\t\t\"name\": \"CCLE-RNAseq\",\n" +
		"\t\t\"tracks\": ["
		)
	
	namesLen = len(namesDict)
	namesDictInd = 0
	for file in namesDict:
		namesDictInd += 1
		pSite = namesDict[file]["primary_site"]
		#if pSite in levelCountOutput:
		#	levelCountOutput[pSite] += 1
		#else:
		#	levelCountOutput[pSite] = 1

		level = levelTracker[pSite]["level"]
		sampName = file.split("_sums.")[1].split(".bw")[0]
		line = namesDict[file]["cell_line"]
		indent = "\t\t\t\t"
		sys.stdout.write("\t\t\t{\n" + 
			indent + "\"type\": \"bigwig\",\n" + 
			indent + "\"assay\": \"RNAseq\",\n" +
			indent + "\"sample\": \"" + str(sampName) + "\",\n" + 
			indent + "\"level1\": \"" + str(pSite) + "\",\n" +
			indent + "\"level2\": \"" + str(level) + "\",\n" +
			indent + "\"name\": \"" + str(line) + "\",\n" +
			indent + "\"file\": \"" + file + "\",\n" +
			"\t\t\t}"
		)
		if(namesDictInd < namesLen):
			print(',')
		else:
			sys.stdout.write('\n')

	print("\t\t]\n" + 
		"\t}\n" +
		"]"
		)

def storeMeta(meta):
	df = pd.DataFrame()
	for chunk in pd.read_table(meta, header = 0, chunksize = 1000, low_memory = False, index_col = None):
		df = pd.concat([df, chunk])
	return(df)
		

if __name__ == "__main__":
	main()

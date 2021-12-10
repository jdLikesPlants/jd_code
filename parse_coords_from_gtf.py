#!/usr/bin/env python3
import io,re,sys
from collections import defaultdict

for line in open(sys.argv[1], 'r'):
	line = line.strip()
	if line.startswith("#"):
		continue
	gene = line.split('"')[5]
	arr = line.split("\t")
	chromo = arr[0]
	left = arr[3]
	right = arr[4]
	if re.search("chr", chromo):
		print(gene + '\t' + chromo + ':' + left + '-' + right)

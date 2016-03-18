#!/usr/bin/env python

from __future__ import print_function

import re
import os
import os.path
import optparse
import sys
import csv
import random
import gc

def GenerateRandoms(listCG, dictMtbr, refFile, mtbrFile):
	listRandom = list(listCG)

	# fix the start and end 'N' s

	listHeadN = []
	listTailN = []
	startN = 0
	for b in listRandom:
		if b[0] == 'N' :
			listHeadN += [b]
			startN += 1
		else:
			break

	for i in range(len(listRandom) - 1, 0, -1):
		if listRandom[i][0] == 'N':
			listTailN += [listRandom.pop()]
		else:
			break

	listRandom = listRandom[startN + 1 : ]

	print("\tshuffling CG list ...")
	random.shuffle(listRandom)

	listRandom = listHeadN + listRandom + listTailN

	# calculate new pos

	print("\twritting new mtbr file ...")
	listRandom[0] += [1]
	for i in range(len(listRandom) - 1):
		listRandom[i + 1] += [listRandom[i][2] + len(listRandom[i][0])]

	# generate random mtbrFile

	mtbrFile.write('chr,pos,rC_n,rC_p,rT_n,rT_p\n')
	for b in listRandom :
		if(b[0] == 'CG'):
			if(b[1] in dictMtbr):
				mtbrItem = dictMtbr[b[1]]
				mtbrFile.write("%s,%ld,%s,%s,%s,%s\n" % ( 
					mtbrItem[0], b[2], 
					mtbrItem[2], mtbrItem[3],
					mtbrItem[4], mtbrItem[5]
				))
			else:
				print('error: CG pos ' + str(b[1]) + ' not found in CG list')

	# generate random refFile

	print('\twritting new reference seq file ...')
	i = 0
	refFile.write('> ' + '\n')
	for b in listRandom :
		refFile.write(b[0])
		if((i + 1) % 80  == 0):
			refFile.write('\n')
		i += len(b[0])

	del listRandom[:]
	gc.collect()

def main():

	# parse the command line options

	usage = 'usage: %prog [options] chr.fa chr.mtbr'
	parser = optparse.OptionParser(usage=usage, version='%prog 0.1.0')
	parser.add_option('-t', '--repeats-times', dest='repeats',
						help='files of repeats times to be generated')

	(options, args) = parser.parse_args()
	if(len(args) != 2):
		parser.print_help()
		sys.exit(0)
	
	refSeqFileName = args[0]
	mtbrFileName = args[1]
	baseFileName = os.path.splitext(os.path.basename(mtbrFileName))[0]
	
	repeats = 1
	if(options.repeats):
		repeats = int(options.repeats)

	if(not os.path.exists(mtbrFileName)):
		print('error: Failed to open mtbr file "', mtbrFileName, '"')
		sys.exit(-1)
	
	if(not os.path.exists(refSeqFileName)):
		print('error: Reference sequence file "', refSeqFileName, '"', ' doest not exist.')
		sys.exit(-1)

	# load reference sequence

	print('[*] Loading reference sequence ...')
	refSeq = ''
	with open(refSeqFileName, 'r') as refSeqFile :
		for line in refSeqFile:
			if(line[0] == '>'):
				continue

			refSeq += line.strip().upper()
	refSeqFile.close()
	
	# generate random elements
	# list format [['A', 123], ['C', 123], ['CG', 123], ['T', 123]]
	
	print('[*] Generating CG list ...')
	reCG = re.compile('[ATNatn]|CG|cg|[CGcg]')
	listCG = reCG.findall(refSeq)
	listCG = [[e, 1] for e in listCG]
	for i in range(len(listCG) - 1):
		listCG[i + 1][1] = listCG[i][1] + len(listCG[i][0])

	# load mtbr

	print('[*] Loading mtbr...')

	# mtbr data structure 
	
	dictMtbr = {}
	with open(mtbrFileName, 'rb') as csvfile :
		mtbrs = csv.reader(csvfile, delimiter = ',')
		next(mtbrs, None)  # skip the headers
		for mtbr in mtbrs:
			dictMtbr[int(mtbr[1])] = mtbr
	csvfile.close()

	# generate random samples

	for i in range(repeats):
		print('generating sample ' + str(i + 1) + '...')
		randomRefFileName = 'random.' + baseFileName + '.' + str(i + 1) + '.fa'
		randomMtbrFileName = 'random.' + baseFileName + '.' + str(i + 1) + '.mtbr'

		randomRefFile = open(randomRefFileName, 'w')
		randomMtbrFile = open(randomMtbrFileName, 'w')
		
		GenerateRandoms(listCG, dictMtbr, randomRefFile, randomMtbrFile)

		randomRefFile.close()
		randomMtbrFile.close()

	print('[*] Complete')

if __name__ == '__main__':
	main()
import sys
import datetime

#usage python parsePileup.py pileupfile outputfile depththreshold qualitythreshold
#example python parsePileup.py mypileup.txt output.txt 7 30
#output explanatin:
#column delimited file
#SCAFFOLD_ID	POSITION	BASE*	DEPTH	BASECOUNT(A C G T)***
#*Top 2 most common bases will be encapsulated in []
#**Repetition of A, C, G, T base count for N individuals from pileup file

class ParsePileup:
	
	infile = ''
	outfile = ''
	depthThreshold = ''
	qualityThreshold = ''

	#initialize
	def __init__(self):
		self.infile = sys.argv[1]
		self.outfile = sys.argv[2]
		self.depthThreshold = int(sys.argv[3])
		self.qualityThreshold = int(sys.argv[4])

	#return dictionary of key A,C,G,T,N	sorted ascendingly by base count
	def getMostCommonBase(self, mylist):	
		d = {'A':0,'C':0,'G':0,'T':0,'N':0}
		for item in mylist:
			try:
				d[item] += 1
			except KeyError:
				pass
		tmp = []
		for key, value in sorted(d.iteritems(), key=lambda (k,v): (v,k)):
			item = key, value
			if value > 0:
				tmp.append(key)
		return tmp
		
	#count each base			
	def countBase(self, mylist):	
		d = {'A':0,'C':0,'G':0,'T':0,'N':0}
		for item in mylist:
			item = item.upper()
			try:
				d[item] += 1
			except KeyError:
				pass
		return d['A'],d['C'],d['G'],d['T'],d['N']
			
	#parse each line in the pileup file		
	def parseString(self, string, quality, ref):
		string = list(string)
		count = 0
		parse = []
		refset = [',','.']
		nucset = ['A','C','G','T','N']
		countString = ''
		while count < len(string):
			item = string[count].upper()
			#case 1: match with reference base
			if item in refset:
				parse.append(ref)
			#case 2: match with other base
			elif item in nucset:
				parse.append(string[count])
			#case 3: sequence start
			elif item == '^':
				count += 1
				quality = string[count]
				count += 1
				startbase = string[count]
				if startbase == '.' or startbase == ',':
					parse.append(ref)
				else:
					parse.append(startbase)
			#case 4: sequence end
			elif item == '$':
				count -= 1
				endbase = string[count]
				if endbase == '.' or endbase == ',':
					parse.append(ref)
				else:
					parse.append(endbase)
				count += 1
			#case 5: indel
			elif item == '+' or item == '-': 				#indel, iterate until not digit
				count += 1
				start = count
				while string[count].isdigit():
					count += 1
				end = count
				span = int(''.join(string[start:end]))	#indel length
				indel = ''.join(string[end:end+span])	#bases for indel
				parse.append(indel)
				count += span
			#others
			else:
				parse.append(string[count])
			count += 1
		#check base quality
		quality = list(quality)
		lowqualCount = 0
		if len(quality) == len(parse):
			for q in range(len(quality)):
				qscore = ord(quality[q]) - 33
				if qscore < self.qualityThreshold: #30
					lowqualCount += 1
					parse[q] = ''					#if quality is lower than threshold, mask the base
		A,C,G,T,N = self.countBase(parse)
		
		#sorted base dictionary
		mcb = self.getMostCommonBase(parse)
	
		dict = {'A':A, 'C':C, 'G':G, 'T':T, 'N':N}		
		countString = '%s\t%s\t%s\t%s\t%s\t' %(A,C,G,T,N)
		return len(parse) - lowqualCount, dict, mcb								#return coverage, base count, sorted base dict
			
	def readfile(self):
		a = datetime.datetime.now()	
		f = open(self.infile, 'r')
		w = open(self.outfile, 'w')
		linenum = 1
		for line in f:
			line = line.strip().split('\t')
			count = 4
			depth = []
			string = ''
			cb = []
			gl = []
			while count < len(line):
				coverage, tstring, mcb = self.parseString(line[count], line[count+1],line[2])
				depth.append(coverage)
				gl.append(tstring)
				for item in mcb:
					cb.append(item)
				count += 3
			depth.sort()
			medianDepth = depth[len(depth)/2]
			print medianDepth, self.depthThreshold
			mcb = self.getMostCommonBase(cb)
			if medianDepth > self.depthThreshold: #7
				for item in gl:
					string += '%s\t%s\t%s\t%s\t' %(item['A'], item['C'], item['G'], item['T'])
				if len(mcb) == 1:
					mcbstring = mcb[0]
				else:
					mcbstring = ', '.join(mcb[-2:])
					mcbstring = '[%s]' %(mcbstring)
					if len(mcb) > 2:
						mcbstring = '%s, %s' %(mcbstring, ', '.join(mcb[:-2]))
				wstring = '%s\t%s\t%s\t%s\t%s\n' % (line[0], line[1], mcbstring, medianDepth, string)
				w.write(wstring)
	
			if linenum%10000 == 0:
				print 'line=%d,'%linenum
				b = datetime.datetime.now()
				c = b - a
				print 'task takes',c.seconds,'s,',c.microseconds,'ms'
			linenum += 1
		f.closed
		w.flush()
		w.closed

	#run	
	def run(self):	
		self.readfile()
		
		
pp = ParsePileup()
pp.run()

# circular genome mapper tool
#!/usr/bin/env python 

from subprocess import call
import re
import itertools
from argparse import ArgumentParser

singleReads = False
fwdReads = False
revReads = False

parser = ArgumentParser()
parser.add_argument("-s", "--singles", dest='singleReads',
                    help="single ended FASTQ file", metavar="FILE")
parser.add_argument("-rl", "--read_length", dest='readLength', default = 250,
                    help="read length", metavar="INT")
parser.add_argument("-r", "--ref_genome", dest='refGenome',
                    help="reference genome in FASTA format", metavar="FILE")
parser.add_argument("-1", "--fwd", dest='fwdReads',
                    help="paired end forward reads FASTQ file", metavar="FILE")
parser.add_argument("-2", "--rev", dest='revReads',
                    help="paired end reverse reads FASTQ file", metavar="FILE")
parser.add_argument("-o", "--output", dest='prefix',
                    help="prefix for output files", metavar="STR")



args = parser.parse_args()
singleReads = args.singleReads
readLength = int(args.readLength) + 1
userRef = open(args.refGenome, "r")
fwdReads = args.fwdReads
revReads = args.revReads
output = args.prefix + "."


# get the length of reference genome from FASTA
x = open(args.refGenome, "r")
bases = ""
for line in x:
	if line[0] == ">":
		pass
	else:
		bases = bases + line
		bases = bases.replace('\n', '')				
refLength = len(bases)
print (refLength)

def extendRef(userRef, readLength):
	# Takes first rl bases of reference sequence and appends to the end
	# create a file to write the extended reference genome to
	extendedRef = open('reference_genome.ext.fasta', 'w')
	# write entire userRef to extendedRef, then return to start of userRef
	extendedRef.write(userRef.read())
	userRef.seek(0)
	# variable to contain bases to be appended
	firstBases = ""
	# until firstBases contains more characters as readLength, each line is appended to it
	while len(firstBases) < readLength:
		f = userRef.readline()
		if f[0] == ">":
			pass
		else:
			firstBases = firstBases + f
			firstBases = firstBases.replace('\n', '')				
	# if firstBases is longer than readLength, it truncated to be the same length
	if len(firstBases) > readLength:
		firstBases = firstBases[:readLength]
	# write firstBases to the end of the reference genome
	extendedRef.write(firstBases)
	extendedRef.close()



def indexRef():
	# Indexes new reference genome using Bowtie2
	call (["bowtie2-build", "reference_genome.ext.fasta" ,output+"ext_ref", "--quiet"])

	
def singlesMapper(filename):
	# maps reads against extended reference genome using Bowtie2
	# -k 2 means that reads that map twice will be added to the SAM file twice
	call (["bowtie2", "-x", output + "ext_ref", "-U" , filename, "-S", output+ "mapped.SAM" , "--quiet"])

def pairMapper(fwd, rev):
	# maps reads against extended reference genome using Bowtie2
	# -k 2 means that reads that map twice will be added to the SAM file twice
	call (["bowtie2", "-x", output + "ext_ref", "-1" , fwd,  "-2" , rev,  "-S", output+ "mapped.SAM"])

def renumber(SAMfile):
	# When we extend the genome, some reads will map twice. Bowtie only includes one read by default,
	# but it could be either. this function will correct the position of reads that map only to the extended part of the genome
	f = open(SAMfile, 'r')
	g = open(output + "renumbered.SAM", 'w')
	for line in f:
		if line[0] == "@":
			g.write(line)
		else:
			x = line.split()
			# also need PNEXT (location of mate) to be changed
			if int(x[7]) > refLength:
				x[7] = str(int(x[7])-refLength)
			if int(x[3]) > refLength:
				x[3] = int(x[3]) - refLength
				x[3] = str(x[3])
				# change MAPQ because multimapping is our fault
				#if x[4] == "1":
				x[4] = "42"
			elif int(x[3]) < readLength:
				x[4] = "42"
			x = "\t".join(x)
			x = x + "\n"
			g.write(x)
	g.close()

		
def cigarLong(cigar):
	# this function takes a cigar in the format 35M3D1I95M and converts it to something like MMMMMMMMMMMDDIIIIMMMMMMMMM
	# this means it can be split like the sequence/quality, before being converted back
	match = (re.split('(\d+)', cigar))
	match.pop(0)
	count = 1
	longCigar = ""
	while count <= len(match):
		longCigar = longCigar + match[count] * int(match[count-1])
		count = count + 2
	return longCigar

def cigarShort(cigar):
	# this function does the opposite of cigarLong - converts cigar in MMMMMMDDIMMMMMMMMMMM format back to proper format
	groups = []
	for _, g in itertools.groupby(cigar):
		groups.append(''.join(g))
	shortCigar = ""
	for item in groups:
		shortCigar = shortCigar + str(len(item)) + item[0]
	return shortCigar
	
	
def readSplitter(SAMfile):
	# parses SAM file searching for reads that span the break
	# need two variables - start position [3] and sequence [9] length 
	# if sequence length + start position are longer than refLength, it needs splitting
	f = open(SAMfile, 'r')
	g = open(output + "split.SAM", 'w')
	for line in f:
		if line[0] == "@":
			pass
			g.write(line)
		else:
			x = line.split()
			# check to see if read spans break
			if int(x[3]) + len(x[9]) - 1 > refLength:
				#Call function
				lineA, lineB = fwdSplit(x)
				# if TLEN is pos and strand is rev then unpair read
				if int(x[8]) > 0 and int(x[1])&16 == 16:
					lineB = "\t".join(unpair(lineB.split()))+"\n"
				# if TLEN is neg and strand is fwd then unpair read
				if int(x[8]) < 0 and int(x[1])&16 != 16:
					lineA = "\t".join(unpair(lineA.split()))+"\n"
				# write both new lines if mapq > 0
				if int(lineA.split()[4]) > 0:
					g.write(lineA)
				if int(lineB.split()[4]) > 0:
					g.write(lineB)
			else:
				#if line does not need splitting, write as is
				if int(line.split()[4]) > 0:
					g.write(line)			
	g.close()
	
	
def fwdSplit(x):
	#This function takes a line (from a fwd strand read), splits it, and returns both lines
	#split the nucleotide sequence into two strings
	seqA = x[9][:len(x[9]) - (int(x[3]) + int(len(x[9])) - refLength)+1]
	if (len(seqA) + int(x[3])) > 5387:
		print("long")
	seqB = x[9][len(x[9]) - (int(x[3]) + int(len(x[9])) - refLength)+1:]
	# split the quality score into two strings
	qualA = x[10][:len(x[9]) - (int(x[3]) + int(len(x[9])) - refLength)+1]
	qualB = x[10][len(x[9]) - (int(x[3]) + int(len(x[9])) - refLength)+1:]
	# take the CIGAR from the list, convert it into long format, split into two, then convert both parts back to std format
	cigar = cigarLong(x[5])
	# need to stop deletions counting towards length
	cigA = cigar[:len(x[9]) - (int(x[3]) + int(len(x[9])) - refLength)+1]
	deletions = cigA.count("D")
	cigA = cigarShort(cigar[:len(x[9]) + deletions - (int(x[3]) + int(len(x[9])) - refLength)+1])
	cigB = cigar[len(x[9]) - (int(x[3]) + int(len(x[9])) - refLength)+1:]
	cigB = cigarShort(cigar[len(x[9]) + deletions - (int(x[3]) + int(len(x[9])) - refLength)+1:])

	# make two new lines, based on the original line but substituting in the new sequence, quality and CIGAR	
	lineA = x
	lineA[9] = seqA
	lineA[10] = qualA
	lineA[5] = cigA
	lineA = "\t".join(lineA)+"\n"
	lineB = x
	lineB[9] = seqB
	lineB[10] = qualB
	lineB[5] = cigB
	# starts at position 1
	lineB[3] = "1"
	lineB = "\t".join(lineB)+"\n"
	# return both lines
	return lineA, lineB

def revSplit(line):
	#This function takes a line (from a reverse strand read), splits it, and returns both lines
	#split the nucleotide sequence into two strings
	# FUNCTION PROBABLY NOT NEEDED
	seqA = x[9][:len(x[9]) - (int(x[3]) + int(len(x[9])) - refLength)+1]
	seqB = x[9][len(x[9]) - (int(x[3]) + int(len(x[9])) - refLength)+1:]
	# split the quality score into two strings
	qualA = x[10][:len(x[9]) - (int(x[3]) + int(len(x[9])) - refLength)+1]
	qualB = x[10][len(x[9]) - (int(x[3]) + int(len(x[9])) - refLength)+1:]
	# take the CIGAR from the list, convert it into long format, split into two, then convert both parts back to std format
	cigar = cigarLong(x[5])
	cigA = cigarShort(cigar[:len(x[9]) - (int(x[3]) + int(len(x[9])) - refLength)+1])
	cigB = cigarShort(cigar[len(x[9]) - (int(x[3]) + int(len(x[9])) - refLength)+1:])
	# make two new lines, based on the original line but substituting in the new sequence, quality and CIGAR	
	lineA = x
	lineA[9] = seqA
	lineA[10] = qualA
	lineA[5] = cigA
	lineA = "\t".join(lineA)+"\n"
	lineB = x
	lineB[9] = seqB
	lineB[10] = qualB
	lineB[5] = cigB
	# starts at position [length of reference genome]
	lineB[3] = str(refLength)
	lineB = "\t".join(lineB)+"\n"
	# return both lines
	return lineA, lineB

def unpair(line):
	#this function turns reads that are part of a pair into single reads
	# Set TLEN to 0
	line[8] = "0"
	# Set PNEXT to 0
	line[7] = "0"
	# Set RNEXT to *
	line[6] = "*"
	# Change flag. Remove 1, 2, 8, 32, 64, 128 if turned on.
	line[1] = str(int(line[1])&3860)
	return line
	
	
# PAIRED END READS - STUFF TO ADD
# If read is FWD, the 3' end remains paired with the reverse
# If read is REV, the 5' end remains paired with the forward
# Remainder is added to SAM file as an unpaired read	

extendRef(userRef, readLength)
indexRef()
# change the next couple of lines to determine if input is a valid file
if type(singleReads) == str:
	singlesMapper(singleReads)
elif type(fwdReads) == str and type(revReads) == str:
	pairMapper(fwdReads, revReads)
else:
	print ("No valid input files")	
renumber(output + "mapped.SAM")
readSplitter(output + "renumbered.SAM")
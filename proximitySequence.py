'''
Take three arguments, the 
fasta file to use, the string to search for, 
and the number of nucleotides on either side to include.


'''
from Bio.Seq import Seq
from collections import defaultdict
import sys
import re

fileName = sys.argv[1]
oligo = sys.argv[2]
totalLength = sys.argv[3]

f = open(fileName, 'r')
o = open("TargetProtector of length_" + totalLength + "_" + fileName, 'w')

list = defaultdict(str)


name =''
oligoLength = len(oligo)
proxLength = (int(totalLength) - oligoLength)/2

#print("Oligo Length is %s and the proxLength is %s" % (oligoLength, proxLength))

for line in f:
		#if your line starts with a > then it is the name of the following sequence
		if line.startswith('>'):
			name = line[1:-1]
			continue #this means skips to the next line
		#this code is only executed if it is a sequence of bases and not a name
		list[name]+=line.strip()

match = re.search(oligo, list[name])
position = match.start()

#print(list[name])
#print(list[name][position:position+oligoLength])
desiredOligo = Seq(list[name][position-proxLength:position+oligoLength+proxLength])
desiredOligoRC = desiredOligo.reverse_complement()
#print(desiredOligo)
#print(desiredOligoRC)
print("Your desired target protector for %s, covers region, %s, with a sequence of %s" % (oligo, desiredOligo, desiredOligoRC))
o.write("Your desired target protector for %s, covers region, %s, with a sequence of %s" % (oligo, desiredOligo, desiredOligoRC))

f.close()
o.close()

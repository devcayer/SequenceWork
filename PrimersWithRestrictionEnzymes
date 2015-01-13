'''
Take one argument, the fileName as a "string.txt" 

store fasta sequence in a dictionary

take restriction site sequences as input, returns primers that will amplify ends of sequence with 
restriction sites appended and a 6 nucleotide spacer TAGATC

'''
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from collections import defaultdict
import sys

fileName = sys.argv[1]
meltCutoff = float(sys.argv[2])
restriction1 = sys.argv[3]
restriction2 = sys.argv[4]

restriction1 = Seq(restriction1)
restriction2 = Seq(restriction2)

fixedSequence = Seq('TAGATC')

end1 = fixedSequence + restriction1
end2 = fixedSequence + restriction2

#print("end1 is %s, and end2 is %s" % (end1, end2))

f = open(fileName, 'r')
o = open("terminalPrimers_" + fileName, 'w')
oIDT = open("IDTterminalPrimers_" + fileName, 'w')

list = defaultdict(str)

name =''

for line in f:
		#if your line starts with a > then it is the name of the following sequence
		if line.startswith('>'):
			name = line[1:-1]
			continue #this means skips to the next line
		#this code is only executed if it is a sequence of bases and not a name
		list[name]+=line.strip()

#print(list)

o.write("Sequence\t5'Primer\tTm\t3'Primer\tTm\n")

for key, value in list.iteritems():

	length5 = 40
	length3 = 40
	primer5 = Seq(value[:length5])
	primer3 = Seq(value[-length3:])
#	print("meltcutoff is %s, primer5 is %s, and primer 3 is %s" % (meltCutoff, primer5, primer3))
#	print("melting temp of primer5 is %s and mt of primer3 is %s" % (mt.Tm_NN(primer5), mt.Tm_NN(primer3)))
#	print(mt.Tm_NN(primer5) > meltCutoff)
	while ((mt.Tm_NN(primer5, Na=50, Tris=10, Mg=2.0, dNTPs=0.5, saltcorr=7) >= meltCutoff)):
		length5 = length5 - 1
		primer5 = Seq(value[:length5])
		meltingTemp5 = mt.Tm_NN(primer5, Na=50, Tris=10, Mg=2.0, dNTPs=0.5, saltcorr=7)
	#	print("length is %s and primer 5 is  %s and primer.mt is %s" % (length5, primer5, meltingTemp5))
	

	while ((mt.Tm_NN(primer3, Na=50, Tris=10, Mg=2.0, dNTPs=0.5, saltcorr=7) >= meltCutoff)):
		length3 = length3 - 1
		primer3 = Seq(value[-length3:])
		primer3 = primer3.reverse_complement()
		meltingTemp3 = mt.Tm_NN(primer3, Na=50, Tris=10, Mg=2.0, dNTPs=0.5, saltcorr=7)
	#	print("length is %s and primer 3 is  %s and primer3.mt is %s" % (length3, primer3, meltingTemp3))

	o.write("%s\t%s\t%0.2f\t%s\t%0.2f\n" % (key, end1+primer5, meltingTemp5, end2+primer3, meltingTemp3))
	#o.write("%s \n A Primer(insert S): %s \n B Primer(insert AS): %s \n C Primer(vector AS): %s \n D Primer(vector S): %s \n" % (key, aPrimer, bPrimer, cPrimer, dPrimer))
	oIDT.write("%s-S;%s;25nm;STD\n" % (key, end1+primer5))
	oIDT.write("%s-AS;%s;25nm;STD\n" % (key, end2+primer3))

f.close()
o.close()
oIDT.close()

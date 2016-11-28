#Anna Maule
#   This is a script to isolate a specific segment of the DNA
# and it writes the sequence on another file, as well as the first line
# of the genome fasta file.

import Bio
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC

a =  int(input("value for begining of the sequence: "))
b =  int(input("value for end of the sequence position:"))
filename = input("File name:")
fasta = input("Fasta file name: ")

for seq_record in SeqIO.parse(fasta,"fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

    f= open(filename,"w")
    f2 = open(fasta,"r")
    line = f2.readline()

    f.write(line)
    f.write(str(seq_record.seq[a:b]))
    
    f.close()
    f2.close()

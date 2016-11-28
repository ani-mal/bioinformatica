import Bio
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
import plotly.plotly as py
import plotly.graph_objs as go
import plotly.tools as tls

username = input("Plotly username: ")
apiKey = input("Plotly API Key: ")

tls.set_credentials_file(username, apiKey)
homo=[]
pantro=[]

for seq_record in SeqIO.parse("homoSapiensMito.fasta","fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

    homo.append(seq_record[1:16568].seq.count('A'))
    homo.append(seq_record[1:16568].seq.count('T'))
    homo.append(seq_record[1:16568].seq.count('C'))
    homo.append(seq_record[1:16568].seq.count('G'))

for seq_record in SeqIO.parse("pantroglodytes.fasta","fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    size = len(seq_record)

    pantro.append(seq_record[1:16568].seq.count('A'))
    pantro.append(seq_record[1:16568].seq.count('T'))
    pantro.append(seq_record[1:16568].seq.count('C'))
    pantro.append(seq_record[1:16568].seq.count('G'))

trace1=go.Scatter(x=['A','T','C','G'], y=[(homo[0]/size),(homo[1]/size),(homo[2]/size),(homo[3]/size)])
trace2=go.Scatter(x=['A','T','C','G'], y=[(pantro[0]/size),(pantro[1]/size),(pantro[2]/size),(pantro[3]/size)])

data=[trace1,trace2]
py.plot(data, filename="genes")

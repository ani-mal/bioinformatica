import Bio
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
import plotly.plotly as py
import plotly.graph_objs as go
import plotly.tools as tls

#username = input("Plotly username: ")
#apiKey = input("Plotly API Key: ")

username= "pawpepe"
apiKey = "07lxwcy4d2"
tls.set_credentials_file(username, apiKey)

bacteria1=[]
bacteria2=[]
bacteria3=[]

for seq_record in SeqIO.parse("sequence.fasta","fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

    maxLen= 48502



    windowSize1 = int(input("Window Size1: "))
    windowSize2 = int(input("Window Size2: "))
    windowSize3 = int(input("Window Size3: "))

    x1 = int(windowSize1)
    x2 = int(windowSize2)
    x3 = int(windowSize3)
    x=1

    for y in range(0,8):
      bacteria1.append(seq_record[x:x1].seq.count('C') + seq_record[x:x1].seq.count('G'))
      x = x1
      x1 += windowSize1

    x=1
    for y in range(0,22):
      bacteria2.append(seq_record[x:x2].seq.count('C') + seq_record[x:x2].seq.count('G'))
      x = x2
      x2 += windowSize2
    x=1
    for y in range(0,44):
      bacteria3.append(seq_record[x:x3].seq.count('C') + seq_record[x:x3].seq.count('G'))
      x = x3
      x3 += windowSize3


trace1=go.Scatter(x=[100,200,300,400,500,600,700,800,900], y=[bacteria1[0],bacteria1[1],bacteria1[2],bacteria1[3],bacteria1[4],bacteria1[5],bacteria1[6],bacteria1[7]])
trace2=go.Scatter(x=[2000,4000,6000,8000,10000,12000,14000,16000], y=[bacteria2[0],bacteria2[1],bacteria2[2],bacteria2[3],bacteria2[4],bacteria2[5],bacteria2[6],bacteria2[7]])
trace3=go.Scatter(x=[5500,11000,16500,22000,27500,33000,38500,44000], y=[bacteria3[0],bacteria3[1],bacteria3[2],bacteria3[3],bacteria3[4],bacteria3[5],bacteria3[6],bacteria3[7]])
data=[trace1,trace2,trace3]
py.plot(data, filename="bacteriacount")

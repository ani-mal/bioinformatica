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

tls.set_credentials_file(username, apiKey)
homo1 = []
homo2 = []
homo3 = []
chimp1 = []
chimp2 = []
chimp3 = []


windowSize1 = int(input("Window Size1: "))
windowSize2 = int(input("Window Size2: "))
windowSize3 = int(input("Window Size3: "))


for seq_record in SeqIO.parse("homoSapiensMito.fasta","fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))



    x1 = int(windowSize1)
    x2 = int(windowSize2)
    x3 = int(windowSize3)
    x=1


    for y in range(0,8):
      homo1.append(seq_record[x:x1].seq.count('C') + seq_record[x:x1].seq.count('G'))
      x = x1
      x1 += windowSize1

    x=1
    for y in range(0,8):
      homo2.append(seq_record[x:x2].seq.count('C') + seq_record[x:x2].seq.count('G'))
      x = x2
      x2 += windowSize2
    x=1
    for y in range(0,8):
      homo3.append(seq_record[x:x3].seq.count('C') + seq_record[x:x3].seq.count('G'))
      x = x3
      x3 += windowSize3


for seq_record in SeqIO.parse("pantroglodytes.fasta","fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

    x1 = int(windowSize1)
    x2 = int(windowSize2)
    x3 = int(windowSize3)
    x=1

    for y in range(0,8):
      chimp1.append(seq_record[x:x1].seq.count('C') + seq_record[x:x1].seq.count('G'))
      x = x1
      x1 += windowSize1

    x=1
    for y in range(0,8):
      chimp2.append(seq_record[x:x2].seq.count('C') + seq_record[x:x2].seq.count('G'))
      x = x2
      x2 += windowSize2
    x=1
    for y in range(0,8):
      chimp3.append(seq_record[x:x3].seq.count('C') + seq_record[x:x3].seq.count('G'))
      x = x3
      x3 += windowSize3

trace1=go.Scatter(x=[1,2,3,4,5,6,7,8], y=[homo1[0],homo1[1],homo1[2],homo1[3],homo1[4],homo1[5],homo1[6],homo1[7]])
trace2=go.Scatter(x=[1,2,3,4,5,6,7,8], y=[chimp1[0],chimp1[1],chimp1[2],chimp1[3],chimp1[4],chimp1[5],chimp1[6],chimp1[7]])
data=[trace1,trace2]
py.plot(data, filename="genes 1")


trace1=go.Scatter(x=[1,2,3,4,5,6,7,8], y=[homo2[0],homo2[1],homo2[2],homo2[3],homo2[4],homo2[5],homo2[6],homo2[7]])
trace2=go.Scatter(x=[1,2,3,4,5,6,7,8], y=[chimp2[0],chimp2[1],chimp2[2],chimp2[3],chimp2[4],chimp2[5],chimp2[6],chimp2[7]])

data=[trace1,trace2]
py.plot(data, filename="genes 2")



trace1=go.Scatter(x=[1,2,3,4,5,6,7,8], y=[homo3[0],homo3[1],homo3[2],homo3[3],homo3[4],homo3[5],homo3[6],homo3[7]])
trace2=go.Scatter(x=[1,2,3,4,5,6,7,8], y=[chimp3[0],chimp3[1],chimp3[2],chimp3[3],chimp3[4],chimp3[5],chimp3[6],chimp3[7]])

data=[trace1,trace2]
py.plot(data, filename="genes 3")

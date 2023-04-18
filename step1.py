from Bio import SeqIO
import numpy as np
import math
import os
import sys
import pandas as pd
import re
import shutil
import argparse
from pprint import pprint
import openpyxl
import textwrap

def en_index():
    parser = argparse.ArgumentParser(
        description='DNA疫苗自动化设计软件V2.0  AutoDVD2 --version 2.0',
        formatter_class=argparse.RawTextHelpFormatter,
        usage='type "python3 %(prog)s --help" for more information'
    )
    parser.add_argument('-p', '--pdb_name', type=str, required=True,
                        help=textwrap.dedent('''\
    Must-option. 
    The user can enter the pathogen protein PDB number or the pathogen protein name. Note that the PDB number needs to be capitalized, 
    and the pathogen protein name must not contain spaces. 
    If you need to enter multiple names, please connect them with "+", for example: "5X5C+5K6G". '''))
    
    parser.add_argument('-c', '--chain', type=str, required=True,
                       help=textwrap.dedent('''\
    Must-option. 
    Note that the name of the selected chain needs to be capitalized. If you need to enter multiple names, you need to connect them with "+", for example: "A+A". 
    Note that the writing order of the chain and PDB number or pathogen protein name should correspond.'''))
    
    parser.add_argument('-if', "--in_fasta_file", type=str, required=True,
                        help=textwrap.dedent('''\
    Must-option. 
    The name format of the fasta file should be "<capitalized PDB number or pathogen protein name>.fasta". 
    Before uploading, please check whether the fasta file only contains sequences of one specific chain.
    If not, you need to manually delete redundant amino acid sequences. If you need to input multiple files, the absolute paths of the file should be connected with "+", 
    for example: "./AutoDVD2/demo/5X5C.fasta+./AutoDVD2/demo/5K6G.fasta”. Note that the writing order of the paths should correspond to the PDB number or the name of the pathogen protein.'''))
    
    parser.add_argument('-m', '--msa_file', type=str,
                        help=textwrap.dedent('''\
    Please enter the absolute path of multiple sequence alignment file. We recommend users to use MEGA to process multiple sequence alignment file.
    Note that the ID of the destination sequence should be named 'Query'.If you need to input multiple files, the absolute paths of the file should be connected with "+". '''))
    
    parser.add_argument('-w', '--way', type=str, required=True,
                        help=textwrap.dedent('''\
    If you choose to upload multiple sequence alignment file, please enter 'msa'. If you want to directly input the fasta file and use 'blastp' in AutoDVD2 
    for multiple sequence alignment, please enter 'blast'.'''))
    
    args = parser.parse_args()
    pprint(args.__dict__)
    return args

def trans0(input,output):    
    fi = open(input,"r")
    fo = open(output,"a")
    wflag =False 
    newline = []
    new=[] 
    for line in fi :
        if "Lambda" in line:
            break 
        if "RecName" in line: 
            wflag = True
            continue
        if wflag == True:
            K = list(line)
            if len(K)>1:
                for i in K : 
                    newline.append(i)
    strlist="".join(newline)
    newlines=str(strlist)  
    fo.write(newlines)
    fi.close()
    fo.close()
    
    fo=open(output,"r")
    for line in fo:
        index=re.search("\s\d",line).span()
        L=list(line)
        a=index[1]
        str1=L[a+2:]
        istr=str(str1)
        str2=istr.replace(' ' ,'-')
        nline=line.replace(istr,str2)
        X=list(nline)
        for i in X:
            new.append(i)
    list1="".join(new)
    news=str(list1)
    fo.close()
    
    fo=open(output,"w")
    fo.write(news)
    fo.close()
   
def trans1(input1,output1):    
    fi = open(input1,"r")
    fo = open(output1,"w")
    new=[]
    news=[]
    new1=[]
    news1=[]
    for line in fi:
        L=list(line)
        if len(L)>1:
            index1=re.search("\s\d",line).span()
            a=index1[1]
            str1=L[a+4:]
            str2=L[:a-3]
            istr1="".join(str1)
            istr2="".join(str2)
            str3=istr1.replace(' ' ,'-')
            str4=">"+istr2+" "+str3
            X=list(str4)
            for i in X:
                new.append(i)
            new.append("\n")
    news="".join(new)
    fo.write(news)
    fo.close()
    fi.close()
    
    ft=open(output1,"r")
    for line in ft:
        A=list(line)
        if len(A)>1:
            index2=re.search("--\d",line).span()
            b=index2[0]
            str5=A[:b]
            for i in str5:
                if i==" ":
                    new1.append("\n")
                else:
                    new1.append(i)
            new1.append("\n")
    news1="".join(new1)  
    ft.close()
    fs=open(output1,"w")
    fs.write(news1)
    fs.close()
    
    fa=open(output1,"r")
    D=[]
    for line in fa:
        if len(list(line))>1:
            for i in list(line): 
                D.append(i)
    str6="".join(D)
    fa.close()
    fb=open(output1,"w")
    fb.write(str(str6))
    fb.close()
    
    
def mkfile(output1,fpath):
    num=0
    linum=0
    li=[]
    pi=[]
    ai=[]
    fi=open(output1,"r")
    for line in fi:
        linum+=1
        pi.append(line)
        if "Query_1" in line:
            num+=1
            li.append(linum)
    fi.close()
    for i in range(num):
        if num-1==i:
            pr="".join(pi[(li[i]-1):])
        else: 
            pr="".join(pi[(li[i]-1):(li[i+1]-1)])
        spr=str(pr)
        ai.append(spr)
   
    os.makedirs(fpath)
    for i in range(0,num):
        ppath=fpath+"/"+str(i)+".msa.fasta"
        fx=open(ppath,"a")
        fx.write(str(ai[i]))
        fx.close()    

def loadAASeq(infile):
    seq = {}
    for i in SeqIO.parse(infile,'fasta'):
        seq[i.id] = list(i.seq)
    return seq

def AAHash():
    AARow = {"A":0,"R":1,"N":2,"D":3,"C":4,"Q":5,"E":6,"G":7,"H":8,"I":9,"L":10,"K":11,"M":12,"F":13,"P":14,"S":15,"T":16,"W":17,"Y":18,"V":19,"-":20,"X":21,'J':22,"B":2,"Z":5}
    return AARow


def Quantile(Vec):
    Vlen = len(Vec)
    SortIndex = Vec.argsort()
    
    VRank = np.zeros(Vlen)
    forwardIndex = list()
    forwardValue = 0
    for i in range(0,Vlen):
        if Vec[SortIndex[i]] == forwardValue:
            forwardIndex.append(SortIndex[i])
            if i == Vlen-1:
                for j in forwardIndex:
                    VRank[j] = i+1
        else:
            if len(forwardIndex) >= 2:
                for j in forwardIndex:
                    VRank[j] = i
            else:
                if i > 0:
                    VRank[SortIndex[i-1]] = i
            if i == Vlen-1:
                VRank[SortIndex[i]] = i+1
            forwardIndex = [SortIndex[i]]
            forwardValue = Vec[SortIndex[i]]
    PercValue = VRank*100/Vlen
    return np.around(PercValue,decimals=2)

def ppssm(inp,protein):
    os.system("psiblast -query "+inp+" -db ./AutoDVD2/blast/db/swissprot/swissprot -num_iterations 4 -out ./"+protein+"_"+chain+".txt -out_ascii_pssm ./"+protein+"_"+chain+".pssm")
    inf="./"+protein+"_"+chain+".pssm"
    fp=open(inf,"r")
    a=len(fp.readlines())-6
    fp.close()
    fp1=open(inf,"r")
    a1=[]
    for line in fp1.readlines()[3:a]:
        b=list(line)
        if b[6] =="A":
            a1.append(str(b[11])+str(b[12]))
        elif b[6] =="R":
            a1.append(str(b[15])+str(b[16]))
        elif b[6]=="N":
            a1.append(str(b[19])+str(b[20]))
        elif b[6]=="D":
            a1.append(str(b[23])+str(b[24]))
        elif b[6]=="C":
            a1.append(str(b[27])+str(b[28]))
        elif b[6]=="Q":
            a1.append(str(b[31])+str(b[32]))
        elif b[6]=="E":
            a1.append(str(b[35])+str(b[36]))
        elif b[6]=="G":
            a1.append(str(b[39])+str(b[40]))
        elif b[6]=="H":
            a1.append(str(b[43])+str(b[44]))
        elif b[6]=="I":
            a1.append(str(b[47])+str(b[48]))
        elif b[6]=="L":
            a1.append(str(b[51])+str(b[52]))
        elif b[6]=="K":
            a1.append(str(b[55])+str(b[56]))
        elif b[6]=="M":
            a1.append(str(b[59])+str(b[60]))
        elif b[6]=="F":
            a1.append(str(b[63])+str(b[64]))
        elif b[6]=="P":
            a1.append(str(b[67])+str(b[68]))
        elif b[6]=="S":
            a1.append(str(b[71])+str(b[72]))
        elif b[6]=="T":
            a1.append(str(b[75])+str(b[76]))
        elif b[6]=="W":
            a1.append(str(b[79])+str(b[80]))
        elif b[6]=="Y":
            a1.append(str(b[83])+str(b[84]))
        else:
            a1.append(str(b[87])+str(b[88]))   
    fp1.close()
    path1="./step1-calcon-"+protein+"_"+chain+"/"+protein+"_"+chain+".pssm"
    shutil.move(inf,path1)
    path2="./"+protein+"_"+chain+".txt"
    path3="./step1-calcon-"+protein+"_"+chain+"/"+protein+"_"+chain+".txt"
    shutil.move(path2,path3)
    return a1

def ddata(inp,protein,fpath,chain):
    path3="./step1-calcon-"+protein+"_"+chain+"/con-"+protein+"_"+chain+".xlsx"
    list0=[]
    list1=[]
    list3=[]
    num1=len(os.listdir(fpath))

# 储存目标氨基酸序列
    for i in range(num1):
        in_msa=open(fpath+"/"+str(i)+".msa.fasta","r")
        nd=in_msa.readlines()[1]
        nd1=nd.rstrip()
        for i in list(nd1):
            if i !="-":
                list0.append(i)
        in_msa.close()

    for i in range(num1):
        in_msa=open(fpath+"/"+str(i)+".msa.fasta","r")
        nseq = "Query_1"
        Seq = loadAASeq(in_msa)
        AA = AAHash()
        seqlenAll = len(Seq[nseq])
        seqlenInt = seqlenAll - Seq[nseq].count('-')
        Count = np.zeros((25,seqlenInt),dtype = float) + 0.000001
        gaps = 0   
        
        for i in range(0,seqlenAll):
            if Seq[nseq][i] == '-':
                gaps += 1
            else:
                for key in Seq.keys():
                    aa = Seq[key][i]
                    row = int(AA[aa])
                    col = i-gaps
                    Count[row,col] += 1
   
        colSum = sum(Count[...,0])
        CountRate = Count/colSum
        Shannon = -(CountRate*np.log2(CountRate)).sum(axis=0)
        Rindex = math.log2(25)-Shannon
        Rrelat = (Rindex - Rindex.min())/(Rindex.max()-Rindex.min())
        Rquant = Quantile(Rindex)
        for i in range(len(Rquant)):
            list1.append(str(Rquant[i]))
        in_msa.close()

    h=len(list0)+1
    for i in range(1,h):
        list3.append(str(i))
       
    list2=ppssm(inp,protein)
    
    list11=[]
    for i in range(h-1):
        list11.append(float(list1[i]))
    list22=[]
    for i in range(len(list2)):
        list22.append(int(list2[i]))
           
    ax=60
    con=[]
    for i in range(h-1):
        if list11[i]>=ax and list22[i]>0:
            con.append(str(i+1))
    path0="./step1-calcon-"+protein+"_"+chain+"/con-pos-"+protein+"_"+chain+".txt"
    fx=open(path0,"a")
    fx.write("\n".join(con))
    
    data=list(zip(list3,list0,list1,list2))
    df=pd.DataFrame(data,columns=["pos","amino-acid","entropy","pssm-score"])
    df.to_excel(path3,engine='openpyxl')
    
def ddata_msa(in_msa,fasta_file,total_path,protein,chain):
    nseq='Query'
    
    Seq = loadAASeq(in_msa)
    AA = AAHash()
    seqlenAll = len(Seq[nseq])
    seqlenInt = seqlenAll - Seq[nseq].count('-')
    Count = np.zeros((25,seqlenInt),dtype = float) + 0.000001
    gaps = 0   
   
    for i in range(0,seqlenAll):
            if Seq[nseq][i] == '-':
                gaps += 1
            else:
                for key in Seq.keys():
                    aa = Seq[key][i]
                    row = int(AA[aa])
                    col = i-gaps
                    Count[row,col] += 1
                    
    colSum = sum(Count[...,0])
    CountRate = Count/colSum
    Shannon = -(CountRate*np.log2(CountRate)).sum(axis=0)
    Rindex = math.log2(25)-Shannon
    Rrelat = (Rindex - Rindex.min())/(Rindex.max()-Rindex.min())
    Rquant = Quantile(Rindex)
    
# 储存目标氨基酸序列
    list0=[]
    faa=open(fasta_file,'r')
    faa1=faa.readlines()[1] 
    for i in list(faa1):
        if i!='\n':
            list0.append(i)
    faa.close()
    
    list1=[]
    for i in range(len(Rquant)):
        list1.append(str(Rquant[i]))    
    
    list1_int=[]
    for i in range(len(list1)):
        list1_int.append(float(list1[i]))
    
    no=[]
    for i in range(len(list0)):
        no.append(i+1)
    
    list2=ppssm(fasta_file,protein)
    list22=[]
    for i in range(len(list2)):
        list22.append(int(list2[i]))
        
    data=list(zip(no,list0,list1,list2))
    df=pd.DataFrame(data,columns=["pos","amino-acid","entropy","pssm-score"])
    xlsx_path=total_path+'/'+'con-'+protein+'_'+chain+'.xlsx'
    df.to_excel(xlsx_path,engine='openpyxl')

    ax=60
    conaa_pos=[]
    for i in range(len(list0)):
        if list1_int[i]>=ax and list22[i]>0:
            conaa_pos.append(str(i+1))

    path_con=total_path+"/"+'con-pos-'+protein+'_'+chain+'.txt'
    fx=open(path_con,"a")
    fx.write("\n".join(conaa_pos))
    fx.close()
   
if __name__ == '__main__':
    args1=en_index()
    chain0=args1.chain.split(sep='+')
    protein0=args1.pdb_name.split(sep='+')
    fasta_file0=args1.in_fasta_file.split(sep='+')
    way=args1.way
    if way=='msa':
        msa_file0=args1.msa_file.split(sep='+')
    for i in range(len(chain0)):
        protein=protein0[i]
        chain=chain0[i]
        fasta_file=fasta_file0[i]       
        path7="./step1-calcon-"+protein+"_"+chain
        os.makedirs(path7)
        
        if way=='msa':
            msa_file=msa_file0[i]
            ddata_msa(msa_file,fasta_file,path7,protein,chain)
        else:    
            os.system("blastp -query "+fasta_file+" -out ./"+protein+"_"+chain+".msa -db ./AutoDVD2/blast/db/swissprot/swissprot -outfmt 4 -evalue 1e-5")
            path6="./"+protein+"_"+chain+".msa"
            path5=path7+"/"+protein+"_"+chain+".msa"
            shutil.move(path6,path5)
            inputx=path7+"/"+protein+"_"+chain+".msa"
            output=path7+"/"+protein+"_"+chain+"-p1.msa"
            input1=output
            output1=path7+"/"+protein+"_"+chain+"-p2.msa"
            fpath=path7+"/split-"+protein+"_"+chain
            trans0(inputx,output)
            trans1(input1,output1)
            mkfile(output1,fpath)
            ddata(fasta_file,protein,fpath,chain)
    
    

import os
import sys
import math
import re
import random
import argparse
from pprint import pprint
import textwrap

def s3_index():
    parser = argparse.ArgumentParser(
        description='DNA疫苗自动化设计软件V2.0  AutoDVD2 --version 2.0',
        formatter_class=argparse.RawTextHelpFormatter,
        usage='type "python3 %(prog)s --help" for more information'
    )    
    parser.add_argument('-p', '--pdb_name', type=str, required=True,
                        help=textwrap.dedent('''\
    The user can enter the pathogen protein PDB number or the pathogen protein name. Note that the PDB number needs to be capitalized, 
    and the pathogen protein name must not contain spaces. 
    If you need to enter multiple names, please connect them with "+", for example: "5X5C+5K6G". '''))
    
    parser.add_argument('-ho', '--host', type=str, required=True,
                        help=textwrap.dedent('''\
    The "./AutoDVD2/codon" folder stores the codon usage tables of 22 common species. Users only need to enter the Latin scientific name of the host to be vaccinated. 
    Note that the scientific name space should be replaced by "-", such as "Sus-scrofa". If there is no appropiate codon usage table in the codon folder, users can obtain it 
    from http://www.kazusa.or.jp/codon/, format it according to files in codon folder, store it in the txt file, name it "<Latin scientific name (space is replaced by -)>.txt", 
    and put it into the codon folder. 
    The Latin scientific name and corresponding English name of the host that can be directly selected are as follows:
    1. Ailuropoda-melanoleuca --> panda            2. Anas-platyrhynchos --> duck                3. Anser-cygnoides-orientalis --> goose         
    4. Bos-taurus --> cattle                       5. Bubalus-bubalis --> water buffalo          6. Camelus-bactrianus --> Bactrian camel
    7. Camelus-dromedarius --> Dromedary camel     8. Canis-lupus --> dog                        9. Capra-hircus --> goat               
    10. Columba-livia --> dove                     11. Equus-asinus --> donkey                   12. Equus-caballus --> horse
    13. Felis-catus --> cat                        14. Gallus-gallus --> chicken                 15. Gorilla-gorilla --> gorilla                
    16. Homo-sapiens --> human being               17. Macaca-mulatta --> macaque                18. Mus-musculus --> mice
    19. Oryctolagus-cuniculus --> rabbit           20. Ovis-aries --> sheep                      21. Rattus-norvegicus --> rat                
    22. Sus scrofa --> pig'''))
    
    parser.add_argument('-e', '--enzymes', type=str, default="BamHI+NheI",
                        help='The name of the enzymes that selected in Double Digest should be connected with "+", and the default is "BamHI+NheI".')
                        
    parser.add_argument('-v', '--vector', type=str, default="./AutoDVD2/vector/pVAX1.fasta",
                        help='Users can upload carrier files by themselves. PVAX1 vector is used by default.')
                        
    parser.add_argument('-us', '--upstream_insert_sequence', type=str, default="GCCGCCACCGCCACCATG",
                        help=textwrap.dedent('''\
    Users can insert base sequence after upstream enzyme digestion site according to the selected carrier.
    "GCCGCCACCGCCACCATG" is used by default.'''))
                        
    parser.add_argument('-ds', '--downstream_insert_sequence', type=str, default="TAG",
                        help=textwrap.dedent('''\
    Users can insert base sequence before downstream enzyme digestion site according to the selected carrier.
    "TAG" is used by default.'''))
    
    args = parser.parse_args()
    pprint(args.__dict__)
    return args
    
def codon(co_table,x):
    fi=open(co_table,"r")
    b1=fi.readlines()
    c1=[]
    c2=[]
    c3=[]
    c4=[]
    for i in range(len(b1)):
        if len(b1[i])>1:
            a1=list(b1[i])
            a2=""
            a2+=a1[4]+a1[5]+a1[6]+a1[7]
            c1.append(float(a2))
    for i in range(len(b1)):
        if len(b1[i])>1:
            a1=list(b1[i])
            a2=""
            a2+=a1[22]+a1[23]+a1[24]+a1[25]
            c2.append(float(a2))
    for i in range(len(b1)):
        if len(b1[i])>1:
            a1=list(b1[i])
            a2=""
            a2+=a1[40]+a1[41]+a1[42]+a1[43]
            c3.append(float(a2))
    for i in range(len(b1)):
        if len(b1[i])>1:
            a1=list(b1[i])
            a2=""
            a2+=a1[58]+a1[59]+a1[60]+a1[61]
            c4.append(float(a2))
    fi.close()
    
    F0={"TTT":c1[0],"TTC":c1[1]}
    L0={"TTA":c1[2],"TTG":c1[3],"CTT":c1[4],"CTC":c1[5],"CTA":c1[6],"CTG":c1[7]}
    I0={"ATT":c1[8],"ATC":c1[9],"ATA":c1[10]}
    M0={"ATG":c1[11]}
    V0={"GTT":c1[12],"GTC":c1[13],"GTA":c1[14],"GTG":c1[15]}
    S0={"TCT":c2[0],"TCC":c2[1],"TCA":c2[2],"TCG":c2[3],"AGT":c4[8],"AGC":c4[9]}
    P0={"CCT":c2[4],"CCC":c2[5],"CCA":c2[6],"CCG":c2[7]}
    T0={"ACT":c2[8],"ACC":c2[9],"ACA":c2[10],"ACG":c2[11]}
    A0={"GCT":c2[12],"GCC":c2[13],"GCA":c2[14],"GCG":c2[15]}
    Y0={"TAT":c3[0],"TAC":c3[1]}
    H0={"CAT":c3[4],"CAC":c3[5]}
    Q0={"CAA":c3[6],"CAG":c3[7]}
    N0={"AAT":c3[8],"AAC":c3[9]}
    K0={"AAA":c3[10],"AAG":c3[11]}
    D0={"GAT":c3[12],"GAC":c3[13]}
    E0={"GAA":c3[14],"GAG":c3[15]}
    C0={"TGT":c4[0],"TGC":c4[1]}
    W0={"TGG":c4[3]}
    R0={"CGT":c4[4],"CGC":c4[5],"CGA":c4[6],"CGG":c4[7],"AGA":c4[10],"AGG":c4[11]}
    G0={"GGT":c4[12],"GGC":c4[13],"GGA":c4[14],"GGG":c4[15]}
    
    F1=sorted(F0.items(),key=lambda kv:kv[1],reverse=True)
    F=[F1[0][0],F1[1][0]]
    L1=sorted(L0.items(),key=lambda kv:kv[1],reverse=True)
    L=[L1[0][0],L1[1][0]]
    I1=sorted(I0.items(),key=lambda kv:kv[1],reverse=True)
    I=[I1[0][0],I1[1][0]]
    M1=sorted(M0.items(),key=lambda kv:kv[1],reverse=True)
    M=[M1[0][0]]
    V1=sorted(V0.items(),key=lambda kv:kv[1],reverse=True)
    V=[V1[0][0],V1[1][0]]
    S1=sorted(S0.items(),key=lambda kv:kv[1],reverse=True)
    S=[S1[0][0],S1[1][0]]
    P1=sorted(P0.items(),key=lambda kv:kv[1],reverse=True)
    P=[P1[0][0],P1[1][0]]
    T1=sorted(T0.items(),key=lambda kv:kv[1],reverse=True)
    T=[T1[0][0],T1[1][0]]
    A1=sorted(A0.items(),key=lambda kv:kv[1],reverse=True)
    A=[A1[0][0],A1[1][0]]
    Y1=sorted(Y0.items(),key=lambda kv:kv[1],reverse=True)
    Y=[Y1[0][0],Y1[1][0]]
    H1=sorted(H0.items(),key=lambda kv:kv[1],reverse=True)
    H=[H1[0][0],H1[1][0]]
    Q1=sorted(Q0.items(),key=lambda kv:kv[1],reverse=True)
    Q=[Q1[0][0],Q1[1][0]]
    N1=sorted(N0.items(),key=lambda kv:kv[1],reverse=True)
    N=[N1[0][0],N1[1][0]]
    K1=sorted(K0.items(),key=lambda kv:kv[1],reverse=True)
    K=[K1[0][0],K1[1][0]]
    D1=sorted(D0.items(),key=lambda kv:kv[1],reverse=True)
    D=[D1[0][0],D1[1][0]]
    E1=sorted(E0.items(),key=lambda kv:kv[1],reverse=True)
    E=[E1[0][0],E1[1][0]]
    C1=sorted(C0.items(),key=lambda kv:kv[1],reverse=True)
    C=[C1[0][0],C1[1][0]]
    W1=sorted(W0.items(),key=lambda kv:kv[1],reverse=True)
    W=[W1[0][0]]
    R1=sorted(R0.items(),key=lambda kv:kv[1],reverse=True)
    R=[R1[0][0],R1[1][0]]
    G1=sorted(G0.items(),key=lambda kv:kv[1],reverse=True)
    G=[G1[0][0],G1[1][0]]
    
    
    y=""
    for i in range(len(x)):
        if x[i]=="F":
            z1=random.choice([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])
            if z1==0:
                y+=F[0]
            else:
                y+=F[1]
        elif x[i]=="L":
            z2=random.choice([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])
            if z2==0:
                y+=L[0]
            else:
                y+=L[1]
        elif x[i]=="I":
            z2=random.choice([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])
            if z2==0:
                y+=I[0]
            else:
                y+=I[1]
        elif x[i]=="M":
            y+=M[0]  
        elif x[i]=="V":
            z2=random.choice([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])
            if z2==0:
                y+=V[0]
            else:
                y+=V[1]
        elif x[i]=="S":
            z2=random.choice([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])
            if z2==0:
                y+=S[0]
            else:
                y+=S[1]  
        elif x[i]=="P":
            z2=random.choice([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])
            if z2==0:
                y+=P[0]
            else:
                y+=P[1]        
        elif x[i]=="T":
            z2=random.choice([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])
            if z2==0:
                y+=T[0]
            else:
                y+=T[1]
        elif x[i]=="A":
            z2=random.choice([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])
            if z2==0:
                y+=A[0]
            else:
                y+=A[1]
        elif x[i]=="Y":
            z1=random.choice([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])
            if z1==0:
                y+=Y[0]
            else:
                y+=Y[1] 
        elif x[i]=="H":
            z1=random.choice([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])
            if z1==0:
                y+=H[0]
            else:
                y+=H[1] 
        elif x[i]=="Q":
            z1=random.choice([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])
            if z1%2==0:
                y+=Q[0]
            else:
                y+=Q[1]            
        elif x[i]=="N":
            z1=random.choice([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])
            if z1%2==0:
                y+=N[0]
            else:
                y+=N[1]    
        elif x[i]=="K":
            z1=random.choice([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])
            if z1%2==0:
                y+=K[0]
            else:
                y+=K[1]          
        elif x[i]=="D":
            z1=random.choice([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])
            if z1%2==0:
                y+=D[0]
            else:
                y+=D[1]       
        elif x[i]=="E":
            z1=random.choice([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])
            if z1%2==0:
                y+=E[0]
            else:
                y+=E[1]       
        elif x[i]=="C":
            z1=random.choice([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])
            if z1%2==0:
                y+=C[0]
            else:
                y+=C[1]       
        elif x[i]=="W":
            y+=W[0]
        elif x[i]=="R":
            z2=random.choice([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])
            if z2==0:
                y+=R[0]
            else:
                y+=R[1]       
        elif x[i]=="G":
            z2=random.choice([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1])
            if z2==0:
                y+=G[0]
            else:
                y+=G[1]
    return y

def trans(co_table,pros,pathaa):   
    f1=open("./AutoDVD2/5FS4.fasta","r")
    a1=f1.readlines()[1]  
    a11=list(a1)
    b=len(a11)-1
    ap205=a11[0:b]
    f1.close()
     
    path0='./step2-ellipro-'+pros+'/'+pros+'_filter-epitopes.txt'
    f2=open(path0,"r")
    a4=f2.readlines()
    for j in range(len(a4)):
        a6=['G','G','G','S']
        a44=list(a4[j])
        for m in range(len(a6)):
            ap205.append(a6[m])
        for k in range(len(a44)-1):   
            ap205.append(a44[k])
    f2.close()

    a3='>AP205'+pros+'_aa\n'
    for i in range(len(ap205)):
        a3+=ap205[i]
    fk=open(pathaa,'a')
    fk.write(a3)
    fk.close()
    xx=codon(co_table,ap205)
    return xx
       
def exclude(co_table,x,xx,en1,en2):
    d={'BamHI':['GGATCC','CCTAGG'],'NheI':['GCTAGC','CGATCG'],'AccI':['GTCGAC','CAGCTG'],'AflIII':['ACATGT','TGTACA'],'AscI':['GGCGCGCC','CCGCGCGG'],
    'AvaI':['CCCGGG','GGGCCC'],'BglII':['AGATCT','TCTAGA'],'BssHII':['GCGCGC','CGCGCG'],'BstXI':['CCAATGCATTGG','GGTTACGTAACC'],'ClaI':['ATCGAT','TAGCTA'],'EcoRI':['GAATTC','CTTAAG'],
    'HaeIII':['GGCC','CCGG'],'HindIII':['AAGCTT','TTCGAA'],'KpnI':['GGTACC','CCATGG'],'MluI':['ACGCGT','TGCGCA'],'NcoI':['CCATGG','GGTACC'],
    'NdeI':['CATATG','GTATAC'],'NotI':['GCGGCCGC','CGCCGGCG'],'NsiI':['ATGCAT','TACGTA'],'PacI':['TTAATTAA','AATTAATT'],'PmeI':['GTTTAAAC','CAAATTTG'],
    'PstI':['CTGCAG','GACGTC'],'PvuI':['CGATCG','GCTAGC'],'SacI':['GAGCTC','CTCGAG'],'SacII':['CCGCGG','GGCGCC'],'SalI':['GTCGAC','CAGCTG'],'ScaI':['AGTACT','TCATGA'],
    'SmaI':['CCCGGG','GGGCCC'],'SpeI':['ACTAGT','TGATCA'],'SphI':['GCATGC','CGTACG'],'StuI':['AGGCCT','TCCGGA'],'XbaI':['TCTAGA','AGATCT'],'XhoI':['CTCGAG','GAGCTC'],'XmaI':['CCCGGG','GGGCCC']
    }
    x1=[]
    when1=0
    when2=0
    co=0
    
    ens=[]
    for i in range(len(d)):
        ens.append(i)

    seqs=list(d.values())
    for i,j in d.items():
        if en1==i:
            when1=ens[co]
        if en2==i:
            when2=ens[co]
        co+=1
                
    c1=xx.count(seqs[when1][0])
    c2=xx.count(seqs[when1][1])
    c3=xx.count(seqs[when2][0])
    c4=xx.count(seqs[when2][1])
    k=0
    while c1!=0 or c2!=0 or c3!=0 or c4!=0:
        x1.append(codon(co_table,x))
        c1=x1[k].count(seqs[when1][0])
        c2=x1[k].count(seqs[when1][1])
        c3=x1[k].count(seqs[when2][0])
        c4=x1[k].count(seqs[when2][1])
        k+=1

    if len(x1)==0:
        return xx            
    else:
        return x1[-1]
    
  
def connect(seq,vec_file,en11,en22,pros,up,down):
    seq1=up
    seq1+=seq+down
    f1=open(vec_file,'r')
    a1=f1.readlines()[1]
    f1.close()
    a2=list(a1)
    d={'BamHI':['GGATCC','CGCGGATCCGCG'],'NheI':['GCTAGC','CTAGCTAGCTAG'],'AccI':['GTCGAC','CCGGTCGACCGG'],'AflIII':['ACATGT','CCCACATGTGGG'],'AscI':['GGCGCGCC','TTGGCGCGCCAA'],
    'AvaI':['CCCGGG','TCCCCCGGGGGA'],'BglII':['AGATCT','GAAGATCTTC'],'BssHII':['GCGCGC','TTGGCGCGCCAA'],'BstXI':['CCAATGCATTGG','CTGCAGAACCAATGCATTGGATGCAT'],'ClaI':['ATCGAT','CCATCGATGG'],
    'EcoRI':['GAATTC','CCGGAATTCCGG'],
    'HaeIII':['GGCC','TTGCGGCCGCAA'],'HindIII':['AAGCTT','CCCAAGCTTGGG'],'KpnI':['GGTACC','CGGGGTACCCCG'],'MluI':['ACGCGT','CGACGCGTCG'],'NcoI':['CCATGG','CATGCCATGGCATG'],
    'NdeI':['CATATG','GGAATTCCATATGGAATTCC'],'NotI':['GCGGCCGC','AAGGAAAAAAGCGCCGCAAAAGGAAAA'],'NsiI':['ATGCAT','CCAATGCATTGGTTCTGCAGTT'],'PacI':['TTAATTAA','CCTTAATTAAGG'],
    'PmeI':['GTTTAAAC','AGCTTTGTTTAAACGGCGCGCCGG'],
    'PstI':['CTGCAG','AACTGCAGAACCAATGCATTGG'],'PvuI':['CGATCG','ATCGATCGAT'],'SacI':['GAGCTC','CGAGCTCG'],'SacII':['CCGCGG','TCCCCGCGGGGA'],
    'SalI':['GTCGAC','ACGCGTCGACGTCGGCCATAGCGGCCGCGGAA'],'ScaI':['AGTACT','AAAAGTACTTTT'],
    'SmaI':['CCCGGG','TCCCCCGGGGGA'],'SpeI':['ACTAGT','GACTAGTC'],'SphI':['GCATGC','ACATGCATGCATGT'],'StuI':['AGGCCT','AAGGCCTT'],'XbaI':['TCTAGA','GCTCTAGAGC'],
    'XhoI':['CTCGAG','CCGCTCGAGCGG'],'XmaI':['CCCGGG','TCCCCCCGGGGGGA']
    }
    x1=[]
    pen1=0
    pen2=0
    co=0
    
    ens=[]      
    for i in range(len(d)):
        ens.append(i)

    seqs=list(d.values())         
    for i,j in d.items():
        if en11==i:
            pen11=ens[co]
        if en22==i:
            pen22=ens[co]
        co+=1
    
    pos1=a1.find(seqs[pen11][0])   
    pos2=a1.find(seqs[pen22][0])
    len1=len(seqs[pen11][0])
    len2=len(seqs[pen22][0])
    enpro1=seqs[pen11][1]
    enpro2=seqs[pen22][1]
    
    if pos1<pos2:
        xpos1=pos1+len1
        xpos2=pos2
        xenpro1=enpro1
        xenpro2=enpro2
    else:
        xpos1=pos2+len2
        xpos2=pos1
        xenpro1=enpro2
        xenpro2=enpro1
    
    sseq='>AP205'+pros+'-synthesis_dna\n'+xenpro1+seq1+xenpro2
    spath="./step3-connection-"+pros+'/'+pros+'_synthesis.fasta'
    fs=open(spath,'a')
    fs.write(sseq)
    fs.close()
       
    b1=a2[:xpos1]
    b2=a2[xpos2:len(a2)-1]
    c1=''
    c2=''
    for i in range(len(b1)):
        c1+=b1[i]
    for i in range(len(b2)):
        c2+=b2[i]
    plasmid=c1+seq1+c2
    
    return plasmid
              
if __name__ == '__main__':
    arg3=s3_index()
    host=arg3.host
    pros=arg3.pdb_name
    vector=arg3.vector
    up=arg3.upstream_insert_sequence
    down=arg3.downstream_insert_sequence
    enzyme0=arg3.enzymes.split(sep='+')
    co_table="./AutoDVD2/codon/"+host+".txt"

    path1="./step3-connection-"+pros
    os.makedirs(path1)
    path2=path1+"/"+pros+"_ap205-epitopes-aa.fasta"
    path3=path1+"/"+pros+"_ap205-epitopes-dna.fasta"
    path5=path1+'/'+pros+'_final-plasmid.fasta'
    
    xx=trans(co_table,pros,path2)
    fi=open(path2,'r')
    rd=fi.readlines()
    rd1=rd[1]
    x=''
    for i in range(len(rd1)):
        x+=rd1[i]
        
    en1=enzyme0[0]
    en2=enzyme0[1]
    xx2=exclude(co_table,x,xx,en1,en2)
    
    f3=open(path3,"a")
    xxx2='>AP205'+pros+'_dna\n'+xx2
    f3.write(xxx2)
    f3.close()    

    pla=connect(xx2,vector,en1,en2,pros,up,down)
    fpla='>AP205'+pros+'_final-plasmid_dna\n'+pla
    f5=open(path5,'a')
    f5.write(fpla)
    f5.close()
    
    
    
  

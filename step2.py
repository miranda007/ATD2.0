import sys
import os
import pandas as pd
import re
import xlrd
import numpy as np
import shutil
import argparse
from pprint import pprint

def el_index():
    parser = argparse.ArgumentParser(
        description='Process introduction.',
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument('-c', '--chain', type=str, required=True,
                        help='链(chain)。')
    parser.add_argument('-p', '--pdb_name', type=str, required=True,
                        help='病毒蛋白PDB号。')
    parser.add_argument('-ip', "--in_PDB_FILE", type=str, required=True,
                        help='PDB 文件路径。')                    
    args = parser.parse_args()
    pprint(args.__dict__)
    return args

def ell(protein,chain,pdb_file):
    initial_epitopes_file0=folder2+"/"+protein+"_"+chain+"-epitopes.txt"
    ellipro_table_file0=folder2+"/"+protein+"_"+chain+"-table.txt"
    os.system("java -jar ./ATD2.0/ElliPro.jar -f "+pdb_file+' -c '+chain+" --output "+initial_epitopes_file0+" --table "+ellipro_table_file0)
    
    sz=os.path.getsize(initial_epitopes_file0)
    if not sz:
        print("\nWarning: No epitopes in ",protein,' ',chain,"!")
        return
        
    inn2="./step1-calcon-"+protein+"_"+chain+"/con-pos-"+protein+"_"+chain+".txt"
    inn3="./step1-calcon-"+protein+"_"+chain+"/split-"+protein+"_"+chain

# 读取氨基酸所有序列    
    list_aa_seq=[]                                                    
    str_aa_seq=""
    num1=len(os.listdir(inn3))
    for i in range(num1):
        in_msa=open(inn3+"/"+str(i)+".msa.fasta","r")
        nd=in_msa.readlines()[1]
        nd1=nd.rstrip()
        for i in list(nd1):
            if i !="-":
                list_aa_seq.append(i)
        in_msa.close()  
    
    for i in range(len(list_aa_seq)):
        str_aa_seq+=str(list_aa_seq[i])
        
    fo=open(initial_epitopes_file0,"r")
    initial_eps=[]
    single_con_aa_inep=[]
    continue_con_aa_pos=[]
    in_ep_pos=[]
    
# 通过ellipro文件“，”位置获得表位序列    
    for line in fo:
        doc_pos=[]
        L=list(line)
        if "Type" in line:
            continue
        elif len(L)<=1:
            break
        else:
            for i in range(len(L)):
                if L[i]==",":
                    doc_pos.append(i)
        initial_ep="".join(L[(doc_pos[4]+1):doc_pos[5]])
        initial_eps.append(initial_ep)
    fo.close()

# 从氨基酸序列中找表位位置    
    for i in range(len(initial_eps)):
        if re.search(initial_eps[i],str_aa_seq)!=None:
            index1=re.search(initial_eps[i],str_aa_seq).span()
            index2=int(index1[0])
            in_ep_pos.append(index2+1)
            in_ep_pos.append(index2+len(initial_eps[i]))

# 找到表位中的保守性氨基酸位置         
    fi=open(inn2,"r")
    con_aa_pos=fi.readlines()
    for i in range(len(con_aa_pos)):
        for j in range(int(len(in_ep_pos)*0.5)):
            if 2*j==len(in_ep_pos):
                break
            if int(in_ep_pos[2*j])<int(con_aa_pos[i]) and int(con_aa_pos[i])<int(in_ep_pos[2*j+1]):
                x1=con_aa_pos[i].replace("\n","")
                single_con_aa_inep.append(x1)
    fi.close()

# 寻找表位中的连续的保守性氨基酸    
    split_continue= []
    for i in range(len(single_con_aa_inep)):
        if not split_continue:
            split_continue.append([single_con_aa_inep[i]])
        elif int(single_con_aa_inep[i - 1]) + 1 == int(single_con_aa_inep[i]):
            split_continue[-1].append(single_con_aa_inep[i])
        else:
            split_continue.append([single_con_aa_inep[i]])

    for i in range(len(split_continue)):
         if len(split_continue[i])>=4 and len(split_continue[i])<=10:
             continue_con_aa_pos.append(split_continue[i])
    lw=len(continue_con_aa_pos)
             
              
    short_eps_no=[]
    start=[]
    end=[]
    sh_eplen=[]
    shep_seq=[]
    shep_conscore=[]
    for i in range(lw):
        short_eps_no.append(i+1)
        start.append(continue_con_aa_pos[i][0])
        sh_ep_len=len(continue_con_aa_pos[i])
        end.append(continue_con_aa_pos[i][sh_ep_len-1])
        sh_eplen.append(sh_ep_len)
        aaa=int(continue_con_aa_pos[i][0])-1
        bbb=int(continue_con_aa_pos[i][sh_ep_len-1])
        ab="".join(list_aa_seq[aaa:bbb])
        shep_seq.append(ab)
    
    if len(shep_seq)==0:
        print("\nWarning: No conserved epitopes in ",protein,' ',chain," !")
        return


# 计算表位平均保守性得分        
    path="./step1-calcon-"+protein+"_"+chain+"/con-"+protein+"_"+chain+".xls"
    fi=xlrd.open_workbook(path)
    fo=fi.sheets()[0]
    col_entro=fo.col_values(3)
    str_all_entro=col_entro[1:]
    float_all_entro=[]
    for i in range(len(str_all_entro)):
        float_all_entro.append(float(str_all_entro[i]))
    for i in range(len(start)):
        op1=int(start[i])-1
        op2=int(end[i])
        c=sum(float_all_entro[op1:op2])/int(sh_eplen[i])
        cc=np.around(c,decimals=2)
        shep_conscore.append(cc)
    
    lipdb=[]
    lichain=[]
    for i in range(len(start)):
       lipdb.append(protein)
       lichain.append(chain)
       
# 计算短表位平均抗原性得分
    shep_antiscore=[]
    ftable=open(ellipro_table_file0,'r')
    full_table=ftable.readlines()
    meaning_lines=full_table[1:]
    aa3_intable=[]
    anti_score=[]
    aa1_intable=''
    for i in range(len(meaning_lines)):
        line_doc=[]
        list_perline=list(meaning_lines[i])
        if meaning_lines[i].find('UNK')!=-1:
            continue
        if len(list_perline)<=1:
            break
        for k in range(len(list_perline)):            
            if list_perline[k]==',':
                line_doc.append(k)            

        aa3_inline=list_perline[(line_doc[2]+1)]+list_perline[(line_doc[2]+2)]+list_perline[(line_doc[2]+3)]
        aa3_intable.append(aa3_inline)
        score_perline=list_perline[(line_doc[3]+1)]+list_perline[(line_doc[3]+2)]+list_perline[(line_doc[3]+3)]+list_perline[(line_doc[3]+4)]+list_perline[(line_doc[3]+5)]
        flscore_perline=float(score_perline)
        anti_score.append(flscore_perline)  
        
    aa3=['ALA','ARG','ASP','CYS','GLN','GLU','HIS','ILE','GLY','ASN','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    aa1=['A','R','D','C','Q','E','H','I','G','N','L','K','M','F','P','S','T','W','Y','V']
    for i in range(len(aa3_intable)):
        for j in range(len(aa3)):
            if aa3_intable[i]==aa3[j]:
                aa1_intable+=aa1[j]
    
    for i in range(len(shep_seq)):
        ep_tapos=aa1_intable.find(shep_seq[i])
        leni=len(shep_seq[i])
        ave_anti0=sum(anti_score[ep_tapos:(ep_tapos+leni)])/leni
        ave_anti1=np.around(ave_anti0,decimals=4)
        if ave_anti1>0.5:
            shep_antiscore.append(ave_anti1)
        else:
            del short_eps_no[i]
            del lipdb[i]
            del lichain[i]
            del start[i]
            del end[i]
            del shep_seq[i]
            del sh_eplen[i]
            del shep_conscore[i]
    ftable.close()   
    
    return short_eps_no,lipdb,lichain,start,end,shep_seq,sh_eplen,shep_antiscore,shep_conscore         
   
        
if __name__ == '__main__':
    args=el_index()
    pros=args.pdb_name
    folder2="./step2-ellipro-"+pros
    os.makedirs(folder2)
    chain0=args.chain.split(sep='+')
    protein0=args.pdb_name.split(sep='+')
    pdb_file0=args.in_PDB_FILE.split(sep='+')   
    short_eps_nox=[]
    lipdbx=[]
    lichainx=[]
    startx=[]
    endx=[]
    shep_seqx=[]
    sh_eplenx=[]
    shep_antiscorex=[]
    shep_conscorex=[]
    for i in range(len(chain0)):
        tup=ell(protein0[i],chain0[i],pdb_file0[i])
        if tup==None:
            continue
        for j in range(len(tup[0])):
            short_eps_nox.append(tup[0][j])
            lipdbx.append(tup[1][j])
            lichainx.append(tup[2][j])
            startx.append(tup[3][j])
            endx.append(tup[4][j])
            shep_seqx.append(tup[5][j])
            sh_eplenx.append(tup[6][j])
            shep_antiscorex.append(tup[7][j])
            shep_conscorex.append(tup[8][j])
    
    filter_scorex=[]
    for i in range(len(lipdbx)):
        fs=sum([0.65*shep_antiscorex[i],0.0035*shep_conscorex[i]])
        filter_scorex.append(fs)  
    
    if len(lipdbx)>1:        
        path11=folder2+"/"+pros+"-filter.xls"
        data1=list(zip(short_eps_nox,lipdbx,lichainx,startx,endx,shep_seqx,sh_eplenx,shep_antiscorex,shep_conscorex,filter_scorex))
        df=pd.DataFrame(data1,columns=["No.",'PDB','Chain',"start","end","sequence","length",'ave-PI',"ave-entropy",'filter-score'])
        df.to_excel(path11)

        total_len=sum(sh_eplenx)
        while total_len>(39-4*len(lipdbx)):   
            b=sorted(enumerate(filter_scorex),key=lambda x:x[1],reverse=False)
            xiabiao=[i[0] for i in b]
            paixu=[i[1] for i in b]
            del shep_seqx[xiabiao[0]]
            del sh_eplenx[xiabiao[0]]
            del lipdbx[xiabiao[0]]
            del filter_scorex[xiabiao[0]]
            total_len=sum(sh_eplenx) 
               
        pathx3=folder2+"/"+pros+"_filter-epitopes.txt"
        fx1=open(pathx3,"a")
        for i in range(len(shep_seqx)):
            fx1.write(str(shep_seqx[i])+"\n")
        fx1.close()



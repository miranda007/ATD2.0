import os
import argparse
from datetime import datetime as dt
from pprint import pprint
import shutil

def all_index():
    parser = argparse.ArgumentParser(
        description='Process introduction.',
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument('-if', "--in_FASTA_FILE", type=str, required=True,
                        help='Before uploading, please check whether the fasta file only contains sequences of one specific chain. If not, you need to manually delete redundant amino acid sequences. If you need to input multiple files, the absolute path of the file should be connected with "+" and capitalized, for example: "/ ATD2.0/demo/5X5C.fasta+./ATD2.0/demo/5K6G.fasta”。 Note that the writing sequence of the document and PDB number should correspond.')
    parser.add_argument('-ip', "--in_PDB_FILE", type=str, required=True,
                        help='If you need to input multiple files, the absolute path of the file should be connected with "+" and capitalized, for example: "/ ATD2.0/demo/5X5C.pdb+./ATD2.0/demo/5K6G.pdb”。 Note that the writing sequence of the document and PDB number should correspond.')                    
    parser.add_argument('-c', '--chain', type=str, required=True,
                        help='If you need to input multiple pathogen proteins, or select multiple chains for the same pathogen protein, the selected chains should be connected with "+" and capitalized, such as "A+A". Note that the writing sequence of the chain and PDB number should correspond.')
    parser.add_argument('-p', '--pdb_name', type=str, required=True,
                        help='If you need to input multiple pathogen proteins, or select multiple chains for the same pathogen protein, the PDB number should be connected with "+" and capitalized, for example, "5X5C+5K6G".')
    parser.add_argument('-e', '--enzymes', type=str, default="BamHI+NheI",
                        help='The name of the enzymes that selected in Double Digest should be connected with "+", and the default is "BamHI+NheI".')
    parser.add_argument('-ho', '--host', type=str, required=True,
                        help='The "./ATD2.0/codon" folder stores the codon usage tables of 22 common species. Users only need to enter the Latin scientific name of the host to be vaccinated. Note that the scientific name space should be replaced by "-", such as "Sus-scrofa". If there is no appropiate codon usage tables in the codon folder, users can download them from http://www.kazusa.or.jp/codon/, name it "<Latin scientific name (space is replaced by -)>. txt", and put it in the codon folder.')
    parser.add_argument('-v', '--vector', type=str, default="./ATD2.0/vector/pVAX1.fasta",
                        help='Users can upload carrier files by themselves. PVAX1 vector is used by default.')                    
    args = parser.parse_args()
    pprint(args.__dict__)
    return args
          
if __name__ == '__main__':
    all_args=all_index()
    os.system('python3 ./ATD2.0/step1.py -c '+all_args.chain+' -if '+all_args.in_FASTA_FILE+' -p '+all_args.pdb_name)
    os.system('python3 ./ATD2.0/step2.py -c '+all_args.chain+' -ip '+all_args.in_PDB_FILE+' -p '+all_args.pdb_name)
    if os.path.exists('./step2-ellipro-'+all_args.pdb_name+'/'+all_args.pdb_name+'_filter-epitopes.txt')==True:
        os.system('python3 ./ATD2.0/step3.py -p '+all_args.pdb_name+' -e '+all_args.enzymes+' -ho '+all_args.host+' -v '+all_args.vector)
    else:
        print('No appropiate epitopes!')
    
    path_sum='./ATD2.0-'+all_args.pdb_name+'_'+all_args.chain
    os.makedirs(path_sum)
    chain0=all_args.chain.split(sep='+')
    protein0=all_args.pdb_name.split(sep='+')
    for i in range(len(chain0)):
        old1='./step1-calcon-'+protein0[i]+'_'+chain0[i]
        new1=path_sum+'/step1-calcon-'+protein0[i]+'_'+chain0[i]
        shutil.move(old1,new1)
    old2='./step2-ellipro-'+all_args.pdb_name
    new2=path_sum+'/step2-ellipro-'+all_args.pdb_name
    shutil.move(old2,new2)
    old3='./step3-connection-'+all_args.pdb_name
    new3=path_sum+'/step3-connection-'+all_args.pdb_name
    shutil.move(old3,new3)

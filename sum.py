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
                        help='FASTA 文件路径。')
    parser.add_argument('-ip', "--in_PDB_FILE", type=str, required=True,
                        help='PDB 文件路径。')                    
    parser.add_argument('-c', '--chain', type=str, required=True,
                        help='链(chain)。')
    parser.add_argument('-p', '--pdb_name', type=str, required=True,
                        help='病毒蛋白PDB号。')
    parser.add_argument('-e', '--enzymes', type=str, default="BamHI+NheI",
                        help='双酶切法所选择的酶名称。')
    parser.add_argument('-ho', '--host', type=str, required=True,
                        help='宿主拉丁文学名。')
    parser.add_argument('-v', '--vector', type=str, default="./ATD2.0/vector/pVAX1.fasta",
                        help='载体文件路径。')                    
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
    
    os.makedirs('./ATD2.0-'+all_args.pdb_name)
    chain0=all_args.chain.split(sep='+')
    protein0=all_args.pdb_name.split(sep='+')
    for i in range(len(chain0)):
        old1='./step1-calcon-'+protein0[i]+'_'+chain0[i]
        new1='./ATD2.0-'+all_args.pdb_name+'/step1-calcon-'+protein0[i]+'_'+chain0[i]
        shutil.move(old1,new1)
    old2='./step2-ellipro-'+all_args.pdb_name
    new2='./ATD2.0-'+all_args.pdb_name+'/step2-ellipro-'+all_args.pdb_name
    shutil.move(old2,new2)
    old3='./step3-connection-'+all_args.pdb_name
    new3='./ATD2.0-'+all_args.pdb_name+'/step-connection-'+all_args.pdb_name
    shutil.move(old3,new3)

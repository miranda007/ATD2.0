import os
import argparse
from datetime import datetime as dt
from pprint import pprint
import shutil
import textwrap

def all_index():
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
    
    parser.add_argument('-ip', "--in_PDB_file", type=str, required=True,
                        help=textwrap.dedent('''\
    Must-option. 
    The name format of the pdb file should be "<capitalized PDB number or pathogen protein name>.pdb". 
    If you need to input multiple files, the absolute paths of the file should be connected with "+", 
    for example: "./AutoDVD2/demo/5X5C.pdb+./AutoDVD2/demo/5K6G.pdb”. 
    Note that the writing order of the paths should correspond to the PDB number or the name of the pathogen protein.'''))                    

    parser.add_argument('-ho', '--host', type=str, required=True,
                        help=textwrap.dedent('''\
    Must-option.
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
          
if __name__ == '__main__':
    all_args=all_index()
    os.system('python3 ./AutoDVD2/step1.py -c '+all_args.chain+' -if '+all_args.in_fasta_file+' -p '+all_args.pdb_name)
    os.system('python3 ./AutoDVD2/step2.py -c '+all_args.chain+' -ip '+all_args.in_PDB_file+' -p '+all_args.pdb_name)
    if os.path.exists('./step2-ellipro-'+all_args.pdb_name+'/'+all_args.pdb_name+'_filter-epitopes.txt')==True:
        os.system('python3 ./AutoDVD2/step3.py -p '+all_args.pdb_name+' -e '+all_args.enzymes+' -ho '+all_args.host+' -v '+all_args.vector
        +' -us '+all_args.upstream_insert_sequence+' -ds '+all_args.downstream_insert_sequence)
    else:
        print('No appropiate epitopes!')
    
    path_sum='./AutoDVD2-'+all_args.pdb_name+'_'+all_args.chain
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

# ATD2.0
AutoDVD2

## Author
Hongli Du; Jerome Lon; Xinwen Zheng

## Email
1602371546@qq.com

# Introduction
AutoDVD2 is a Linux-based command-line-driven tool for automating the construction of DNA vaccines. AutoDVD2 can quickly construct DNA vaccine plasmids that can express conserved linear B cell epitopes. AutoDVD2 follows the basic construction framework of AutoDVD, but has improved the amino acid conservation analysis, antigen epitope screening and codon optimization. Both AutoDVD conservative analysis and codon optimization need to be connected to the Internet to capture the results of web page analysis. AutoDVD2 uses a new set of algorithms without networking, which greatly improves the stability and speed of the software.It integrates automated methods into an easy-to-use tool to assist researchers without programming skills.

The main functions of the software are as follows:
1. Linear B cell epitope screening
2. Epitope Conservation Analysis
3. Codon Optimization
4. Insertion sequence design suitable for double digestion vector construction
5. DNA Vaccine Plasmid Construction

# Installation
## 1.Obtain AutoDVD2 source code:  
1) Register github account;  
2) Run the following command line at the terminal:  
```    
sudo apt install git   
```
3) Download Git Credential Manager at: https://github.com/GitCredentialManager/git-credential-manager/blob/release/docs/install.md  
   Select Tarball install at: https://github.com/GitCredentialManager/git-credential-manager/releases/download/v2.0.886/gcm-linux_amd64.2.0.886.tar.gz  
4) Run the following command line at the terminal:  
```  
sudo tar -xvf <path-to-tarball> -C /usr/local/bin  
sudo git-credential-manager configure  
sudo git config --global credential.credentialStore secretservice  
sudo git config --global credential.credentialStore gpg  
sudo git config --global credential.credentialStore cache  
sudo git config --global credential.cacheOptions "--timeout 300"  
sudo git config --global credential.credentialStore plaintext  
sudo git config --global user.email <email account>  
sudo git config --global user.name “<user name of github>”  
sudo git config --global --unset http.proxy  
sudo git config --global --unset https.proxy  
sudo git clone https://github.com/miranda007/AutoDVD2.git  
```
## 2.Install dependent environment:
Run the following command line at the terminal:  
```
sudo apt install openjdk-11-jre-headless  
sudo apt install pip  
pip install -r ./AutoDVD2/requirements.txt  
```
## 3.Install blast:
1) Download blast at: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz  
2) Run the following command line at the terminal:  
```
sudo chown -R $USER /home/$USER
sudo chown -R $USER <path-to-blast>
sudo tar -xvf <path-to-blast>
sudo chown -R $USER ./ncbi-blast-2.13.0+
sudo mv ./ncbi-blast-2.13.0+ ./AutoDVD2/blast
sudo vim ~/.bashrc
```
3) Open ~/.bashrc, insert " export PATH=./AutoDVD2/blast/bin:$PATH" at the last line.  
4) Run the following command line at the terminal:  
```
source ~/.bashrc
```
## 4.Install swissprot:
1) Download swissprot at: https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz  
2) Run the following command line at the terminal:  
```
sudo mkdir ./AutoDVD2/blast/db
sudo chown -R $USER ./AutoDVD2/blast/db
sudo mkdir ./AutoDVD2/blast/db/swissprot
sudo chown -R $USER ./AutoDVD2/blast/db/swissprot
cd ./AutoDVD2/blast/db/swissprot
sudo tar -xvf <path-to-swissprot>
```

# Test  
(Enter the following command to find the relevant output file in the $HOME/$USER directory)  
```
python3 ./AutoDVD2/sum.py -c A+A -p 5X5C+5K6G -if ./AutoDVD2/demo/5X5C.fasta+./AutoDVD2/demo/5K6G.fasta -ip ./AutoDVD2/demo/5X5C.pdb+./AutoDVD2/demo/5K6G.pdb -ho Sus-scrofa
```
# Usage
Run the following command line at the terminal:  
```
python3 ./AutoDVD2/sum.py <parameters>
```

# Parameters
(Items 1-5 are mandatory parameters, and items 6-9 are optional parameters)  
1.``` -p   --pdb_name```    
  The user can enter the pathogen protein PDB number or the pathogen protein name. Note that the PDB number needs to be capitalized, and the pathogen protein name must not contain spaces. If you need to enter multiple names, please connect them with "+", for example: "5X5C+5K6G".  
  
2. ```-c   --chain```   
  Note that the name of the selected chain needs to be capitalized. If you need to enter multiple names, you need to connect them with "+", for example: "A+A". Note that the writing order of the chain and PDB number or pathogen protein name should correspond.  
  
3.``` -if   --in_fasta_file```  
  The name format of the fasta file should be "<capitalized PDB number or pathogen protein name>.fasta". Before uploading, please check whether the fasta file only contains sequences of one specific chain.If not, you need to manually delete redundant amino acid sequences. If you need to input multiple files, the absolute paths of the file should be connected with "+", for example: "./AutoDVD2/demo/5X5C.fasta+./AutoDVD2/demo/5K6G.fasta”. Note that the writing order of the paths should correspond to the PDB number or the name of the pathogen protein.  
  
4.``` -ip   --in_PDB_file```    
  The name format of the pdb file should be "<capitalized PDB number or pathogen protein name>.pdb". If you need to input multiple files, the absolute paths of the file should be connected with "+", for example: "./AutoDVD2/demo/5X5C.pdb+./AutoDVD2/demo/5K6G.pdb”. Note that the writing order of the paths should correspond to the PDB number or the name of the pathogen protein.  
  
5.``` -ho   --host```  
  The "./AutoDVD2/codon" folder stores the codon usage tables of 22 common species. Users only need to enter the Latin scientific name of the host to be vaccinated. Note that the scientific name space should be replaced by "-", such as "Sus-scrofa". If there is no appropiate codon usage table in the codon folder, users can obtain it from http://www.kazusa.or.jp/codon/, format it according to files in codon folder, store it in the txt file, name it "Latin scientific name (space is replaced by -).txt",  and put it into the codon folder.   
  The Latin scientific name and corresponding English name of the host that can be directly selected are as follows:  
1. Ailuropoda-melanoleuca --> panda             2. Anas-platyrhynchos --> duck                  
3. Anser-cygnoides-orientalis --> goose         4. Bos-taurus --> cattle                        
5. Bubalus-bubalis --> water buffalo            6. Camelus-bactrianus --> Bactrian camel  
7. Camelus-dromedarius --> Dromedary camel      8. Canis-lupus --> dog                         
9. Capra-hircus --> goat                        10. Columba-livia --> dove                       
11. Equus-asinus --> donkey                     12. Equus-caballus --> horse  
13. Felis-catus --> cat                         14. Gallus-gallus --> chicken                   
15. Gorilla-gorilla --> gorilla                 16. Homo-sapiens --> human being                 
17. Macaca-mulatta --> macaque                  18. Mus-musculus --> mice    
19. Oryctolagus-cuniculus --> rabbit            20. Ovis-aries --> sheep                        
21. Rattus-norvegicus --> rat                   22. Sus scrofa --> pig  
 
6.``` -e   --enzymes```  
The name of the enzymes that selected in Double Digest should be connected with "+", and the default is "BamHI+NheI".  

7.``` -v   --vector```  
Users can upload carrier files by themselves. PVAX1 vector is used by default.  
                        
8.``` -us   --upstream_insert_sequence```
Users can insert base sequence after upstream enzyme digestion site according to the selected carrier. "GCCGCCACCGCCACCATG" is used by default.  
                        
9.``` -ds   --downstream_insert_sequence```  
Users can insert base sequence before downstream enzyme digestion site according to the selected carrier. "TAG" is used by default.  
  
10.``` -w   --way```  
If you choose to upload multiple sequence alignment file, please enter 'msa'. If you want to directly input the fasta file and use 'blastp' in AutoDVD2 for multiple sequence alignment, please enter 'blast'.  
  
11.``` -m   --msa_file```  
Please enter the absolute path of multiple sequence alignment file. We recommend users to use MEGA to process multiple sequence alignment file. Note that the ID of the destination sequence should be named 'Query'.If you need to input multiple files, the absolute paths of the file should be connected with "+".  

# Updates
The software will be updated on github from time to time. If you have any problems, you can send an email to 1602371546@qq.com.

# ATD2.0
AutoDVD2

## Author
Xinwen Zheng

## Email
1602371546@qq.com

# Introduction
AutoDVD2, which is also abbreviated as ATD2.0, is a Linux-based command-line-driven tool for automating the construction of DNA vaccines. AutoDVD2 can quickly construct DNA vaccine plasmids that can express conserved linear B cell epitopes. AutoDVD2 follows the basic construction framework of AutoDVD, but has improved the amino acid conservation analysis, antigen epitope screening and codon optimization. Both AutoDVD conservative analysis and codon optimization need to be connected to the Internet to capture the results of web page analysis. AutoDVD2 uses a new set of algorithms without networking, which greatly improves the stability and speed of the software.It integrates automated methods into an easy-to-use tool to assist researchers without programming skills.

The main functions of the software are as follows:
1. Linear B cell epitope screening
2. Epitope Conservation Analysis
3. Codon Optimization
4. Insertion sequence design suitable for double digestion vector construction
5. DNA Vaccine Plasmid Construction

# Usage
After installing dependent packages such as blast, you need to enter the following parameters:
(Items 1-5 are mandatory parameters, and items 6-7 are optional parameters)
1. -p [PDB number of pathogen protein]  
If you need to input multiple pathogen proteins, or select multiple chains for the same pathogen protein, the PDB number should be connected with "+" and capitalized, for example, "5X5C+5K6G".
2. -c [Specific chain of pathogen protein]  
If you need to input multiple pathogen proteins, or select multiple chains for the same pathogen protein, the selected chains should be connected with "+" and capitalized, such as "A+A". Note that the writing sequence of the chain and PDB number should correspond.
3. -if [absolute path of pathogen protein fasta file]  
Before uploading, please check whether the fasta file only contains sequences of one specific chain. If not, you need to manually delete redundant amino acid sequences. If you need to input multiple files, the absolute path of the file should be connected with "+" and capitalized, for example: "/ ATD2.0/demo/5X5C.fasta+./ATD2.0/demo/5K6G.fasta”。 Note that the writing sequence of the document and PDB number should correspond.
4. -ip [Absolute path of pathogen protein pdb file]  
If you need to input multiple files, the absolute path of the file should be connected with "+" and capitalized, for example: "/ ATD2.0/demo/5X5C.pdb+./ATD2.0/demo/5K6G.pdb”。 Note that the writing sequence of the document and PDB number should correspond.
5. -ho [Latin scientific name of host]  
The "./ATD2.0/codon" folder stores the codon usage tables of 22 common species. Users only need to enter the Latin scientific name of the host to be vaccinated. Note that the scientific name space should be replaced by "-", such as "Sus-scrofa". If there is no appropiate codon usage tables in the codon folder, users can download them from http://www.kazusa.or.jp/codon/, name it "<Latin scientific name (space is replaced by '-')>. txt", and put it in the codon folder. 
6. -e [Name of enzyme selected by double digestion]  
The name of the enzyme is connected with "+", and the default is "BamHI+NheI"
7. -v [absolute path of carrier file]  
Users can upload carrier files by themselves. PVAX1 carrier is used by default

# Vaccine construction instance  
(Enter the following command to find the relevant output file in the $HOME/$USER directory)  
python3 ./ATD2.0/sum.py -c A+A -p 5X5C+5K6G -if ./ATD2.0/demo/5X5C.fasta+./ATD2.0/demo/5K6G.fasta -ip ./ATD2.0/demo/5X5C.pdb+./ATD2.0/demo/5K6G.pdb -ho Sus-scrofa

## Updates
The software will be updated on github from time to time. If you have any problems, you can send an email to 1602371546@qq.com.

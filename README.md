# SNPMAP
SNPMAP uses NCBI database files to perform SNP annotation for proteins of interest 

Installation instructions to enable execution of python script snpmap.py from command line in Unix-based operating systems

In the Unix terminal, cd to the directory containing the snpmap script and modify it's permissions by entering:

$ chmod +x snpmap.py

Add /path/to/snpmap to PATH (e.g. for Bash login shell by entering):

$ export PATH="$PATH:/path/to/snpmap"

Please try the example sequence.gbx.xml and snp_result.xml files distributed with this script: cd into the ‘test’ directory and enter the following command then follow the on-screen instructions:

$ snpmap.py


INFORMATION ABOUT SNPMAP OUTPUT

The script will create a new folder corresponding to the protein of interest. Within this folder, the script will generate the following files:

- genomic.fasta: a fasta file with the wild-type and mutagenized genomic DNA locus

- cds.fasta: a fasta file with the wild-type and mutagenized coding sequence

- protein: a fasta file with the wild-type protein sequence

- snpmap.log: a tab-delimited text file containing information about the refSNP ID, the mutation type, the flanking genomic sequence, the position of the mutation in the CDS (#CDS) and protein (#Codon), the base variation at each SNP position (DNA) and the corresponding amino acid variation at each position (AA). Note that the order of the amino acids in the AA column is random; please use the snplist.txt file to identify wild-type/reference amino acids.

- snplist.txt: a tab-delimited text file containing columns of differently formatted input for SNP prediction servers (e.g. Provean, VEP, PolyPhen2)

- snpmap_output.txt: a tab-delimited text file containing the coding sequence and protein sequence annotated with the SNP information

- snpmap_idx.txt: a text file containing residue numbers corresponding to different sequences from the alignment files (optional)


OBTAINING SNP INFO FROM NCBI

NCBI:
(Resources -> Variation -> dbSNP)

http://www.ncbi.nlm.nih.gov/snp
 Search for <gene>
Filters:
- Organism: Homo sapiens
- Variation class: snp
- Annotation: protein
- Function Class: missense

Send to: File -> Format: XML -> Sort by: Chromosome Base Position


OBTAINING GENE INFO FROM NCBI

NCBI:
(Resources -> DNA and RNA -> Reference Sequence)

http://www.ncbi.nlm.nih.gov/refseq/
 Search for <gene>

Select chromosome result 
Filters:
- Organism: Homo sapiens
- genomic DNA/RNA
 Display GenBank (full)

If too size too large, click 'Display features'
Under 'Change region shown' click selected region and put in start and end positions of gene
N.B. 
Do not make changes to 'Customise view'
Send: Complete Record -> File -> INSDSeq XML -> Create File


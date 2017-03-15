#!/usr/bin/env python

print 'Single Nucleotide Polymorphism (SNP) mapper Version 2'
# Version 2 recognises the new *.gbx.xml INSDSeq file format

# Import required python modules and functions
import os
import math
import sys
from xml.etree import ElementTree as ET
from string import maketrans

# Input data files (mandatory)
gene_filename = ''
while len(gene_filename) == 0:
	try:
		gene_filename = raw_input('Input gene INSDSeq (*.gbx.xml): ')
		if gene_filename.find('.gbx.xml',-8) == -1:
			gene_filename="%s.gbx.xml" % (gene_filename) 
		# Open DNA sequence in extensible markup language (XML) with NCBI INSDSeq specification 
		INSDSeq = ET.ElementTree(file=gene_filename)
		INSDSeq = INSDSeq.getiterator('GBSeq')
	except KeyboardInterrupt:
		print '\nINTERRUPTED!!'		
		sys.exit(0)
	except:
		print "Error: Failed to open filename '%s'" % (gene_filename) 
		gene_filename = ''   
snp_filename = ''
while len(snp_filename) == 0:
	try:
		snp_filename = raw_input('Input SNP result (*.xml): ')
		if snp_filename.find('.xml',-4) == -1:
			snp_filename="%s.xml" % (snp_filename) 
		# Open SNP search result stored in extensible markup language (XML).
		Rs = ET.ElementTree(file=snp_filename)
		Rs = Rs.getiterator('Rs')
	except KeyboardInterrupt:
		print '\nINTERRUPTED!!'	
		sys.exit(0)
	except:
		print "Error: Failed to open filename '%s'" % (snp_filename) 
		snp_filename = ''  


# Create function to find the reverse complement of DNA/RNA sequences
def revcmpl(seq):
	if seq.find('u') != -1:
		seq = seq.replace('u','t')
	elif seq.find('U') != -1:
		seq = seq.replace('U','T')
	from string import maketrans
	seq = seq[::-1] 
	table=maketrans('atcgwsmkrybvdhnxATCGWSMKRYBVDHNX','tagcwskmyrvbhdnxTAGCWSKMYRVBHDNX')
	seq = seq.translate(table)
	return seq

def transeq(seq):
	# Convert RNA to DNA sequence
	seq = seq.lower().replace('\n', '').replace(' ', '')
	if seq.find('u') != -1:
		seq = seq.replace('u','t')

	# Make standard codon table
	bases = ['t','c','a','g']
	codons = [a+b+c for a in bases for b in bases for c in bases]
	amino_acids ='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
	codon_table = dict(zip(codons, amino_acids))

	# Translate DNA sequence into peptide sequence
	peptide = ''
	for i in xrange(0, len(seq), 3):
		aa = codon_table.get(seq[i: i+3], 'X')
		peptide += aa
	return peptide

def snplist():
    # SNPMAP output formatting for Effect Prediction Servers

    # Import required python modules
    import os
    from sys import platform

    # Read protein sequence
    protein_file = open('protein.fasta','r')
    protein_file.next() # skip header line
    sequence = protein_file.next()

    # Load snpmap results
    snpmap_log = open('snpmap.log','r')
    snpmap_log.next() # skip header line
    varpos = []
    count = 0
    while True:
        tmp = snpmap_log.next()
        if tmp[0:9] == 'Completed': 
            break
        else:
            tmp = tmp.split('\t')
            if tmp[4] != 'NA':
                varpos.append(tmp)
                count += 1

    # Loop through variable positions
    if platform == "linux" or platform == "linux2":
        # linux
        NCBI_RefSeq_ID = os.getcwd().rpartition('/')[2]
    elif platform == "darwin":
        NCBI_RefSeq_ID = os.getcwd().rpartition('/')[2]
        # OS X
    elif platform == "win32":
        # Windows...
        NCBI_RefSeq_ID = os.getcwd().rpartition('\\')[2]
    snplist = []
    for i in range(count):
        rsID = varpos[i][2]
        pos = varpos[i][4]
        ref = sequence[int(pos)-1]
        var = varpos[i][6].replace(",","")
        var = var.replace("\n","")
        if var.find(ref) == -1:
            print "Warning: Reference residue not found at position %s" % pos
        var = var.replace(ref,'')
        for j in range(len(var)):
            snplist.append([rsID,NCBI_RefSeq_ID,pos,ref,var[j]])

    # Prepare output
    output=""
    for i in range(len(snplist)):
        output += snplist[i][0]+"\t"
        output += snplist[i][1]+" "
        output += snplist[i][2]+" "
        output += snplist[i][3]+" "
        output += snplist[i][4]+"\t"
        output += snplist[i][1]+":p."
        output += snplist[i][3]
        output += snplist[i][2]
        output += snplist[i][4]+"\n"
    output += '\nhttp://provean.jcvi.org/protein_batch_submit.php?species=human'
    output += '\nhttp://genetics.bwh.harvard.edu/pph2/bgi.shtml'
    output += '\nhttp://www.ensembl.org/Tools/VEP\n'
    snplist_file = open('snplist.txt','w')
    snplist_file.write(output)
    snplist_file.close()

# Extract genomic sequence
seq = INSDSeq[0].find('GBSeq_sequence').text.upper()

# Extract coordinates delineating the consensus coding sequence
tmp = ET.tostring(INSDSeq[0].find('GBSeq_feature-table'))
N = tmp.count('<GBFeature_key>CDS</GBFeature_key>')
cds_data = {}
protein_seq = {}
seq_direction = []
for i in range(N):
	tmp = tmp.partition('<GBFeature_key>CDS</GBFeature_key>')[2]
	if tmp.partition('join(')[0].find('complement') == -1:
		seq_direction.append('+')
	else:
		seq_direction.append('-')
	tmp = tmp.partition('join(')[2]
	coord = tmp.partition(')')[0]
	coord = coord.partition(',')
	start = []
	end = []
	l = len(coord[0])
	while l != 0:
		start.append(coord[0].partition('..')[0])
		end.append(coord[0].partition('..')[2])
		coord = coord[2].partition(',')
		l = len(coord[0])
	start = [int(j) for j in start]
	start = [j-1 for j in start]
	end = [int(j) for j in end]
	if seq_direction[i] == '-':
		tmp_start = end
		tmp_start = [len(seq)-j for j in tmp_start]
		tmp_start = tmp_start [::-1] 
		tmp_end = start
		tmp_end = [len(seq)-j for j in tmp_end]
		tmp_end = tmp_end [::-1] 
		start = tmp_start
		end = tmp_end
	cds_data.update({'CDS ' + str(i+1):(start,end)})
	tmp = tmp.partition('<GBQualifier_name>protein_id</GBQualifier_name>')[2]
	tmp = tmp.partition('<GBQualifier_value>')[2]
	protein_seq.update({'Protein ' + str(i+1):tmp.partition('<')[0]})
	tmp = tmp.partition('<GBQualifier_name>translation</GBQualifier_name>')[2]
	tmp = tmp.partition('<GBQualifier_value>')[2]
	protein_seq['Protein ' + str(i+1)] = (protein_seq['Protein ' + str(i+1)],tmp.partition('<')[0])
	print str(i+1) + '\t' + str(len(protein_seq['Protein ' + str(i+1)][1])) + ' aa' + '\t' + protein_seq['Protein ' + str(i+1)][0]
pidx = ''
while len(str(pidx)) == 0:
	try:
		pidx = input('Input index of protein for annotation: ')
		start = cds_data['CDS ' + str(pidx)][0]
		end = cds_data['CDS ' + str(pidx)][1]
	except KeyboardInterrupt:
		print '\nINTERRUPTED!!'		
		sys.exit(0)
	except:
		print "Error: Index value is out of range"
		pidx = ''   
start = cds_data['CDS ' + str(pidx)][0]
end = cds_data['CDS ' + str(pidx)][1]

# Reverse-complment genomic and cDNA mutant coding sequences if applicable
if seq_direction[pidx-1] == '-':
	seq = revcmpl(seq)

# Fetch basic SNP information
rsId = []
status = []
source = []
dbSnpBuild = []
position_seq = []
direction = []
mutation = []
impact = []
for i in range(len(Rs)):
    rsId.append('rs'+str(Rs[i].attrib['rsId']))
    status.append(Rs[i].attrib['snpType'])
    source.append(Rs[0].find('PrimarySequence').attrib['source'])
    dbSnpBuild.append(Rs[0].find('PrimarySequence').attrib['dbSnpBuild'])
    mutation.append(Rs[i].find('Sequence')[1].text)
    if Rs[i].find('Phenotype') is not None:
        if Rs[i].find('Phenotype').find('ClinicalSignificance').text:
            impact.append(Rs[i].find('Phenotype').find('ClinicalSignificance').text)
        else:
            impact.append('')
    else:
        impact.append('')

# Scan each direction using 3' flanking sequence of SNP
for i in range(len(Rs)):
    # Try to find 3' flanking sequence on forward direction
    position_seq.append(seq.find(Rs[i].find('Sequence')[2].text))
    direction.append('+')  
    if position_seq[i] == -1:  
        # Try to find 3' flanking sequence on reverse direction
        position_seq[i] = seq.find(revcmpl(Rs[i].find('Sequence')[0].text)) 
        direction[i] = '-'
    if position_seq[i] != -1:
        position_seq[i] = position_seq[i] - 1

# In cases where there is not an exact match to the 3' flanking 
# sequence, scan each direction using 5' flanking sequence of SNP
for i in range(len(Rs)):
    if position_seq[i] == -1:
        # Determine length of 5' flanking sequence on forward direction
        l = len(Rs[i].find('Sequence')[0].text)
        # Try to find 5' flanking sequence on forward direction
        position_seq[i] = seq.find(Rs[i].find('Sequence')[0].text)
        direction[i] = '+'
        if position_seq[i] == -1:  
            # Determine length of 5' flanking sequence on reverse direction
            l = len(Rs[i].find('Sequence')[2].text)
            # Try to find 5' flanking sequence on reverse direction
            position_seq[i] = seq.find(revcmpl(Rs[i].find('Sequence')[2].text))
            direction[i] = '-'
        if position_seq[i] != -1:
            position_seq[i] = position_seq[i] + l

# Determine single-letter IUPAC abbreviation for variation at each SNP position
log = sorted(zip(position_seq,direction,mutation,rsId,impact))
SNP = []
err_count = 0
for i in range(len(Rs)):
    SNP.append(mutation[i])
    SNP[i] = SNP[i].replace('/','')
    SNP[i] = SNP[i].upper()
    if direction[i] == '-':
        SNP[i] = revcmpl(SNP[i])
IUPAC = {'AT':'W','CG':'S','AC':'M','GT':'K','AG':'R','CT':'Y','CGT':'B','AGT':'D','ACT':'H','ACG':'V','ACGT':'N'}
for i in range(len(Rs)):
	try:
		SNP[i] = IUPAC [SNP[i]] 
	except:
		SNP[i] = 'X' 
		err_count = err_count + 1	

# Create mutated genomic sequence
mutant_seq = ''
for i in range(len(Rs)):
    if log[i][0] != -1:
        if len(mutant_seq) == 0:
            mutant_seq = seq[0:log[i][0]] + sorted(zip(position_seq,SNP))[i][1]
        else:
            # Include the following if statement in case of duplicate entries 
            if log[i][0] != log[i-1][0]:
                mutant_seq = mutant_seq + seq[log[i-1][0]+1:log[i][0]] + sorted(zip(position_seq,SNP))[i][1]          
            if i == len(Rs)-1:
                mutant_seq = mutant_seq + seq[log[i][0]+1:len(seq)]
if len(mutant_seq) != len(seq):
    err_count = err_count +1

# Create mutant coding sequence
cds = ''
mutant_cds = ''
position_cds = []
for i in range(len(start)):
    cds = cds + seq[start[i]:end[i]]
    mutant_cds = mutant_cds + mutant_seq[start[i]:end[i]]
    position_cds.append(range(start[i],end[i]))
position_cds = [item for sublist in position_cds for item in sublist]
rsId_cds = ['' for i in range(len(mutant_cds))]
impact_cds = ['' for i in range(len(mutant_cds))]

# Map SNP positions to mutant coding sequence and annotate
table=maketrans('ACGTUWSMKRYBDHVNX','00000111111111110')
snpmap = mutant_cds.translate(table)
snpmap = list(snpmap)
snpmap = [int(i) for i in snpmap]
snpidx_cds = []
snpidx_pep = []
for i in range(len(Rs)): 
	try: 
		rsId_cds[position_cds.index(position_seq[i])] = rsId[i]
		impact_cds[position_cds.index(position_seq[i])] = impact[i]
		snpidx_cds.append(1+position_cds.index(position_seq[i]))
		snpidx_pep.append(int(math.ceil(snpidx_cds[i]/3.0)))
	except:
		snpidx_cds.append('NA')
		snpidx_pep.append('NA')
IUPAC_inv = dict(zip(IUPAC.values(),IUPAC.keys()))
AA=[]
snptype=[]
for i in range(len(Rs)):
	if snpidx_cds[i] == 'NA':
		AA.append('')
	else:
		for j in range(len(IUPAC_inv[SNP[i]])):
			tmp = list(cds)
			tmp[snpidx_cds[i]-1] = IUPAC_inv[SNP[i]][j]
			tmp=''.join(tmp)
			if j == 0:
				AA.append(transeq(tmp)[snpidx_pep[i]-1])
			elif AA[i].find(transeq(tmp)[snpidx_pep[i]-1]) == -1:
				AA[i] += ',%s' % (transeq(tmp)[snpidx_pep[i]-1])
	if len(set(AA[i].replace(',',''))) == 1:
		snptype.append('synonymous')
	elif len(set(AA[i].replace(',',''))) > 1:
		snptype.append('nonsynonymous')
	else:
		snptype.append('        ')
DNA_cds = ['' for i in range(len(mutant_cds))]
AA_cds = ['' for i in range(len(mutant_cds))]
snptype_cds = ['' for i in range(len(mutant_cds))]
status_cds = ['' for i in range(len(mutant_cds))]
for i in range(len(Rs)): 
	try:
		DNA_cds[position_cds.index(position_seq[i])] = SNP[i]
		AA_cds[position_cds.index(position_seq[i])] = AA[i]
		snptype_cds[position_cds.index(position_seq[i])] = snptype[i]
		status_cds[position_cds.index(position_seq[i])] = status[i]
	except:
		pass
AA_ref = list(protein_seq['Protein ' + str(pidx)][1])
protein_size = len(protein_seq['Protein ' + str(pidx)][1])
[AA_ref.insert(i+1,'..') for i in range(protein_size-1)[::-1]]
AA_ref=''.join(AA_ref)
snpmap = zip([i+1 for i in range(len(mutant_cds))],list(cds),list(AA_ref),snpmap,rsId_cds,status_cds,DNA_cds,AA_cds,snptype_cds,impact_cds)
snpmap = str(snpmap)
snpmap = snpmap[2:-2].replace('), (','\n')
snpmap = snpmap.replace(', ','\t')
snpmap = snpmap.replace("'","")

# Reproduce flanking sequence in original orientation to cross-check with dbSNP
seq_segment = []
for i in range(len(Rs)):
    if position_seq[i] != -1:
        if direction[i] == '+':
            seq_segment.append(seq[position_seq[i]-25:position_seq[i]] + ' ' + mutation[i] + ' ' + seq[position_seq[i]+1:position_seq[i]+26])
        elif direction[i] == '-':
            seq_segment.append(revcmpl(seq[position_seq[i]+1:position_seq[i]+26]) + ' ' + mutation[i] + ' ' + revcmpl(seq[position_seq[i]-25:position_seq[i]]))
    else:
        seq_segment.append('')

# Recreate log file to include remaining SNP information
log = sorted(zip(position_seq,direction,rsId,snpidx_cds,snpidx_pep,SNP,AA,snptype,seq_segment,impact))

# Input protein alignment file in FASTA format to map SNP positions to multiple related sequences
aln_filename = ''
try:
	# Load alignment sequences
	aln_filename = raw_input('Input protein multiple fasta alignment (optional): ')
	aln_file=open(aln_filename,'r').read()
	aln_seqID=[]
	aln_seq=[]
	tmp=aln_file.partition('>')
	while len(tmp[2]) > 0:
		aln_seqID.append(tmp[2].partition('\n')[0])
		tmp = tmp[2].partition('\n')[2].partition('>')
		aln_seq.append(tmp[0].replace('\n',''))

	# Check that the selected protein exactly matches one of those in the alignment.
	match = []
	for j in range(len(aln_seqID)):
		if aln_seq[j].replace('-','') == protein_seq['Protein ' + str(pidx)][1]:
			match.append('yes')
		else:
			match.append('no')
	match=match.index('yes')

	# Obtain SNP positions in selected protein sequence and remove duplicate entries
	aln_synonymous = [log[i][4] for i in range(len(Rs)) if log[i][4] != 'NA' and log[i][7] == 'synonymous']
	aln_nonsynonymous = [log[i][4] for i in range(len(Rs)) if log[i][4] != 'NA' and log[i][7] == 'nonsynonymous']
	aln_synonymous=sorted(list(set(aln_synonymous)))
	aln_nonsynonymous=sorted(list(set(aln_nonsynonymous)))

	# Correct SNP positions in alignment with reference to the selected protein sequence
	gap_count = 0
	tmp_synonymous = []
	for k in range(len(aln_seq[match])):
		if aln_seq[match][k] == '-':
			gap_count += 1
		else:
			if [aln_synonymous.index(k-gap_count) if k-gap_count in aln_synonymous else -1][0] > -1:
				tmp_synonymous.append(aln_synonymous[aln_synonymous.index(k-gap_count)] + gap_count)
	aln_synonymous = tmp_synonymous
	gap_count = 0
	tmp_nonsynonymous = []
	for k in range(len(aln_seq[match])):
		if aln_seq[match][k] == '-':
			gap_count += 1
		else:
			if [aln_nonsynonymous.index(k-gap_count) if k-gap_count in aln_nonsynonymous else -1][0] > -1:
				tmp_nonsynonymous.append(aln_nonsynonymous[aln_nonsynonymous.index(k-gap_count)] + gap_count)
	aln_nonsynonymous = tmp_nonsynonymous

	# Find SNP positions for all aligned sequences 
	snpmap_idx=''
	for j in range(len(aln_seqID)):
		synonymous_idx = ''
		for i in range(len(aln_synonymous)):
			if aln_seq[j][aln_synonymous[i]] != '-':
				synonymous_idx = synonymous_idx + str(len(aln_seq[j][0:aln_synonymous[i]].replace('-',''))) + '+'
		nonsynonymous_idx = ''
		for i in range(len(aln_nonsynonymous)):
			if aln_seq[j][aln_nonsynonymous[i]] != '-':
				nonsynonymous_idx = nonsynonymous_idx + str(len(aln_seq[j][0:aln_nonsynonymous[i]].replace('-',''))) + '+'
		synonymous_idx = synonymous_idx[0:-1]
		nonsynonymous_idx = nonsynonymous_idx[0:-1]
		snpmap_idx = snpmap_idx + '>' + aln_seqID[j] + '_synonymous' + '\n' + synonymous_idx + '\n'
		snpmap_idx = snpmap_idx + '>' + aln_seqID[j] + '_nonsynonymous' + '\n' + nonsynonymous_idx + '\n'
	snpmap_idx = snpmap_idx[0:-1]
except KeyboardInterrupt:
	print '\nINTERRUPTED!!'		
	sys.exit(0)
except:
	snpmap_idx='No alignment was loaded for mapping SNPs to alternative related sequences'
	print snpmap_idx

# Reformat mutation log for export 
if err_count == 0:
    summary = '\nCompleted successfully!'
elif err_count == 1:
    summary = '\nCompleted with %0.f error' % err_count
else:
    summary = '\nCompleted with %0.f errors' % err_count
iupac = '\nDNA IUPAC ambiguity codes: Y = C/T; R = A/G; W = A/T; S = G/C; K = T/G; M = C/A; D = A/G/T; V = A/C/G; H = A/C/T; B = C/G/T'
log = str(log)
log = log[2:-2]
log = log.replace("'","")
log = log.replace("), (",'\n')
log = log.replace(", ",'\t')
log = 'Loc\tStrand\tRefSNP ID\t#CDS\t#Codon\tDNA\tAA\tMutation Type\tFlanking sequence (original strand)\n' + log + summary + iupac
if seq_direction[pidx-1] == '-':
	log = log + '\nNote: The reverse complement of the input sequence was used in order to match the ORF orientation'

# Save log file and mutant DNA output
cwd = os.getcwd()
try:
	output_name = raw_input('Output folder [%s]: ' % (protein_seq['Protein ' + str(pidx)][0])) 
except KeyboardInterrupt:
	print '\nINTERRUPTED!!'		
	sys.exit(0)
if len(output_name) == 0:
	output_name = protein_seq['Protein ' + str(pidx)][0]
try: 
    os.chdir(output_name)
except:
    os.mkdir(output_name)
    os.chdir(output_name)
mutant_genomic_file = open('genomic.fasta','w')
mutant_genomic_file.write('>genomic\n' + seq + '\n\n>mutant_genomic\n' + mutant_seq)
mutant_genomic_file.close() 
mutant_cds_file = open('cds.fasta','w')
mutant_cds_file.write('>cds\n'+cds+'\n\n>mutant_cds\n'+mutant_cds)
mutant_cds_file.close()
protein_file = open('protein.fasta','w')
protein_file.write('>' + protein_seq['Protein ' + str(pidx)][0] + '\n' + protein_seq['Protein ' + str(pidx)][1])
protein_file.close() 
log_file = open('snpmap.log','w')
log_file.write(log)
log_file.close()
print log
snpmap_output_file = open('snpmap_output.txt','w')
snpmap_output_file.write('#CDS\tCDS\tPeptide\tSNP\tRefSNP ID\tStatus      \tDNA\tAA\tSNP type\tImpact\n'+snpmap)
snpmap_output_file.close()
snpmap_idx_file = open('snpmap_idx.txt','w')
snpmap_idx_file.write(snpmap_idx)
snpmap_idx_file.close()
snplist()
os.chdir(cwd)


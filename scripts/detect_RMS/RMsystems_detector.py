import os.path
import argparse as args

parser = args.ArgumentParser(description='Restriction-Modification systems detector (lsb@sanger.ac.uk, July 2017)', usage = '%(prog)s [options]')
parser.add_argument("-i", "--input_gff", help="GFF file.", required=True)
parser.add_argument("-p", "--pfam_dom", help="List of PFAM domains to search for.", required=True)
parser.add_argument("-m", "--pfam_hmm", help="Result of hmmscan on a Pfam database.", required=False)
parser.add_argument("-d", "--distance", help="Maximum distance (in number of genes) to cluster R and M (default = 10).", default=10, required=False)
parser.add_argument("-o", "--outfile", help="Prefix for output files.", required=True)

args = parser.parse_args()

## Functions ##

def check_input_data(args):
	gffParse = args.input_gff
	pfamDom = args.pfam_dom
	pfamHMM = args.pfam_hmm
	if not os.path.exists(gffParse):
		raise Exception("GFF file not found.\n")
	else:
		print("\n## GFF file: "+gffParse)
	if not os.path.exists(pfamDom):
		raise Exception("PFAM domain list not found.")
	else:
		print("## PFAM domain list: "+pfamDom)
	if pfamHMM is not None:
		if not os.path.exists(pfamHMM):
			raise Exception("PFAM hmmscan output not found.")
		else:
			print("## PFAM hmmscan output file: "+pfamHMM)
	else:
		print("## PFAM hmmscan output file: "+str(pfamHMM)+"\n")

def parse_GFF(gff):
	gffParse = {}
	gffFile = {}
	geneOrder = []
	sequence = {}
	tmpname = ''
	tmpseq = [] 
	with open(gff, 'rU') as gfile:
		for line in gfile:
			if not line.startswith('#') and 'CDS' in line:
				if '\t' in line:
					linesplit = line.rstrip().split('\t')
					annot = linesplit[8]
					annotsplit = annot.split(';')
					locus_tag = [x for x in annotsplit if 'locus_tag=' in x][0]
					locus_tag = locus_tag.replace('locus_tag=', '').replace('"', '')
					geneOrder.append(locus_tag)
					gffFile[locus_tag] = line
					pfam = [x for x in annotsplit if 'Pfam:' in x]
					if len(pfam)>0:
						pfam = pfam[0].split(':')
						gffParse[locus_tag] = pfam[len(pfam)-1].split('.')[0]
				else:
					gline = line.strip()
					if gline.startswith('>'):
						if len(tmpseq)>0:
							sequence[tmpname] = ''.join(tmpseq)
						tmpseq = []
						tmpname = gline.replace('>', '')
					else:
						if tmpname != '':
							tmpseq.append(gline)
		sequence[tmpname] = ''.join(tmpseq)
	return gffParse, geneOrder, gffFile, sequence

def parse_Pfam_list(pfam):
	# First column MUST be the Pfam domain
	pfamList = {}
	with open(pfam, 'rU') as plist:
		for line in plist:
			if not line.startswith('Pfam'):
				linesplit = line.rstrip().split('\t')
				pfamList[linesplit[0].split('.')[0]] = linesplit[1:]
	return pfamList

def parse_Pfam_HMM(pfam):
	# MUST be created using --tblout (per-seq parseable table)
	pfamHMM = {}
	with open(pfam, 'rU') as pHMM:
		for line in pHMM:
			if not line.startswith('#'):
				linesplit = line.rstrip().split(' ')
				linesplit = [x for x in linesplit if x != '']
				target = linesplit[0]
				accession = linesplit[1].split('.')[0]
				gene = linesplit[2]
				if gene in pfamHMM:
					pfamHMM[gene].append(accession)
				else:
					pfamHMM[gene] = [accession] 
	return pfamHMM

## Main ##
check_input_data(args)

# Get Pfam list of target domains
pfamList = parse_Pfam_list(args.pfam_dom)
#print pfamList

# Parse GFF file and get annotated Pfam information
gffParse = parse_GFF(args.input_gff)

# Get results:
result = {}
RMresult = {}
if args.pfam_hmm is not None:
	# Parse hmmscan output file
	pfamHMM = parse_Pfam_HMM(args.pfam_hmm)
	#print pfamHMM
	for gene in pfamHMM:
		pfams = pfamHMM[gene]
		for p in pfams:
			if p in pfamList:
				if gene in result:
					result[gene].append(p)
					RMresult[gene].append(pfamList[p][2])
				else:
					result[gene] = [p]
					RMresult[gene] = [pfamList[p][2]]

# Add the info from annotated Pfams in the GFF
for gene in gffParse[0]:
	p = gffParse[0][gene]
	if p in pfamList:
		if gene in result:
			if p not in result[gene]:
				result[gene].append(p)
				RMresult[gene].append(pfamList[p][2])
		else:
			result[gene] = [p]
			RMresult[gene] = [pfamList[p][2]]

# Detect RMs (R/M <= '-d' genes apart)
geneOrder = gffParse[1]
endonucleases = []
methylases = []
for gene in RMresult:
	t = ','.join(sorted(list(set(RMresult[gene]))))
	if 'M' in t:
		methylases.append(gene)
	if 'R' in t:
		endonucleases.append(gene)

indices = sorted([geneOrder.index(x) for x in RMresult])

rmsystems = []
for i in indices:
	ini = i-args.distance
	end = i+args.distance
	match = [x for x in indices if x in range(ini,end+1)]
	# does it overlap with any already existing RMs?
	if len(match)>0:
		if len(rmsystems)==0:
			rmsystems = [match]
		else:
			for rm in rmsystems:
				check = []
				for x in match:
					if x in rm:
						check.append(True)
					else:
						check.append(False)
				if True not in check:
					todo = 'no_overlap'
				elif False not in check:
					todo = 'overlap'
				else:
					todo = 'partial_overlap'
			# Modify the list outside the loop
			if todo == 'no_overlap':
				rmsystems.append(match)
			elif todo == 'partial_overlap':
				[rmsystems[len(rmsystems)-1].append(x) for x in match if x not in rmsystems[len(rmsystems)-1]]

# Get only grouped genes
RMclean = {}
count = 0
for rm in rmsystems:
	if len(rm)>1:
		count += 1
		RMclean['RM-'+str(count)] = rm

# Report results:

# Detected PFAM domains
outfile = args.outfile+'_PFAM.txt'
with open(outfile, 'w') as out:
	for res in result:
		gene = res
		pfams = ','.join(result[gene])
		pfamsType = ','.join(sorted(list(set(RMresult[gene]))))
		out.write(gene+'\t'+pfamsType+'\t'+pfams+'\n')

# Summary file
outfile = args.outfile+'_summary.txt'
grouped_loci = []
for x in RMclean:
	for y in RMclean[x]:
		grouped_loci.append(y)
ungrouped_methylases = [x for x in methylases if geneOrder.index(x) not in grouped_loci]
ungrouped_endonucleases = [x for x in endonucleases if geneOrder.index(x) not in grouped_loci]
with open(outfile, 'w') as out:
	out.write('## Number of RM systems: '+str(len(RMclean))+'\n')
	out.write('## Number of total MTases: '+str(len(methylases))+'\n')
	out.write('## Number of MTases in RM: '+str(len(methylases)-len(ungrouped_methylases))+'\n')
	out.write('## Number of ungrouped MTases: '+str(len(ungrouped_methylases))+'\n')
	out.write('## Number of total REases: '+str(len(endonucleases))+'\n')
	out.write('## Number of REases in RM: '+str(len(endonucleases)-len(ungrouped_endonucleases))+'\n')
	out.write('## Number of ungrouped REases: '+str(len(ungrouped_endonucleases))+'\n\n')
	for rm in RMclean:
		labs = []
		for x in RMclean[rm]:
			labs.append(geneOrder[x])
		out.write(rm+'\t'+'|'.join(labs)+'\n')

# Gff file
outfile = args.outfile+'_RMs.gff'
gffFile = gffParse[2]
sequence = gffParse[3]
with open(outfile, 'w') as out:
	out.write('##gff-version 3\n')
	for seq in sequence:
		out.write('##sequence-region '+seq+' 1 '+str(len(sequence[seq]))+'\n')
	# RM systems
	for x in RMclean:
		for y in RMclean[x]:
			locus_tag = geneOrder[y]
			unit = ','.join(sorted(list(set(RMresult[locus_tag]))))
			if unit == 'M':
				col = 4 # blue
			elif unit == 'R':
				col = 2 # red
			else:
				col = 11 # brown
			toprint = gffFile[locus_tag].rstrip()+';note='+x+'('+unit+' unit)'+';colour='+str(col)+'\n'
			out.write(toprint)
	# Ungrouped MTases
	for x in ungrouped_methylases:
		col = 9
		toprint = gffFile[x].rstrip()+';note=ungrouped_methylase;colour='+str(col)+'\n'
		out.write(toprint)
	# Ungrouped REases
	for x in ungrouped_endonucleases:
		col = 16
		toprint = gffFile[x].rstrip()+';note=ungrouped_endonuclease;colour='+str(col)+'\n'
		out.write(toprint)

print("Finished.\n\n")

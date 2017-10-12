import os.path
import argparse as args

parser = args.ArgumentParser(description='Motif/Methylation enrichment by gene (lsb@sanger.ac.uk, August 2017)', usage = '%(prog)s [options]')
parser.add_argument("-i", "--input_gffs", help="Path to input gffs", required=True)
parser.add_argument("-r", "--roary_csv", help="gene_presence_absence.csv file from Roary", required=True)
parser.add_argument("-p", "--path_to_IPDs", help="Path to folder with IPD .txt files", required=True)
parser.add_argument("-b", "--methylated_bases", help="Methylated bases (number starting with 1, i.e. -b 3 for 6mA in CCACC. Multiple bases can be specified separated by commas.)", required=True)
parser.add_argument("-t", "--threshold", help="IPD value to use as threshold (default = 3)", required=False, default=3)

args = parser.parse_args()
pathtogffs = args.input_gffs
pathtoipds = args.path_to_IPDs
input_gffs = os.listdir(pathtogffs)
input_gffs = [x for x in input_gffs if not x.startswith('.')]
methbases = args.methylated_bases.split(',')
methbases = [int(x) for x in methbases]
threshold = float(args.threshold)

## Functions ##

def parse_GFF(gff):
	gffFile = {} 
	with open(gff, 'r') as gfile:
		for line in gfile:
			if not line.startswith('#'):
				if '\t' in line:
					linesplit = line.rstrip().split('\t')
					coords = linesplit[3]+'..'+linesplit[4]
					annot = linesplit[8]
					annotsplit = annot.split(';')
					locus_tag = [x for x in annotsplit if 'ID=' in x][0]
					locus_tag = locus_tag.replace('ID=', '').replace('"', '')
					gffFile[locus_tag] = coords
	return gffFile

def parse_Roary(tab):
	allclusters = []
	with open(tab, 'r') as csv:
		header = csv.readline().rstrip().replace('"', '').split(',')
		strains = header[14:]
		roaryDic = {}
		for line in csv:
			linesplit = line.rstrip().replace('"', '').split(',')
			cluster = linesplit[0]
			allclusters.append(cluster)
			straininfo = linesplit[14:]
			for s in range(len(straininfo)):
				if straininfo[s] != '':
					if strains[s] not in roaryDic:
						roaryDic[strains[s]] = {}
						roaryDic[strains[s]][cluster] = straininfo[s]
					else:
						roaryDic[strains[s]][cluster] = straininfo[s]
	return [roaryDic, allclusters]

def checkMotifs(f, gene_start, gene_end):
	mot_in_gene = []
	mot_methyl = []
	with open(f, 'r') as ipds:
		for line in ipds:
			linesplit = line.rstrip().split('\t')
			if linesplit[1] != 'IPDstrand':
				start = int(linesplit[2])
				end = int(linesplit[3])
				# Check if motif is inside the gene
				if start > gene_start and end < gene_end:
					mot_in_gene.append(start)
					# Check if motif has any IPD value > threshold
					ipdvalues = [float(x) for x in linesplit[8:]]
					above = False
					for base in methbases:
						if ipdvalues[base-1]>threshold:
							above = True
					if above is True:
						mot_methyl.append(start)
	# Account for potential redundancy
	mot_in_gene = list(set(mot_in_gene))
	mot_methyl = list(set(mot_methyl))
	return [len(mot_in_gene), len(mot_methyl)]

## Main ##

# Parse Roary table
roarytab = args.roary_csv #'gene_presence_absence.csv'
roary = parse_Roary(roarytab)
roarytab = roary[0]
allclusters = roary[1]

motif_files = os.listdir(pathtoipds)
motif_files = [x for x in motif_files if x.endswith('_IPDs.txt')]

byStrain = {}
for strain in roarytab:
	print strain+'...'
	strainname = strain.replace('Ngono_','').replace('Ngonorrhoeae_','').replace('_with_plasmids','')
	gffile = [x for x in input_gffs if strain in x]
	gffpath = pathtogffs+'/'+gffile[0]
	gff = parse_GFF(gffpath)
	byGene = {}
	for gene in gff:
		coords = gff[gene].split('..')
		gene_start = int(coords[0])
		gene_end = int(coords[1])
		filename = [x for x in motif_files if strainname in x][0]
		f = pathtoipds+'/'+filename
		check = checkMotifs(f, gene_start, gene_end)
		byGene[gene] = check # number of motifs, number of methylated motifs
	byStrain[strain] = byGene

## Build final table

strains = sorted(byStrain.keys())

outTab_mot = {}
outTab_meth = {} 
for cluster in allclusters:
	mot_row = []
	meth_row = []
	for st in range(len(strains)):
		strain = strains[st]
		if cluster in roarytab[strain]:
			genes = roarytab[strain][cluster]
			genes = genes.split('\t')
			# Account for paralogs, sum values
			mot = 0
			met = 0
			for g in genes:
				genedata = byStrain[strain][g]
				mot = mot+genedata[0]
				met = met+genedata[1]
			mot_row.append(str(mot))
			meth_row.append(str(met))
		else:
			mot_row.append('')
			meth_row.append('')
	outTab_mot[cluster] = mot_row
	outTab_meth[cluster] = meth_row

## Write output files
outfile = pathtoipds+'/'+'motif_enrichment.txt'
with open(outfile, 'w') as out:
	header = 'cluster'+','+','.join(strains)+'\n'
	out.write(header)
	for cluster in outTab_mot:
		outline = cluster+','+','.join(outTab_mot[cluster])+'\n'
		out.write(outline)

## Write output files
outfile = pathtoipds+'/'+'methylation_enrichment.txt'
with open(outfile, 'w') as out:
	header = 'cluster'+','+','.join(strains)+'\n'
	out.write(header)
	for cluster in outTab_meth:
		outline = cluster+','+','.join(outTab_meth[cluster])+'\n'
		out.write(outline)


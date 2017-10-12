import os.path
import argparse as args

parser = args.ArgumentParser(description='Converts InterProScan GFF output file into a GENE<TAB>GO_TERMS_LIST format (lsb@sanger.ac.uk, August 2017)', usage = '%(prog)s [options]')
parser.add_argument("-i", "--input", help="Input GFF file", required=True)
parser.add_argument("-o", "--output", help="Output file prefix", required=True)

args = parser.parse_args()

geneDic = {}
with open(args.input, 'r') as gff:
	for line in gff:
		if not line.startswith('#') and '\t' in line:
			linesplit = line.rstrip().split('\t')
			gene = linesplit[0]
			annot = linesplit[8].split(';')
			# Check if GO term
			catchGO = [x for x in annot if x.startswith('Ontology_term')]
			if len(catchGO)==1:
				catchGO = catchGO[0].replace('Ontology_term=', '').replace('"', '')
				catchGO = catchGO.split(',')
				if gene in geneDic:
					if len(geneDic[gene])>0:
						add = [geneDic[gene].append(x) for x in catchGO if x not in geneDic[gene]]
				else:
					geneDic[gene] = catchGO

# Write output file #
outfile = args.output+'.GO.txt'
with open(outfile, 'w') as out:
	for i in geneDic:
		out.write(i+'\t'+','.join(geneDic[i])+'\n')

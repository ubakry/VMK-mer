"""
Copyright 2020 | Ahmed Omar, Asmaa Ali, Mohamed AboelEla, and Usama Bakry
https://github.com/ubakry/VMK-mer

This file is used for running VMK-mer directly without installing it.

VMK-mer is a free tool under the terms of the GNU General Public License v3.0.
VMK-mer is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License v3.0 for 
more details. You should have received a copy of the GNU General Public License v3.0 along with VMK-mer. If
not, see <http://www.gnu.org/licenses/>.
"""

## VMK-mer Description
print('==================================================================================')
print('VMK-mer v1.0 | by A. Omar, A. Ali, M. AboelEla, and U. Bakry')
print('Check https://github.com/ubakry/VMK-mer for updates.')
print('==================================================================================')

# Importing libraries
import argparse
import pysam
import pandas as pd
import sys

args = None

# ----------------------------------------------------------------------
def get_args():
	"""
	get_args() function creates required and optional arguments and returns a copy of the argument list to be used within the script.
	"""
	parser = argparse.ArgumentParser(
		description="VMK-mer - Standalone tool that converts mutations in VCF file into\nk-mer sequences that are affected by these mutations.",
		epilog="This is where you might put example usage"
	)

	# required argument
	parser.add_argument('-f', action="store", required=True, help='Input fasta file (*.fasta)')
	parser.add_argument('-v', action="store", required=True, help='Input vcf file (*.vcf)')
	parser.add_argument('-k', action="store", required=True, type=int, help='Length of k-mer (e.g. 5)')

	# optional arguments
	parser.add_argument('-o', action="store", help='The output directory - default: current directory', default='.')
	parser.add_argument('--outfmt', action="store", help='The output file type [TSV or XML] - default: both', default='both')
	parser.add_argument('--outfile', action="store", help='The output file name without extention - default: vmkmer-results', default='vmkmer-results')
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	arguments = vars(parser.parse_args())
	return arguments
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
def write_in_tsv(mut,mut_dict):
	"""
	This function generate the output in tsv format.
	"""
	with open(args['o']+'/'+args['outfile']+'.tsv','a') as fd:
		fd.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(mut_dict['chr'], 
					mut_dict['pos'], mut_dict['id'], mut_dict['ref'], mut_dict['alt'], mut_dict['refk'], mut_dict['mutk']))
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
def write_in_xml(mut,mut_dict):
	"""
	This function generate the output in xml format.
	"""
	with open(args['o']+'/'+args['outfile']+'.xml','a') as fd:
		fd.write('<{}>\n\t<Chr>{}</Chr>\n\t<Pos>{}</Pos>\n\t<Mutation-ID>{}</Mutation-ID>\n\t<Ref-Allele>{}</Ref-Allele>\n\t<Mut-Allele>{}</Mut-Allele>\n\t<Ref-Kmers>{}</Ref-Kmers>\n\t<Mut-Kmers>{}</Mut-Kmers>\n</{}>\n'.format(mut,mut_dict['chr'],mut_dict['pos'], mut_dict['id'], mut_dict['ref'], mut_dict['alt'], mut_dict['refk'], mut_dict['mutk'],mut))
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
def add(mut,mut_dict):
	"""
	This function write the output in the required format which the user has selected using the --outfmt argument.
	"""
	if args['outfmt'].upper() == 'TSV':
		write_in_tsv(mut,mut_dict)
	elif args['outfmt'].upper() == 'XML':
		write_in_xml(mut,mut_dict)
	elif args['outfmt'] == 'both':
		write_in_tsv(mut,mut_dict)
		write_in_xml(mut,mut_dict)
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
def snp(record, genome, k):
	
	# Extracting the sequence surrounding the position of the mutation. 
	seq = genome.fetch(record.chrom, record.pos-k, record.pos+k-1)
	
	# This statement handles the issue of empty "alt" field in the current record of vcf file.
	if record.alts == None:
		ref_kmers = []
		mutant_kmers = []
		alt='N'
		mut_seq = seq[:k-1]+alt+seq[k:]

		for i in range(1,k+1):
			ref_kmers.append(seq[i-1:k+i-1])
			mutant_kmers.append(mut_seq[i-1:k+i-1])
	else:
		# loop on each alt value in the alt field of the current record.
		for alt in record.alts:
			ref_kmers = []
			mutant_kmers = []
			mut_seq = seq[:k-1]+alt+seq[k:]

			# loop on the length of kmer entered by the user to create the suitable number of kmers.
			# For (SNP) The number of kmers = the length of kmer.
			for i in range(1,k+1):

				ref_kmers.append(seq[i-1:k+i-1])
				mutant_kmers.append(mut_seq[i-1:k+i-1])

	snp_dict={'chr': record.chrom, 'pos': record.pos, 'id': record.id, 'ref':record.ref , 'alt':alt, 'refk':ref_kmers, 'mutk':mutant_kmers}
	add("SNP",snp_dict)
# ----------------------------------------------------------------------	

# ----------------------------------------------------------------------	
def insertion(record, genome, k):
	
	seq = genome.fetch(record.chrom, record.pos-k, record.pos+k-1)
	
	# loop on each alt value in the alt field of the current record.
	for alt in record.alts:
		ref_kmers = []
		mutant_kmers = []
		mut_seq = seq[:k]+alt[1:]+seq[k:]

		# loop on the length of kmer entered by the user to create the suitable number of kmers.
		# For (insertion) The number of Reference kmers = the length of kmer 
		for i in range(1,k+1):
			ref_kmers.append(seq[i-1:k+i-1])

		# and the number of Mutant kmers = length of kmer+length of the inserted nucleotides.
		for i in range(1, k+len(alt)):
			mutant_kmers.append(mut_seq[i-1:k+i-1])
                  
		ins_dict={'chr': record.chrom, 'pos': record.pos, 'id': record.id, 'ref':record.ref , 'alt':alt, 'refk':ref_kmers, 'mutk':mutant_kmers}
		add("Insertion",ins_dict)
# ----------------------------------------------------------------------	

# ----------------------------------------------------------------------
def deletion(record, genome, k):
	
	# loop on each alt value in the alt field of the current record.
	for alt in record.alts:
		ref_kmers = []
		mutant_kmers = []
		seq = genome.fetch(record.chrom, record.pos-k, record.pos+k-2+len(record.ref))
		mut_seq = seq[:k+1]+seq[len(record.ref)+k:]
		
		# loop on the length of kmer entered by the user to create the suitable number of kmers.
		# For (insertion) The number of Reference kmers = length of kmer+length of the inserted nucleotides.
		for i in range(1, k+len(record.ref)):
			ref_kmers.append(seq[i-1:k+i-1])
			
		# and the number of Mutant kmers = the length of kmer.
		for i in range(1,k+1):
			mutant_kmers.append(mut_seq[i-1:k+i-1])

		del_dict={'chr': record.chrom, 'pos': record.pos, 'id': record.id, 'ref':record.ref , 'alt':alt, 'refk':ref_kmers, 'mutk':mutant_kmers}
		add("Deletion",del_dict)
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
def progress(iterations, total):
    bar_len = 50
    fill_len = int(round(bar_len * iterations / float(total)))
    percent = round(iterations / float(total) * 100.0, 1)
    bar = "#" * fill_len + "-" * (bar_len - fill_len)
    sys.stdout.write("%s [%s] %s%s   %s\r" %
                     ("Extracting K-mers", bar, percent, "%", "Complete"))
    sys.stdout.flush()


def recordcount(v, vcf):
    f = open(args["v"], 'r')
    header = str(vcf.header).splitlines()
    total = sum(1 for line in f)-len(header)+1
    return total
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
class FileFormatError(IOError):
	'''raise this error when input file format is not acceptable'''

# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
def main():
    
	# Exception handling for input files format.
	if args['f'].endswith(".fa") or args['f'].endswith(".fasta"):
		genome = pysam.FastaFile(args['f'])  # open fasta file
		print('[	  OK       ] Reading Fasta file is done.')
	else:
		raise FileFormatError("\n[	  ERROR    ] Input File is not in Fasta format.")

	if args['v'].endswith(".vcf"):
		vcf = pysam.VariantFile(args['v'])  # open vcf file
		print('[	  OK       ] Reading vcf file is done.')
	else:
		raise FileFormatError("\n[	  ERROR    ] Input File is not in VCF format.")
		
	k = args['k']  # Length of kmer

	# Handling different output formats.
	if args['outfmt'].upper() == 'TSV':
		with open(args['o']+'/'+args['outfile']+'.tsv','w') as fd:
			fd.write('## VMK-mer version: v1.0\n## Output file: {}\n## Reference fasta file: {}\n## VCF file: {}\n'.format(args['o']+'/'+args['outfile'],args['f'],args['v']))
		head_dict={'chr': 'Chr', 'pos': 'Pos', 'id': 'Mutation-ID', 'ref':'Ref-Allele' , 'alt':'Mut-Allele', 'refk':'Ref-Kmers', 'mutk':'Mut-Kmers'}
		write_in_tsv("Head",head_dict)
        
	elif args['outfmt'].upper() == 'XML':
		with open(args['o']+'/'+args['outfile']+'.xml','w') as fd:
			fd.write('## VMK-mer version: v1.0\n## Output file: {}\n## Reference fasta file: {}\n## VCF file: {}\n'.format(args['o']+'/'+args['outfile'],args['f'],args['v']))
	
	elif args['outfmt'] == 'both':
		with open(args['o']+'/'+args['outfile']+'.tsv','w') as fd:
			fd.write('## VMK-mer version: v1.0\n## Output file: {}\n## Reference fasta file: {}\n## VCF file: {}\n'.format(args['o']+'/'+args['outfile'],args['f'],args['v']))
		with open(args['o']+'/'+args['outfile']+'.xml','w') as fd:
			fd.write('## VMK-mer version: v1.0\n## Output file: {}\n## Reference fasta file: {}\n## VCF file: {}\n'.format(args['o']+'/'+args['outfile'],args['f'],args['v']))
	


	print('[	PROCESS    ] Extracting mutant kmers, please wait...')

	total = recordcount("v", vcf)
	iterations = 0
	for record in vcf:
        
		# Handling the mutation type included in the info tag of some VCF file as "TSA".
		if 'TSA' in record.info.keys():
			
			mutation_type = str(record.info['TSA'])
			if mutation_type == "SNV":
				snp(record, genome, k)
				#pass
                
			elif mutation_type == "insertion":
				insertion(record, genome, k)
				#pass
                
			elif mutation_type == "deletion":
				deletion(record, genome, k)
				#pass

        	# Handling the mutation type included in the info tag of some VCF file as "VT".
		elif 'VT' in record.info.keys():

			mutation_type = str(record.info['VT'][0])
			if mutation_type == "SNP":
				snp(record, genome, k)
				#pass
                
			elif mutation_type == "INDEL":
				if len(record.alts[0]) > len(record.ref):
					insertion(record, genome, k)
					#pass
                    
				elif len(record.alts[0]) < len(record.ref):
					deletion(record, genome, k)
					#pass
		iterations += 1
		progress (iterations, total)
	print('\n[	  OK       ] All kmers have been extracted successfully.')
#------------------------------------------------------------------------------

if __name__ == '__main__':
	args = get_args()
	main()

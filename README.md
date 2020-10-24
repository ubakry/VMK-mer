# VMK-mer : Standalone tool that converts mutations in VCF file into the k-mer sequences that were affected by these mutations.

## Introduction
A new approach for genome-wide association studies (GWAS) relies of sequencing data instead of microarrays. In this approach, some tools build their association studies based on k-mers frequency that change between healthy and diseased individuals. Instead of counting the k-mers throughout the whole genome, we are building a tools that generates k-mers, of any size, only around the sites of mutations, whether SNPs or indels. This will reduce the required computational power needed for such GWAS studies. Moreover, the output can be used for other disease-networks studies.

## Aim
Create a tool that converts VCF file mutations to list of k-mers. For each mutation in the VCF file, it should find the k-mers that got affected by this specific mutation. It should be able to handle SNPs, insertions, deletions, and multiple mutations at the same locus. Finally, the tool should be easy-to-use with multiple options for output file format, output directory, and k-mer length.

## Manual
This tools was built with Python and requires the installation of the following two Python libraries: [pandas](https://pandas.pydata.org/) and [pysam](https://pysam.readthedocs.io/en/latest/installation.html).

### Main arguments
The following are the required arguments to run VMK-mer:

- `-f <fasta_file>`: input fasta file. It should be the same reference used for variant calling.
- `-v <vcf_file>`: input VCF file. VMK-mer should handle all VCF formats till v.4.3.
- `-k <int>`: the k-mer size used of the mutated k-mers.

### Optional arguments
The following are extra arguments that can be used with VMK-mer:

- `-o <path>`: Output file path (directory). default is the current working directory.
- `--outfmt <TSV|XML>`: specifies the output file format (`TSV` or `XML`). The default mode will keep both files.
- `--outfile <str>`: The basename of the output file; _i.e._ without extention - default: `vmkmer-results`.
- `-h|--help`:  show the help message (manual) of the tool and exit.
- `--version`:   show program's version number and exit.


### Using VMK-mer
```bash
python vmkmer.py -f <fasta_file> -v <vcf_file> -k <k-mer size (5)> [-o <output_file_path>] [--outfile <base output file name>] [--outfmt <output file format (TSV or XML)>]

Main arguments
  -f F        Input fasta file (*.fasta or *.fa)
  -v V        Input vcf file (*.vcf)
  -k K        Length of k-mer (e.g. 5)

Optional arguments:
  -h, --help  show the help message and exit
  -o O        The output directory - default: current directory
  --outfmt    output file format (TSV or XML). default value would produce both files.
  --outfile OUTFILE  The basename of the output file; _i.e._ without extention - default: "vmkmer-results"
  --version   show program's version number and exit

```


### VMK-mer WORKFLOW 
<p align="center">
  <img src="https://github.com/ubakry/VMK-mer/blob/main/misc/vmkmer-workflow.jpg"  width="90%" height="90%">
</p>


### VMK-mer TEST CASES
This tool was tested on the files of Human Chromosome 01, 21, 22 and Y.
You can download these files using the following commands:

```bash
For chr01:
wget ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
wget ftp://ftp.ensembl.org/pub/release-96/variation/vcf/homo_sapiens/homo_sapiens-chr1.vcf.gz

For chr21:
wget ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
wget ftp://ftp.ensembl.org/pub/release-96/variation/vcf/homo_sapiens/homo_sapiens-chr21.vcf.gz

For chr22:
wget ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
wget ftp://ftp.ensembl.org/pub/release-96/variation/vcf/homo_sapiens/homo_sapiens-chr22.vcf.gz

For chrY:
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
```

### Feedback and bug reports

Your comments, bug reports, and suggestions are very important to improve VMK-mer for you. You can leave your comments and bug reports at our GitHub issues.


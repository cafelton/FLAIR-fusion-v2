# FLAIR-fusion-v2
New version of long-read fusion isoform detection

THIS IS CURRENTLY STILL IN DEVELOPMENT, please post any issues to the wiki

Basic requirements: python 3+, pysam
FLAIR installed using conda and conda environment activated - https://flair.readthedocs.io/en/latest/requirements.html

This is long-read fusion and fusion isoform detection meant for Oxford Nanopore or PacBio sequencing data. It is tailored to remove artifactual chimeras due to alignment and library prep errors.
One advantage of this tool is the ability to detect alternative splicing in gene fusions.

This tool works best with well-annotated species such as human and mouse, as it primarily uses the annotated genes to detect gene fusions. Future releases will include more detction of fusions in unannotated regions.
It has been tested with multiple human and mouse gencode releases (https://www.gencodegenes.org/human/release_38.html), so please try to match the formatting of annotations from other sources to gencode.

preprocessing: This tool works with highly error-prone data, but that will slow it down so we reccommend removing very short reads (<300bp) before running. DO NOT filter for read quality, as that can remove real chimeras. This is particularly important for data known to have many short artifacts, such as single-cell data.

FLAIR-fusion preprocessing: FLAIR-fusion needs to generate two reference files - a transcriptome intron location reference and a transcriptome homology reference. This takes a while, so make sure you only use the -q option to generate these files once per genome/transcriptome reference and afterwards provide the file locations wil -e and -p


```bash
usage: python[3+] fusionfindingpipeline.py -r reads.[fq/fa] -t transcriptome.fa -g genome.fa -a annotation.gtf [-m OR -s readsAlignedToTranscriptome.bam] [-q OR -e path.tsv -p path.tsv] [other options] -i

FLAIR-fusion 2.0 parse options

options:
  -h, --help            show this help message and exit
  
  -g G, --genome G      path to genome
  
  -d D, --scratchFolder D
                        path to scratch folder for writing large files in preprocessing
                        
  -t T, --transcriptome T
                        path to transcriptome (.fa)
                        
  -r R, --reads R       .fa or fq file
  
  -s S, --alignedReads S
                        .bam file that has a matching index
                        
  -e E, --intronCoords E
                        path to intron to genome coords file (.tsv)
                        
  -p P, --paralogReference P
                        path to intron to genome coords file (.tsv)
                        
  -b B, --buffer B      length of buffer for calling alignments as too close on genomic scale (bp)
  
  -l L, --readSupport L
                        number of reads required to call fusion
                        
  -a A, --anno A        path to anno.gtf
  
  -i, --callIsoforms    whether to detect fusion isoforms
  
  -q, --runPreprocessing
                        whether to run preprocessing steps (intron to genome and homology reference making)
                        
  -m, --alignTranscriptome
                        whether to align reads to transcriptome, if this is not selected, need to give .bam file with -s option
  -o O, --output O      output file name base, if not specified, will be derived from reads file name. This will prefix all output files.
```

OUTPUTS

There are currently many output files, future releases will have a more trimmed version of output files. The most important ones are as follows:

filePrefix.transcriptomeAligned-fusionReadCounts.tsv: This is the fusions identified and their read support.

filePrefix.transcriptomeAligned-rejectedChimerasAfterParaRemoved.tsv: This is all chimeras that FLAIR-fusion threw out in filtering, good for troubleshooting

filePrefix-fusionOnly.genomeAligned-flair.collapse.isoforms.bed: These are the final fusion isoforms detected, each fusion will represent two lines of the .bed file, one line for the alignment to each locus. These lines will have the same name. This is also the file to look at for final predictions of fusion breakpoints.

filePrefix.syntheticAligned-flair.collapse.combined.isoform.read.map.txt: These are the final isoforms with all reads supporting each isoform. The total reads supporting all isoforms of the fusion will likely be less than the number in the ReadCounts.tsv file, as some reads are lost in the isoform identification process. If you want more precision on which reads support which isoforms (and likely more read support for each isoform), feel free to run FLAIR-quantify using the .flair.collapse.isoforms.fa file and the filePrefix-fusionsOnly.[fa/fq] file.

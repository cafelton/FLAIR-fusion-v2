import sys, os, argparse
from collections import Counter
from datetime import date
import subprocess
from statistics import median,stdev
import time


path = os.path.dirname(os.path.realpath(__file__))

parser = argparse.ArgumentParser(description='FLAIR-fusion 2.0 parse options',
                                 usage='python[3+] fusionfindingpipeline.py -r reads.[fq/fa] -t transcriptome.fa -g genome.fa -a annotation.gtf [-m OR -s readsAlignedToTranscriptome.bam] [-q OR -e path.tsv -p path.tsv] [other options] -i')
parser.add_argument('-o', '--output', action='store', dest='o',
                    help='output file name base, if not specified, will be derived from reads file name. This will prefix all output files.')
# parser.add_argument('-f', '--flair', action='store', dest='f', default="",
#                     help='flair path')
parser.add_argument('-g', '--genome', action='store', dest='g',
                    default="",
                    help='path to genome')
parser.add_argument('-d', '--scratchFolder', action='store', dest='d',
                    default="",
                    help='path to scratch folder for writing large files in preprocessing')
# parser.add_argument('-k', '--remapSize', action='store', dest='k', default=0, type=int, help='size of area to remap - only remaps if this is specified')
parser.add_argument('-t', '--transcriptome', action='store', dest='t',
                    default="",
                    help='path to transcriptome (.fa)')
parser.add_argument('-r', '--reads', action='store', dest='r', default="", help='.fa or fq file')
parser.add_argument('-s', '--alignedReads', action='store', dest='s', default="",
                    help='.bam file that has a matching index')
parser.add_argument('-e', '--intronCoords', action='store', dest='e',
                    default="",
                    help='path to intron to genome coords file (.tsv)')
parser.add_argument('-p', '--paralogReference', action='store', dest='p',
                    default="",
                    help='path to intron to genome coords file (.tsv)')
parser.add_argument('-b', '--buffer', action='store', dest='b', default='50000',
                    help='length of buffer for calling alignments as too close on genomic scale (bp)')
parser.add_argument('-l', '--readSupport', action='store', dest='l', default='3',
                    help='number of reads required to call fusion')
parser.add_argument('-a', '--anno', action='store', dest='a',
                    default="",
                    help='path to anno.gtf')
parser.add_argument('-i', '--callIsoforms', action='store_true', dest='i',
                    help='whether to detect fusion isoforms')
parser.add_argument('-q', '--runPreprocessing', action='store_true', dest='q',
                    help='whether to run preprocessing steps (intron to genome and homology reference making)')
parser.add_argument('-m', '--alignTranscriptome', action='store_true', dest='m',
                    help='whether to align reads to transcriptome, if this is not selected, need to give .bam file')
# /private/groups/brookslab/reference_annotations/
args = parser.parse_args()
overallstart = time.time()
prefix = '.'.join(args.r.split('.')[:-1])
if args.o: prefix = args.o

if len(args.d) > 0 and args.d[-1] != '/': args.d += '/'

if not os.path.isfile(args.r):
    raise Exception('reads file does not exist')
elif args.r.split('.')[-1] not in ['fa', 'fasta', 'fq', 'fastq']:
    raise Exception('reads must be in fasta or fastq format')

if args.q: ##run preprocessing
    if len(args.t) == 0 or len(args.a) == 0:
        raise Exception('please provide transcriptome.fa and annotation.gtf')
    elif args.t.split('.')[-1] not in ['fa', 'fasta'] or args.a.split('.')[-1] != 'gtf':
        raise Exception('transcriptome must be .fa or .fasta and annotation must be in gtf format')
    elif not os.path.isfile(args.t):
        raise Exception('transcriptome file does not exist')
    elif not os.path.isfile(args.a):
        raise Exception('annotation file does not exist')
    else:
        start = time.time()
        subprocess.call([sys.executable, path + '/simulateReadsFromIdentity.py', args.t, args.d])
        print('simulating errors on transcriptome done')
        process = subprocess.Popen('minimap2 -a -N4 ' + args.t + ' ' + args.d + 'sim-avg-100x-' + args.t.split('/')[-1] +
            ' | samtools view -bS - | samtools sort - -o '+ args.d +'sim-avg-100x-' + '.'.join(args.t.split('/')[-1].split('.')[:-1]) + '.transcriptomeAligned.bam; samtools index ' + args.d + 'sim-avg-100x-' + '.'.join(args.t.split('/')[-1].split('.')[:-1]) + '.transcriptomeAligned.bam',
            stdout=subprocess.PIPE, shell=True)
        print(process.communicate()[0].strip())
        subprocess.call([sys.executable, path + 'clusterAlignedParalogs-transcriptome-pysam.py', args.a, args.d + 'sim-avg-100x-' + '.'.join(args.t.split('/')[-1].split('.')[:-1]) + '.transcriptomeAligned.bam'])
        print('creating homology graph done')
        subprocess.call([sys.executable, path + 'transcriptToGenomeCoords.py', args.a])
        print('creating introns to genome done')
        args.e = 'transcriptome_introns_to_genome_coords_' + '.'.join(args.a.split('/')[-1].split('.')[:-1]) + '.tsv'
        args.p = args.t.split('/')[-1].split('.')[0] +  "TranscriptomeGeneToNeighbors-filteredReadLen.tsv"
        print('total preprocessing time: ', time.time()-start)

print(prefix)
if args.m: #align reads to transcriptome
    if not os.path.isfile(args.t):
        raise Exception('transcriptome file does not exist')
    elif args.t.split('.')[-1] not in ['fa', 'fasta']:
        raise Exception('transcriptome must be .fa or .fasta and annotation must be in gtf format')
    else:
        start = time.time()
        process = subprocess.Popen('minimap2 -a -N4 ' + args.t + ' ' + args.r +
                ' | samtools view -bS - | samtools sort - -o ' + prefix + '.transcriptomeAligned.bam; samtools index ' + prefix + '.transcriptomeAligned.bam',
                stdout=subprocess.PIPE, shell=True)
        print(process.communicate()[0].strip())
        args.s = prefix + '.transcriptomeAligned.bam'
        print('alignment to transcriptome done')
        print('transcriptome alignment time: ', time.time() - start)


if args.s == '': args.s = prefix + '.transcriptomeAligned.bam'

if not os.path.isfile(args.e):
    raise Exception('intron to genome file does not exist')
if not os.path.isfile(args.p):
    raise Exception('homology file does not exist')
if not os.path.isfile(args.s):
    raise Exception('aligned .bam file does not exist')
if not os.path.isfile(args.s + '.bai'):
    raise Exception('bam file index does not exist, index your file please')

start = time.time()
subprocess.call([sys.executable, path + '/removeParalogsGetChim-07-18-23.py', '-r', args.r, '-s', args.s, '-e', args.e, '-p', args.p, '-b', args.b, '-l', args.l, '-a', args.a, '-o', prefix])
print('base fusion finding done')
print('total fusion finding time: ', time.time()-start)

if args.i: #want fusion isoforms and further filtering
    if not os.path.isfile(args.g):
        raise Exception('genome file does not exist')
    start = time.time()
    subprocess.call([sys.executable, path + '/make_synthetic_fusion_reference-06-27-2023.py', '-g', args.g, '-a', args.a, '-r', prefix + 'chimericBreakpoints.tsv', '-o', prefix])
    print('synthetic fusion genome and annotation creation done')

    process = subprocess.Popen('minimap2 -ax splice --secondary=no -G 1000k ' + prefix + '-syntheticFusionGenome.fa ' + prefix +
                               '-fusionOnly.' + args.r.split('.')[-1] + ' | samtools view -bS - | samtools sort - -o ' + prefix + '-fusionOnly.syntheticAligned.bam;' +
                               ' bamToBed -bed12 -i ' + prefix + '-fusionOnly.syntheticAligned.bam > ' + prefix + '-fusionOnly.syntheticAligned.bed',
                               stdout=subprocess.PIPE, shell=True)
    print(process.communicate()[0].strip())
    print('realignment of fusion reads to synthetic genome done')

    subprocess.call(
        ['flair', 'correct', '-q', prefix + '-fusionOnly.syntheticAligned.bed',
         '-g', prefix + '-syntheticFusionGenome.fa',
         '-f', prefix + '-syntheticReferenceAnno.gtf',
         '--output', prefix + '-fusionOnly.syntheticAligned-flair'])

    subprocess.call(
        ['flair', 'collapse', '-q', prefix + '-fusionOnly.syntheticAligned-flair_all_corrected.bed',
         '-r', prefix + '-fusionOnly.' + args.r.split('.')[-1],
         '-g', prefix + '-syntheticFusionGenome.fa',
         '--gtf', prefix + '-syntheticReferenceAnno.gtf',
         '--annotation_reliant', 'generate', '--generate_map', '--check_splice',
         '--output', prefix + '-fusionOnly.syntheticAligned-flair.collapse'])
    print('flair collapse done')
    # print(path + '/convertSyntheticToGenomeBed.py', prefix + '-fusionOnly.syntheticAligned-flair.collapse.isoforms.bed', prefix + '-fusionOnly ' + args.r.split('.')[-1])
    subprocess.call([sys.executable, path + '/convertSyntheticToGenomeBed.py', prefix + '-fusionOnly.syntheticAligned-flair.collapse.isoforms.bed', prefix + '-fusionOnly.' + args.r.split('.')[-1]])
    print('synthetic converted to genome positions')
    process = subprocess.Popen(
        'minimap2 -ax splice -N 4 ' + args.g + ' ' + prefix +
        '-fusionOnly-isoSupport.' + args.r.split('.')[-1] + ' | samtools view -bS - | samtools sort - -o ' + prefix + '-fusionOnly-isoSupport.genomeAligned.bam;' +
        ' samtools index ' + prefix + '-fusionOnly-isoSupport.genomeAligned.bam',
        stdout=subprocess.PIPE, shell=True)
    print(process.communicate()[0].strip())
    #minimap2 -ax splice --secondary=no -G 1000k GRCm39.primary_assembly.genome.fa vollmers-mouse-r10-r2c2-all-fusionOnly-isoSupport.fasta | samtools view -bS - | samtools sort - -o vollmers-mouse-r10-r2c2-all-fusionOnly-isoSupport.bam
    print('fusion isoform finding done')
    print('total isoform finding time: ', time.time() - start)

print('total overall time: ', time.time() - overallstart)


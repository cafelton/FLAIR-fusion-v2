import sys, os, argparse
from statistics import median,stdev
from datetime import date

parser = argparse.ArgumentParser(description='FLAIR-fusion 2.0 parse options', usage='python3 realignToFilteredGenome2.py  ')
parser.add_argument('-r', '--chimBp', action='store', dest='r', default="", help='.fa or fq file')
parser.add_argument('-g', '--genome', action='store', dest='g', default="", help='path to genome')
parser.add_argument('-a', '--anno', action='store', dest='a', default="", help='path to anno.gtf')
parser.add_argument('-o', '--output', action='store', dest='o',
                    help='output file name base, if not specified, will be derived from reads file name. This will prefix all output files.')
args = parser.parse_args()

prefix = '.'.join(args.r.split('.')[:-2])
if args.o: prefix = args.o

def revComp(seq):
    newseq = ''
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N':'N'}
    for char in seq[::-1]:
        newseq += comp[char]
    return newseq

#####DONELoad in transcriptome and genome breakpoints and process them into one list of breakpoint locations
####DONEGo through transcript reference and load in gene start/end locations
####DONE    figure out if any predicted breakpoints are in the same intron node and collapse them
####Cut genes at all predicted breakpoints, label with gene names and cut locations, make synthetic fasta
####make synthetic annotation for these sequences

###check if any reads map to wrong orientation of fusion loci

allBP = {}
fgeneslist = set()
###['fusionName', 'geneName', 'orderInFusion', 'geneChr', 'breakpointCoord', 'outerEdgeCoord', 'readSupport']
for line in open(args.r):#'31-01-2023DRR059313-transcriptomeChimericBreakpoints-correctDir.tsv'):
    if line[:6] != 'fusion':
        line = line.split('\t')
        fusion,gene,isfive,thisChr,bp,outer = line[0], line[1], line[2],line[3], int(line[4]), int(line[5])
        fusion = tuple(fusion.split('--'))
        if fusion not in allBP: allBP[fusion] = {"5'gene":[], "3'gene":[]}
        # if gene not in allBP[fusion]: allBP[fusion][gene] = []
        allBP[fusion][isfive].append((gene, thisChr, bp, outer))
        for g in fusion:
            fgeneslist.add(g)

##load in genomic sequence
genome = {}
last = None
for line in open(args.g):#"/private/groups/brookslab/reference_sequence/GRCh38.primary_assembly.genome.fa"):
    if line[0] == '>':
        last = line.lstrip(">").split(" ")[0]
        genome[last] = []
    else: genome[last].append(line.rstrip('\n'))
for c in genome:
    genome[c] = "".join(genome[c])





fgenes = {}
for g in fgeneslist:
    fgenes[g] = {'bounds':(0,0)}#, 'splicesites':[]}
####To make synthetic transcriptome:
####DONE Get transcript/exon annotation for fusion genes
####    all annotation is recorded in plain left-right direction
####Filter this annotation to 5'/3' ends based on each breakpoint
####    Convert annotation values to be 0-based depending on start of gene (5' end) or breakpoint location (3' end)
####        Make sure to flip - strand values accordingly
####When making synthetic references, simulatneously make gtf annotation file - make sure to convert 3' side values based on
transcripts = {}

for line in open(args.a):#'/private/groups/brookslab/reference_annotations/gencode.v38.annotation.gtf'):
    if line[0] != '#':
        line = line.split('\t')
        if line[2] == 'gene' or line[2] == 'exon':
            genename = line[8].split('; gene_name "')[1].split('"')[0]
            genename += '*' + line[8].split('gene_id "')[1].split('"')[0]
            if genename in fgenes:
                if line[2] == 'gene':
                    ###learned that can't assume that transcript appears in anno only once - two diff ENSG can have same hugo name
                    if fgenes[genename]['bounds'] == (0,0):
                        fgenes[genename]['bounds'] = (line[0], int(line[3])-1, int(line[4]), line[6])
                    else:
                        fgenes[genename]['bounds'] = (line[0], min([int(line[3]) - 1, fgenes[genename]['bounds'][1]]), max([int(line[4]), fgenes[genename]['bounds'][2]]), line[6])
                elif line[2] == 'exon':
                    # fgenes[genename]['splicesites'].append(int(line[3]))
                    # fgenes[genename]['splicesites'].append(int(line[4]))
                    if genename not in transcripts: transcripts[genename] = {}
                    tname = line[8].split('; transcript_name "')[1].split('"')[0]
                    if tname not in transcripts[genename]: transcripts[genename][tname] = []
                    if line[6] == '+': transcripts[genename][tname].append((int(line[3])-1, int(line[4])))
                    else: transcripts[genename][tname].insert(0,(int(line[3])-1, int(line[4])))

out = open(prefix + '-syntheticFusionGenome.fa', 'w')#'syntheticFusionGenomeAttempt4.fa', 'w')
annoOut = open(prefix + '-syntheticReferenceAnno.gtf', 'w')#'syntheticReferenceAnnoAttempt1.gtf', 'w')
bpOut = open(prefix + '-syntheticBreakpointLoc.bed', 'w')#'syntheticFusionBreakpointLoc.bed', 'w')
for fusion in allBP:
    labels, sequence = [], []
    seqlen = 0
    isosByEnd = {"5'gene":{}, "3'gene":{}}
    startLoc = 0
    for end in ["5'gene", "3'gene"]:
        gene, thisChr = allBP[fusion][end][0][0], allBP[fusion][end][0][1]
        medianBp, medianOuter = median([x[2] for x in allBP[fusion][end]]), median([x[3] for x in allBP[fusion][end]])
        if medianBp < medianOuter:
            finalBp = min([x[2] for x in allBP[fusion][end]])
            finalOuter = fgenes[gene]['bounds'][2]
            sequence.append(genome[thisChr][finalBp:finalOuter])
        else:
            finalBp = max([x[2] for x in allBP[fusion][end]])
            finalOuter = fgenes[gene]['bounds'][1]
            sequence.append(genome[thisChr][finalOuter:finalBp])
        ###TEMP
        seqlen += abs(finalOuter-finalBp)

        # print(gene, 'fasta', finalOuter, finalBp, startLoc, startLoc + abs(finalOuter-finalBp))
        # print(gene, fgenes[gene]['bounds'][3], fgenes[gene]['bounds'][1], fgenes[gene]['bounds'][2], finalBp)
        labels.append('.'.join([str(x) for x in [gene, thisChr, finalBp, finalOuter]]))
        if fgenes[gene]['bounds'][3] == '-':
            sequence[-1] = revComp(sequence[-1])
        ###NEED TO ADD ALTERNATIVE ANNOTATION FOR ALTERNATIVE BREAKPOINTS, MAKE EXTRA TRANSCRIPT ANNOTATION
        for tname in transcripts[gene]:
            isosByEnd[end][tname] = []
            # if 'CCDC6' in tname:
            #     print(tname, end, fgenes[gene]['bounds'][3], medianBp, medianOuter, finalBp, finalOuter, transcripts[gene][tname])
            if fgenes[gene]['bounds'][3] == '+':
                for exon in transcripts[gene][tname]:
                    if end == "5'gene":
                        if exon[1] < finalBp:
                            isosByEnd[end][tname].append((exon[0]-fgenes[gene]['bounds'][1], exon[1]-fgenes[gene]['bounds'][1]))
                            startLoc = finalBp - fgenes[gene]['bounds'][1]
                    else:
                        if exon[0] > finalBp:
                            isosByEnd[end][tname].append(((exon[0]-finalBp)+startLoc, (exon[1]-finalBp)+startLoc))
            else:
                for exon in reversed(transcripts[gene][tname]):
                    if end == "5'gene":
                        if exon[0] > finalBp:
                            isosByEnd[end][tname].append((fgenes[gene]['bounds'][2]-exon[1], fgenes[gene]['bounds'][2]-exon[0]))
                            startLoc = fgenes[gene]['bounds'][2] - finalBp
                    else:
                        if exon[1] < finalBp:
                            isosByEnd[end][tname].append(((finalBp-exon[1])+startLoc, (finalBp-exon[0])+startLoc))
            # if len(isosByEnd[end][tname]) > 0:
            #     print(tname,fgenes[gene]['bounds'][3],end,isosByEnd[end][tname][-1])
    # print(isosByEnd)
    for end in ["5'gene", "3'gene"]:
        seen = []
        for iso in list(isosByEnd[end].keys()):
            if isosByEnd[end][iso] in seen:
                isosByEnd[end].pop(iso)
                # print(iso)
            else: seen.append(isosByEnd[end][iso])
    bpOut.write('\t'.join(['--'.join(labels), str(startLoc), str(startLoc), 'breakpoint']) + '\n')
    out.write('>' + '--'.join(labels) + '\n')
    out.write(''.join(sequence) + '\n')
    # annoOut.write('\t'.join(['--'.join(labels), 'SYNTHFUSION', 'gene', '1', str(len(''.join(sequence))+1), '.', '+', '.','gene_id "' + '--'.join(fusion) + '"']) + '\n')
    annoOut.write('\t'.join(['--'.join(labels), 'SYNTHFUSION', 'gene', '1', str(seqlen+1), '.', '+', '.','gene_id "' + '--'.join(fusion) + '"']) + '\n')

    for fiveIso in isosByEnd["5'gene"]:
        if len(isosByEnd["5'gene"][fiveIso]) > 0:
            for threeIso in isosByEnd["3'gene"]:
                if len(isosByEnd["3'gene"][threeIso]) > 0:
                    annoOut.write('\t'.join(['--'.join(labels), 'SYNTHFUSION', 'transcript', str(isosByEnd["5'gene"][fiveIso][0][0]+1), str(isosByEnd["3'gene"][threeIso][-1][-1]), '.', '+', '.','; '.join(['gene_id "' + '--'.join(fusion) + '"', 'transcript_id "' + '--'.join([fiveIso, threeIso]) + '"'])]) + '\n')
                    for exon in isosByEnd["5'gene"][fiveIso]:
                        annoOut.write('\t'.join(['--'.join(labels), 'SYNTHFUSION', 'exon', str(exon[0]+1),str(exon[1]), '.', '+', '.', '; '.join(['gene_id "' + '--'.join(fusion) + '"', 'transcript_id "' + '--'.join([fiveIso, threeIso]) + '"'])]) + '\n')
                    for exon in isosByEnd["3'gene"][threeIso]:
                        annoOut.write('\t'.join(['--'.join(labels), 'SYNTHFUSION', 'exon', str(exon[0]+1),str(exon[1]), '.', '+', '.','; '.join(['gene_id "' + '--'.join(fusion) + '"', 'transcript_id "' + '--'.join([fiveIso, threeIso]) + '"'])]) + '\n')



# out.close()
annoOut.close()
bpOut.close()


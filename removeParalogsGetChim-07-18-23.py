import re, time, sys, os, pysam, argparse
from collections import defaultdict
from statistics import median


def binarySearch(arr, t):
    if t <= arr[0]: return arr[0]
    if t >= arr[-1]: return arr[-1]
    i, j, mid = 0, len(arr) - 1, 0
    while i < j:
        mid = int((i + j) / 2)
        if arr[mid] == t:
            return arr[mid]
        elif t < arr[mid]:
            if mid > 0 and t > arr[mid - 1]:
                if abs(arr[mid] - t) < abs(arr[mid - 1] - t):
                    return arr[mid]
                else:
                    return arr[mid - 1]
            j = mid
        else:
            if mid < len(arr) - 1 and t < arr[mid + 1]:
                if abs(arr[mid] - t) < abs(arr[mid + 1] - t):
                    return arr[mid]
                else:
                    return arr[mid + 1]
            i = mid + 1


path = os.path.dirname(os.path.realpath(__file__))

parser = argparse.ArgumentParser(description='FLAIR-fusion 2.0 parse options',
                                 usage='python3 removeParalogsGetChim.py  [options]')
parser.add_argument('-r', '--reads', action='store', dest='r', default="", help='.fa or fq file')
parser.add_argument('-s', '--alignedReads', action='store', dest='s', default="",
                    help='.bam file that has a matching index')
parser.add_argument('-e', '--intronCoords', action='store', dest='e',
                    default="",
                    help='path to intron to genome coords file (.tsv)')
parser.add_argument('-p', '--paralogReference', action='store', dest='p',
                    default="",
                    help='path to intron to genome coords file (.tsv)')
parser.add_argument('-b', '--buffer', action='store', dest='b', default=50000,
                    help='length of buffer for calling alignments as too close on genomic scale')
parser.add_argument('-l', '--readSupport', action='store', dest='l', default=3,
                    help='number of reads required to call fusion')
parser.add_argument('-a', '--anno', action='store', dest='a',
                    default="",
                    help='path to anno.gtf')
parser.add_argument('-o', '--output', action='store', dest='o',
                    help='output file name base, if not specified, will be derived from reads file name. This will prefix all output files.')
args = parser.parse_args()
args.b = int(args.b)
args.l = int(args.l)
prefix = '.'.join(args.s.split('.')[:-1])
if args.o: prefix = args.o

paralogs = {}
for line in open(args.p):
    line = line.rstrip().split('\t')
    paralogs[line[0]] = set(line[1].split(','))
print('paralog reference processed')

genePos = {}
# cdsPos = {}
for line in open(args.a):
    if line[0] != '#':
        line = line.split('\t')
        if line[2] == 'gene':
            genename = line[8].split('"; gene_name "')[1].split('"')[0] + '*' + \
                       line[8].split('gene_id "')[1].split('"')[0]
            genePos[genename] = (line[0], int(line[3]), int(line[4]), line[6])
            # cdsPos[genename] = [int(line[3]), int(line[4])]
        # elif line[2] == 'start_codon':
        #     genename = line[8].split('"; gene_name "')[1].split('"')[0] + '*' + \
        #                line[8].split('gene_id "')[1].split('"')[0]
        #     if line[6] == '+': cdsPos[genename][0] = int(line[3])
        #     else: cdsPos[genename][1] = int(line[4])
        # elif line[2] == 'stop_codon':
        #     genename = line[8].split('"; gene_name "')[1].split('"')[0] + '*' + \
        #                line[8].split('gene_id "')[1].split('"')[0]
        #     if line[6] == '+': cdsPos[genename][1] = int(line[4])
        #     else: cdsPos[genename][0] = int(line[3])

print('gene pos reference loaded')

###start of gene annotated as 0.start+/-500.start, end annotated as tend.end.end+/-500
intronLocs, intronToGenome = {}, {}
for line in open(args.e):
    line = line.split('\t')
    if line[1] not in intronLocs:
        intronLocs[line[1]] = []
        intronToGenome[line[1]] = {}
    for intron in line[3].split(','):
        pos = [int(i) for i in intron.split('.')]
        intronLocs[line[1]].append(pos[0])
        intronToGenome[line[1]][pos[0]] = (pos[1], pos[2])
    intronLocs[line[1]] = sorted(intronLocs[line[1]])
print('intron to genome reference loaded')

flagbinary = {'0': 0, '2048': 0, '16': 1, '2064': 1, '256': 0, '272': 1}

# aligncount = {}
# alignlen = {}
timestart = time.time()
alignlen = defaultdict(lambda: defaultdict(list))
##NEW VERSION WITH PYSAM
##this now uses bam, not sam file
samfile = pysam.AlignmentFile(args.s,
                              "rb")  # "sim-avg-10x-gencode38-fusion-sim-test-06-12-2023.transcriptomeAligned.sorted.bam", "rb")
for s in samfile:
    if s.is_mapped:
        readname = s.query_name
        geneinfo = s.reference_name.split('|')
        genename, tname = geneinfo[5] + '*' + geneinfo[1], geneinfo[4]
        # if genename == 'Gm49410*ENSMUSG00000115979.2' or genename == 'Apol7a*ENSMUSG00000010601.15':
        tstart, tend = s.reference_start, s.reference_end
        # fastqstart, fastqend = s.query_alignment_start, s.query_alignment_end ###THIS DOESN'T WORK, IGNORES HARD-CLIPPED BASES
        cigar = s.cigartuples
        readlen = s.infer_read_length()
        if cigar[0][0] == 0:
            fastqstart = 0  # 'M':
        else:
            fastqstart = cigar[0][1]
        if cigar[-1][0] == 0:
            fastqend = readlen
        else:
            fastqend = readlen - cigar[-1][1]
        if s.is_reverse:
            fastqstart, fastqend = readlen - fastqstart, readlen - fastqend
        thisalignlen = s.get_cigar_stats()[0][0]
        alignlen[readname][genename].append([fastqstart, fastqend, thisalignlen, tname, tstart, tend])
print(time.time() - timestart)
print('alignment file processed')

###save start and end positions on fastq read
##while going through chim, get outside start and end positions for each gene for each read
##for each read find fastq dist between genes
##get median fastq dist for each read

##DONE get distance between chim genes, remove genes not a sufficent distance apart


chimToReads = {}
totChimReads = 0
for r in alignlen:
    if len(alignlen[r].keys()) > 1:
        totChimReads += 1
        chimname = frozenset(alignlen[r].keys())
        if chimname not in chimToReads: chimToReads[chimname] = set()
        chimToReads[chimname].add(r)
print('chimeras compressed')
print('total aligned reads, chimeric fraction', len(alignlen.keys()), totChimReads / len(alignlen.keys()))
print('total chimeras, chimeric reads', len(chimToReads.keys()), totChimReads)


# @profile
def removeParalogs(chimToReads, alignlen):
    paraRemovedChimToReads = {}
    tempParaSets = {}
    for chim in chimToReads:
        if len(chimToReads[chim]) > 0:  ##removed #more than one read support
            genes = list(chim)
            # print('genes', genes)
            paralogSets1 = []
            for g in genes:
                theseParas = {g, }
                shortG = g.split('*')[0]
                if shortG in paralogs:
                    for g2 in genes:
                        if g2.split('*')[0] in paralogs[shortG]: theseParas.add(g2)
                paralogSets1.append(theseParas)
            # print('paras', paralogSets1)
            # print([x.split('*')[0]])
            # print(paralogSets1)
            paraSets = []
            while len(paralogSets1) > 0:
                if len(paralogSets1) > 1:
                    indexToRemove = [0]
                    combSet = paralogSets1[0]
                    for i in range(1, len(paralogSets1)):
                        if len(paralogSets1[0] & paralogSets1[i]) > 0:
                            combSet = combSet | paralogSets1[i]
                            indexToRemove.append(i)
                    paraSets.append(combSet)
                    paralogSets1 = [x for i, x in enumerate(paralogSets1) if i not in indexToRemove]
                else:
                    paraSets.append(paralogSets1[0])
                    paralogSets1 = []
            if len(paraSets) > 1:
                frozenParaSets = frozenset([frozenset(x) for x in paraSets])
                if frozenParaSets not in tempParaSets:
                    tempParaSets[frozenParaSets] = chimToReads[chim]
                else:
                    tempParaSets[frozenParaSets] = tempParaSets[frozenParaSets] | chimToReads[chim]
    # del aligncount
    del chimToReads

    ###combine different chimeras with paralogs, very inefficient, not sure it can be faster though

    parasetsl = list(tempParaSets.keys())
    parasetslist = []
    for i in parasetsl:
        totgenes = 0
        for j in i:
            totgenes += len(j)
            # if 'ZNF765' in j: print('parasetslist', i)
        parasetslist.append((totgenes, i))
    parasetslist.sort(reverse=True)
    # for i in parasetslist: print(i)
    del parasetsl

    tempParaSets2 = {}
    while len(parasetslist) > 0:
        if len(parasetslist) > 1:
            indexToRemove = [0]
            combSetOuter = parasetslist[0][1]
            for i in range(2):  ##traverse twice to make sure we didn't miss anything the first time around
                for i in range(1, len(parasetslist)):
                    if i not in indexToRemove:
                        combSet = []
                        for paragenes1 in combSetOuter:
                            for paragenes2 in parasetslist[i][1]:
                                if len(paragenes1 & paragenes2) > 0:
                                    combSet.append(paragenes1 | paragenes2)
                        if len(combSet) == len(parasetslist[0][1]):
                            combSet2 = []
                            for paragenes1 in combSetOuter:
                                for paragenes2 in combSet:
                                    if len(paragenes1 & paragenes2) > 0:
                                        combSet2.append(paragenes1 | paragenes2)
                            combSetOuter = frozenset(combSet2)
                            indexToRemove.append(i)
            tempParaSets2[combSetOuter] = set()
            for i in indexToRemove:
                tempParaSets2[combSetOuter] = tempParaSets2[combSetOuter] | tempParaSets[parasetslist[i][1]]
            parasetslist = [x for i, x in enumerate(parasetslist) if i not in indexToRemove]
        else:
            tempParaSets2[parasetslist[0][1]] = tempParaSets[parasetslist[0][1]]
            parasetslist = []
    del tempParaSets

    ###this condenses paralogs to one gene using alignment length
    for frozenParaSets in tempParaSets2:
        maxLen = 0
        paraSets = list(frozenParaSets)
        for group in paraSets:
            if len(group) > maxLen: maxLen = len(group)
            # if 'ZNF765' in group: print('tempParaSets2', paraSets)
        if len(paraSets) > 1 and maxLen == 1:
            chimname = frozenset([list(x)[0] for x in paraSets])
            if chimname not in paraRemovedChimToReads:
                paraRemovedChimToReads[chimname] = tempParaSets2[frozenParaSets]
            else:
                paraRemovedChimToReads[chimname] = paraRemovedChimToReads[chimname] | tempParaSets2[frozenParaSets]
        elif len(paraSets) > 1:
            chimname = []
            # print(paraSets)
            for group in paraSets:
                if len(group) == 1:
                    chimname.append(list(group)[0])
                else:
                    geneToAlignLen = {}
                    for gene in group:
                        geneToAlignLen[gene] = []
                        for read in tempParaSets2[frozenParaSets]:
                            if gene in alignlen[
                                read]:  ###this is for when we have combined chim groups and not all reads will align to all paralogs in the chim
                                geneToAlignLen[gene].append(max([x[2] for x in alignlen[read][
                                    gene]]))  # use max here because this is transcriptomic alignment so we get best isoform alignment
                        geneToAlignLen[gene] = median(
                            geneToAlignLen[gene])  # median of len of all read alignments to this gene in this chimera
                    geneToAlignLen = list(geneToAlignLen.items())
                    geneToAlignLen.sort(key=lambda x: x[-1], reverse=True)
                    # print(geneToAlignLen)
                    chimname.append(geneToAlignLen[0][0])
            chimname = frozenset(chimname)
            if chimname not in paraRemovedChimToReads:
                paraRemovedChimToReads[chimname] = tempParaSets2[frozenParaSets]
            else:
                paraRemovedChimToReads[chimname] = paraRemovedChimToReads[chimname] | tempParaSets2[frozenParaSets]
    # del tempParaSets2
    return paraRemovedChimToReads
paraRemovedChimToReads = removeParalogs(chimToReads, alignlen)

readsAfterParaRemoved = 0
for c in paraRemovedChimToReads:
    readsAfterParaRemoved += len(paraRemovedChimToReads[c])
print('chim, reads after removing paralogs', len(paraRemovedChimToReads.keys()), readsAfterParaRemoved)

rejectOut = open(prefix + '-rejectedChimerasAfterParaRemoved.tsv', 'w')

genomeCloseRemovedChimToReads = {}
for chim in paraRemovedChimToReads:
    if len(paraRemovedChimToReads[chim]) > 1:
        genes = list(chim)
        # if 'ZNF765' in genes: print('paraRemovedChimToReads', genes)
        groupsGenomeDist = []
        done = set()
        # sameChrGenesSet = set()
        for g1 in genes:
            for g2 in genes:
                if g1 != g2 and frozenset((g1, g2)) not in done:
                    if g1 in genePos and g2 in genePos and genePos[g1][0] == genePos[g2][0]:
                        # sameChrGenesSet.add((g1,g2))
                        intervals = sorted([(genePos[g1][1], genePos[g1][2]), (genePos[g2][1], genePos[g2][2])])
                        if intervals[1][0] < intervals[0][1] + args.b:
                            groupsGenomeDist.append(frozenset({g1,
                                                               g2}))  ###buffer was originally just 1kb, moved default up to 50kb b/c of testing on mouse    #genomedistoverlap = True ####check if the genes are separated by at least 1kb
                        else:
                            groupsGenomeDist.append(frozenset({g1, }))
                            groupsGenomeDist.append(frozenset({g2, }))
                    else:
                        groupsGenomeDist.append(frozenset({g1, }))
                        groupsGenomeDist.append(frozenset({g2, }))
                    done.add(frozenset((g1, g2)))
        finalGroupsGenomeDist = []
        groupsGenomeDist = sorted(list(set(groupsGenomeDist)), key=len, reverse=True)
        # print('initial groups', groupsGenomeDist)
        if len(groupsGenomeDist) != len(genes):
            while len(groupsGenomeDist) > 0:
                if len(groupsGenomeDist) > 1:
                    indexToRemove = [0]
                    combSet = groupsGenomeDist[0]
                    for i in range(1, len(groupsGenomeDist)):
                        if len(combSet & groupsGenomeDist[i]) > 0:
                            combSet = combSet | groupsGenomeDist[i]
                            indexToRemove.append(i)
                    finalGroupsGenomeDist.append(combSet)
                    groupsGenomeDist = [x for i, x in enumerate(groupsGenomeDist) if i not in indexToRemove]
                else:
                    finalGroupsGenomeDist.append(groupsGenomeDist[0])
                    groupsGenomeDist = []
        else:
            finalGroupsGenomeDist = groupsGenomeDist
        # print('finalgroups', finalGroupsGenomeDist)
        if len(finalGroupsGenomeDist) > 1:
            if len(finalGroupsGenomeDist) < len(genes):
                # print('genome overlap', finalGroupsGenomeDist)
                chimname = []
                # print(paraSets)
                for group in finalGroupsGenomeDist:
                    if len(group) == 1:
                        chimname.append(list(group)[0])
                    else:
                        geneToAlignLen = {}
                        for gene in group:
                            geneToAlignLen[gene] = []
                            for read in paraRemovedChimToReads[chim]:
                                if gene in alignlen[
                                    read]:  ###this is for when we have combined chim groups and not all reads will align to all paralogs in the chim
                                    geneToAlignLen[gene].append(max([x[2] for x in alignlen[read][
                                        gene]]))  # use max here because this is transcriptomic alignment so we get best isoform alignment
                            geneToAlignLen[gene] = median(geneToAlignLen[
                                                              gene])  # median of len of all read alignments to this gene in this chimera
                        geneToAlignLen = list(geneToAlignLen.items())
                        geneToAlignLen.sort(key=lambda x: x[-1], reverse=True)
                        # print(geneToAlignLen)
                        chimname.append(geneToAlignLen[0][0])
            else:  ##no overlapping genes on genome
                chimname = genes
            chimname = frozenset(chimname)
            if chimname not in genomeCloseRemovedChimToReads:
                genomeCloseRemovedChimToReads[chimname] = paraRemovedChimToReads[chim]
            else:
                genomeCloseRemovedChimToReads[chimname] = genomeCloseRemovedChimToReads[chimname] | \
                                                          paraRemovedChimToReads[chim]
        else:
            rejectOut.write('--'.join(genes) + '\t' + 'genomeDist' + '\n')
        # else: print('genomedist', genes)
del paraRemovedChimToReads
readsAfterGenomeRemoved = 0
for c in genomeCloseRemovedChimToReads:
    readsAfterGenomeRemoved += len(genomeCloseRemovedChimToReads[c])
print('chim, reads after removing genome dist', len(genomeCloseRemovedChimToReads.keys()), readsAfterGenomeRemoved)
readsAfterFastqDistRemoved, chimAfterFastqDistRemoved = 0, 0
readsAfterReadSupRemoved, chimAfterReadSupRemoved = 0, 0

out = open(prefix + '-fusionReadCounts.tsv',
           'w')  # 'sim-nice-10x-gencode38-fusion-sim-test-06-12-2023.transcriptomeAligned-readCounts-para-genomedist-fastqdist-removed-keep-1rs-combine-chim.tsv', 'w')
out3 = open(prefix + 'genomeChunksToCut.bed', 'w')
out4 = open(prefix + 'chimericBreakpoints.tsv', 'w')
out6 = open(prefix + 'fusionWithSupportingReads.tsv', 'w')

# fusionReads = set()
readToFusion = {}
# check fastq distance between alignments
for chimname in genomeCloseRemovedChimToReads:
    genes = list(chimname)
    fastqdistoverlap = False
    fastqdistpaircomp = {}
    if len(genomeCloseRemovedChimToReads[chimname]) >= args.l:  # default read support = 3
        chimAfterReadSupRemoved += 1
        readsAfterReadSupRemoved += len(genomeCloseRemovedChimToReads[chimname])
        for g1 in genes:
            for g2 in genes:
                pair = frozenset([g1, g2])
                if g1 != g2 and pair not in fastqdistpaircomp:
                    fastqdistpaircomp[pair] = []
                    for read in genomeCloseRemovedChimToReads[chimname]:
                        intervals = []
                        for gene in pair:
                            if gene in alignlen[read]:
                                coords = [100000000000000000000000, 0]  # with chr
                                for alignment in alignlen[read][gene]:
                                    # if pair == frozenset({'SLC22A10', 'SLC22A25'}): print(read, alignment)
                                    tempcoord = sorted([alignment[0], alignment[
                                        1]])  ##to account for alignments on reverse strand where it might be higher number, smaller number
                                    if tempcoord[0] < coords[0]: coords[0] = tempcoord[0]
                                    if tempcoord[1] > coords[1]: coords[1] = tempcoord[1]
                                intervals.append(coords)
                        # if pair == frozenset({'SLC22A10', 'SLC22A25'}): print(read, 'intervals', intervals)
                        if len(intervals) > 1:
                            intervals.sort()
                            fastqdistpaircomp[pair].append(intervals[1][0] - intervals[0][1])
        for pair in fastqdistpaircomp:
            # print(pair, fastqdistpaircomp[pair])
            if len(fastqdistpaircomp[pair]) == 0:
                fastqdistoverlap = True
            elif abs(median(fastqdistpaircomp[pair])) > 20:
                fastqdistoverlap = True
        if not fastqdistoverlap:
            geneToGenomePos = {}
            geneToOuterTPos = {}
            for gene in genes:
                geneToGenomePos[gene] = {'bp': [], 'outer': []}
                geneToOuterTPos[gene] = []
            for read in genomeCloseRemovedChimToReads[chimname]:
                allBestTPos = []
                for gene in genes:
                    if gene in alignlen[read]:  ##[locs[0], locs[1], thisalignlen, tname, tstart, tend]
                        bestTAlign = (
                        [(0, 0), (0, 0)], 0, None, gene)  ##([(fstart,tstart),(fend,tend)],alignlen, tname, gene)
                        # if gene == 'FAM174C': print(gene, read, alignlen[read][gene])
                        for alignment in alignlen[read][gene]:
                            tname = alignment[3]
                            leftSSDist, rightSSDist = abs(binarySearch(intronLocs[tname], alignment[4])-alignment[4]), abs(binarySearch(intronLocs[tname], alignment[5])-alignment[5])
                            ##Combine absolute alignment length with distance from splice sites to pick best alignment
                            aligncomp = alignment[2] - (leftSSDist+rightSSDist)
                            if aligncomp > bestTAlign[1]:  # pick longest alignment
                                bestTAlign = (
                                sorted([(alignment[0], alignment[4]), (alignment[1], alignment[5])]), aligncomp,
                                alignment[3], gene)
                        if bestTAlign[2]:
                            allBestTPos.append(bestTAlign)
                if len(allBestTPos) == len(genes):  ##only consider reads that align to both of the final genes
                    allBestTPos.sort()
                    ###THIS IS ALL NOT OPTIMIZED FOR THREE GENE + FUSIONS, ASSUMES ONE SIDE OF EACH GENE IS 'OUTER'
                    for i in range(1, len(allBestTPos)):  ##for each breakpoint between two genes
                        for j in [1, 0]:  # this hits each gene independently
                            tname = allBestTPos[i - j][2]
                            gene = allBestTPos[i - j][3]
                            bpCoord = allBestTPos[i - j][0][j][1]
                            otherCoord = allBestTPos[i - j][0][1 - j][1]
                            ##next, get bp coord on transcriptome
                            closestSS = binarySearch(intronLocs[tname], bpCoord)
                            bpIntronMed = int(sum(intronToGenome[tname][closestSS]) / 2)  # (pos1, pos2)
                            # print(tname, closestSS, bpIntronMed, genePos[gene])

                            ##this if statement means read only counted if within gene boundaries
                            if (closestSS != 0 and closestSS != intronLocs[tname][-1]) or (genePos[gene][1] <= bpIntronMed <= genePos[gene][2]):
                                geneToGenomePos[gene]['bp'].append(bpIntronMed)
                                if otherCoord > bpCoord:
                                    outerMed = int(sum(intronToGenome[tname][max(intronLocs[
                                                                                     tname])]) / 2)  ###intronLocs[tname][0] is 0, intronLocs[tname][-1] is end of transcript
                                    geneToOuterTPos[gene].append(max(intronLocs[tname]))
                                else:
                                    outerMed = int(sum(intronToGenome[tname][min(intronLocs[tname])]) / 2)
                                    geneToOuterTPos[gene].append(min(intronLocs[tname]))
                                geneToGenomePos[gene]['outer'].append(outerMed)
                            # else: print(tname, intronLocs[tname], closestSS, bpIntronMed, genePos[gene])
            temp = [len(y) for x, y in geneToOuterTPos.items()]
            if min(temp) > 1:
                geneOrder = sorted([(median(geneToOuterTPos[g]), g) for g in geneToOuterTPos.keys()])
                fusionname = '--'.join([x[1] for x in geneOrder])
                chimoutlines, bedoutlines = '', ''
                bpatendofgene = False

                # print(geneToGenomePos)

                for i in range(len(geneOrder)):
                    gene = geneOrder[i][1]
                    if median(geneToGenomePos[gene]['bp']) > median(geneToGenomePos[gene]['outer']) or sum(
                            geneToGenomePos[gene]['bp']) / len(geneToGenomePos[gene]['bp']) > sum(
                            geneToGenomePos[gene]['outer']) / len(geneToGenomePos[gene]['outer']):
                        coord = [max(geneToGenomePos[gene]['bp']), min(geneToGenomePos[gene]['outer'])]
                    else:
                        coord = [min(geneToGenomePos[gene]['bp']), max(geneToGenomePos[gene]['outer'])]
                    # if genePos[gene][1]-100 <= coord[0] <= genePos[gene][1]+100 or genePos[gene][2]-100 <= coord[0] <= genePos[gene][2]+100:
                    #     bpatendofgene = True
                    #     break
                    # else:
                    # fusionName	geneName	orderInFusion	geneChr	breakpointCoord	outerEdgeCoord	readSupport
                    fusionend = "5'gene" if i == 0 else "3'gene"

                    chimoutlines += '\t'.join(
                        [fusionname, gene, fusionend, genePos[gene][0], str(coord[0]), str(coord[1]),
                         str(len(genomeCloseRemovedChimToReads[chimname]))]) + '\n'
                    coord.sort()

                    ###DONT LIKE SAVING END AS TEXT, SHOULD JUST SAVE GENE ORDER AS NUMBERS - FIRST GENE IS 0, etc

                    bedoutlines += '\t'.join(
                        [genePos[gene][0], str(coord[0]), str(coord[1]), fusionname + '-.-' + gene]) + '\n'

                    # for i in range(len(geneToGenomePos[gene]['bp'])):
                    #     coord = sorted([geneToGenomePos[gene]['bp'][i], geneToGenomePos[gene]['outer'][i]])
                    #     out3.write('\t'.join([genePos[gene][0], str(coord[0]), str(coord[1]), '-'.join(genes) + '--' + gene]) + '\n')
                # if not bpatendofgene:
                out4.write(chimoutlines)
                out3.write(bedoutlines)

                chimAfterFastqDistRemoved += 1
                readsAfterFastqDistRemoved += len(genomeCloseRemovedChimToReads[chimname])

                out.write('-'.join(genes) + '\t' + str(len(genomeCloseRemovedChimToReads[chimname])) + '\n')
                # fusionReads = fusionReads | genomeCloseRemovedChimToReads[chimname]
                out6.write('-'.join(genes) + '\t' + ','.join(genomeCloseRemovedChimToReads[chimname]) + '\n')
                for i in genomeCloseRemovedChimToReads[chimname]: readToFusion[i] = '-'.join(genes)
            else:
                rejectOut.write('--'.join(genes) + '\t' + 'edgeOfGene' + '\n')
        else:
            rejectOut.write('--'.join(genes) + '\t' + 'fastqDist' + '\n')
    else:
        rejectOut.write('--'.join(genes) + '\t' + 'readSup' + '\n')
        # print('fastqdist', genes)
print('chim, reads after removing low read support', chimAfterReadSupRemoved, readsAfterReadSupRemoved)
print('chim, reads after removing fastq dist', chimAfterFastqDistRemoved, readsAfterFastqDistRemoved)

# print(readToFusion)

if chimAfterFastqDistRemoved > 0:
    if args.r.split('.')[-1] == 'fastq' or args.r.split('.')[-1] == 'fq':
        out5 = open(prefix + '-fusionOnly.fastq', 'w')
        last = False
        for line in open(args.r):
            if line[0] == '@':
                if line.lstrip('@').split(' ')[0] in readToFusion:
                    last = True
                else:
                    last = False
            if last: out5.write(line)
    elif args.r.split('.')[-1] == 'fasta' or args.r.split('.')[-1] == 'fa':
        out5 = open(prefix + '-fusionOnly.fasta', 'w')
        last = False
        for line in open(args.r):
            if line[0] == '>':
                if line.rstrip().lstrip('>').split(' ')[0] in readToFusion:
                    last = True
                else:
                    last = False
            if last: out5.write(line)


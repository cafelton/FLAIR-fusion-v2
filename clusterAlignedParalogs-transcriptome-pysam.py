import sys
from collections import Counter
import pysam

###for mapping gencode reference simulated by badreads to genome to figure out spurious chimeras
##first convert badreads .fastq file with random read naming to name read by reference
###ENST00000387460.2|ENSG00000210195.2|-|-|MT-TT-201|MT-TT|66|Mt_tRNA|

###https://stackoverflow.com/questions/10301000/how-to-find-connected-components
def get_all_connected_groups(graph):
    already_seen = set()
    result = []
    for node in graph:
        if node not in already_seen:
            connected_group, already_seen = get_connected_group(graph, node, already_seen)
            result.append(connected_group)
    return result


def get_connected_group(graph, node, already_seen):
        result = []
        nodes = set([node])
        while nodes:
            node = nodes.pop()
            already_seen.add(node)
            # nodes = nodes or graph[node] - already_seen
            nodes.update(graph[node] - already_seen)
            result.append(node)
        return result, already_seen

geneLen = {}
# for line in open('/private/groups/brookslab/cafelton/fusions-code/gencode.v38.annotation-short.gtf'):
#     if line[0] != '#':
#         line = line.split('\t')
#         gene_name = line[-1].split('; gene_name "')[1].split('"')[0]
#         geneLen[gene_name] = abs(int(line[4])-int(line[3]))

lastgenelen = 0
tToG = {}
for line in open(sys.argv[1]): ###'/private/groups/brookslab/cafelton/fusions-code/gencode.v38.annotation-short.gtf'):
    if line[0] != '#':
        line = line.split('\t')
        if line[2] == 'gene':
            gene_name = line[-1].split('; gene_name "')[1].split('"')[0]
            geneLen[gene_name] = abs(int(line[4])-int(line[3]))
            lastgenelen = abs(int(line[4])-int(line[3]))
        if line[2] == 'transcript':
            t_name = line[-1].split('; transcript_name "')[1].split('"')[0]
            gene_name = line[-1].split('; gene_name "')[1].split('"')[0]
            tToG[t_name] = gene_name


# readMap = {}
# last = False
# for line in open(sys.argv[2]):###'simulatedGencode38-100x-avg.fastq'):
#     if line[0] == '@':
#         line = line.split(' ')
#         readname = line[0].lstrip('@')
#         if len(line) > 1 and line[1] != 'junk_seq' and line[1] != 'random_seq' and (len(line) > 2 and line[2] != 'chimera'): ##remove simulated chimeric reads, annotated as @43ab8b0e-b3d1-eb19-424e-a394b477e551 CCDC6-201--RET-202,-strand,0-3716 chimera EML4-201--ALK-205,-strand,0-3048 length=6390 error-free_length=6812 read_identity=82.52%
#                             ##normal = @d1374e79-ddca-6b71-88ee-08363a983c43 ENST00000481672.5|ENSG00000076344.16|OTTHUMG00000064893.7|OTTHUMT00000139331.1|RGS11-207|RGS11|772|retained_intron|,+strand,0-772 length=743 error-free_length=796 read_identity=81.15%
#             # print(line)
#             gname = line[1].split('|')[5]
#             # line[0] += '--' + tname
#             # out.write(' '.join(line))
#             last = readname
#             readMap[readname] = [gname]
#         else:
#             readMap[readname] = 'bad'
#             last = False
#             # rejectFile.write(' '.join(line))
#     elif last:
#         readMap[last].append(len(line))
#     #     out.write(line)
#     # else: rejectFile.write(line)
# # print(Counter(rejects))


# print("read in readMap")
gene_graph = {}
edge_weights = {}
align_count = {}
###dict[correct_GENE]: {aligned_GENE_1:countOfAlignments, aligned_GENE_2:c, etc}
##2665564e-81bb-724e-7c78-58a2a32a19d9    272     ENST00000259470.6|ENSG00000136943.12|OTTHUMG00000020314.3|OTTHUMT00000053301.3|CTSV-201|CTSV|4359|protein_coding|
###need to cluster paralogs somehow


##>Gm26206-201--len107--ident92.08%--63
samfile = pysam.AlignmentFile(sys.argv[2], "rb")#"sim-avg-10x-gencode38-fusion-sim-test-06-12-2023.transcriptomeAligned.sorted.bam", "rb")
for s in samfile:
    if s.is_mapped:
        readname = s.query_name
        geneinfo = s.reference_name.split('|')
        alignGene = geneinfo[5]
        trueGene = readname.split('--')[0]
        #if '-' in trueGene: trueGene = '-'.join(trueGene.split('-')[:-1])
        trueGene = tToG[trueGene]
        #print(readname, trueGene, geneinfo)
        readlen = s.infer_read_length()
        if readlen > 350 or readlen > 0.8 * geneLen[trueGene]:
            if trueGene != alignGene:
                edgeName = frozenset([trueGene, alignGene])
                if edgeName not in edge_weights: edge_weights[edgeName] = 0
                edge_weights[edgeName] += 1

            if trueGene not in gene_graph: gene_graph[trueGene] = {trueGene}
            # if trueGene not in align_count: align_count[trueGene] = {}
            # if alignGene not in align_count[trueGene]: align_count[trueGene][alignGene] = 0
            # align_count[trueGene][alignGene] += 1
            gene_graph[trueGene].add(alignGene)
            if alignGene not in gene_graph: gene_graph[alignGene] = {alignGene}

#
# for line in open('simulatedGencode38-100x-avg.alignedTranscriptome-full.sam'):
#     if line[0] != '@':
#         line = line.split('\t')
#         if line[2] != '*':
#             readname, alignGene = line[0], line[2].split('|')[5]
#             if readMap[readname] != 'bad': #and alignGene[:3] != 'chr':
#                 trueGene = readMap[readname][0]
#                 if readMap[readname][1] > 350 or readMap[readname][1] > 0.8 * geneLen[trueGene]:
#                     if trueGene != alignGene:
#                         edgeName = frozenset([trueGene, alignGene])
#                         if edgeName not in edge_weights: edge_weights[edgeName] = 0
#                         edge_weights[edgeName] += 1
#
#                     if trueGene not in gene_graph: gene_graph[trueGene] = {trueGene}
#                     # if trueGene not in align_count: align_count[trueGene] = {}
#                     # if alignGene not in align_count[trueGene]: align_count[trueGene][alignGene] = 0
#                     # align_count[trueGene][alignGene] += 1
#                     gene_graph[trueGene].add(alignGene)
#                     if alignGene not in gene_graph: gene_graph[alignGene] = {alignGene}
print("made gene graph")


# count_genes_with_x_alignments = {}
# totAlign, correctAlign = 0, 0
# correct_align_frac_by_num_genes_aligned_to = {}
# for trueGene in align_count:
#     tot, corr = 0, 0
#     for alignGene in align_count[trueGene]:
#         if alignGene == trueGene:
#             correctAlign += align_count[trueGene][alignGene]
#             corr += align_count[trueGene][alignGene]
#         totAlign += align_count[trueGene][alignGene]
#         tot += align_count[trueGene][alignGene]
#     alignNum = len(align_count[trueGene])
#     if alignNum not in correct_align_frac_by_num_genes_aligned_to: correct_align_frac_by_num_genes_aligned_to[alignNum] = []
#     if alignNum not in count_genes_with_x_alignments: count_genes_with_x_alignments[alignNum] = 0
#     count_genes_with_x_alignments[alignNum] += 1
#     correct_align_frac_by_num_genes_aligned_to[alignNum].append(round(float(corr)/tot, 3))
# print('count_genes_with_x_alignments')
# out = open("simulatedGencode38-100x-avg.Transcriptome-count_genes_with_x_alignments.csv", 'w')
# for i in range(max(count_genes_with_x_alignments)):
#     if i+1 in count_genes_with_x_alignments: out.write(str(i+1) + ',' + str(count_genes_with_x_alignments[i+1]) + '\n')
#     else: out.write(str(i+1) + ',' + '0' + '\n')
# out.close()
# print('fraction of correct alignments', correctAlign, totAlign, round(float(correctAlign)/totAlign, 3))
# out = open("simulatedGencode38-100x-avg.Transcriptome-correct_align_frac_by_num_genes_aligned_to.txt", 'w')
# for i in correct_align_frac_by_num_genes_aligned_to:
#     out.write(str(i) + '\t' + ','.join([str(x) for x in correct_align_frac_by_num_genes_aligned_to[i]]) + '\n')
# out.close()

# paralog_clusters = get_all_connected_groups(gene_graph)
# print("got connected groups")
# ##checking if genes are showing up multiple times in the clustering - not sure if this is actually a problem yet tho
# # print(Counter([item for sublist in paralog_clusters for item in sublist]).most_common()[:20])
# ###seems like not
#
# paralog_clusters.sort(key=len, reverse=True)
#
# out = open("simulatedGencode38-100x-avg.TranscriptomeEdgeGraph-biggestcluster.csv", 'w')
# biggestCluster = set(paralog_clusters[0])
# for edge in edge_weights:
#     e2 = list(edge)
#     if e2[0] in biggestCluster:
#         out.write(e2[0] + ',' + e2[1] + ',' + str(edge_weights[edge]) + '\n')
# out.close()

out = open(sys.argv[2].split('/')[-1].split('.')[0] +  "TranscriptomeGeneToNeighbors-filteredReadLen.tsv", 'w')
for gene in gene_graph:
    if len(gene_graph[gene]) > 1:
        out.write(gene + '\t' + ','.join(list(gene_graph[gene]-{gene})) + '\n')
out.close()


# out = open("simulatedGencode38-100x-avg.TranscriptomeParalogClusters.txt", 'w')
# for clust in paralog_clusters:
#     if len(clust) > 1:
#         out.write('\t'.join(clust) + '\n')
# out.close()

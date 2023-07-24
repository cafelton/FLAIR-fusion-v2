
import sys
###get intron chain for isoforms of interested genes (genomic)
###transform coordinates of matching region into transcriptomic coordinates
###get ends of sam alignments to transcript, then get set of isoforms that match that region

# genes = {}
# for line in open('fusion-genes-drr.txt'):
#     genes[line.rstrip()] = {}
genes = {}
transcripts = {}
transcriptlen = {}
for line in open(sys.argv[1]):#'gencode.vM32.primary_assembly.annotation.gtf'): ##"/private/groups/brookslab/reference_annotations/gencode.v37.annotation.gtf"):
    if line[0] != '#':
        line = line.split('\t')
        if line[2] == 'exon':# and genename in genes:
            genename = line[8].split('; gene_name "')[1].split('"')[0]
            isoname = line[8].split('; transcript_name "')[1].split('"')[0]
            if genename not in genes: genes[genename] = {}
            if isoname not in transcriptlen: transcriptlen[isoname] = 0
            transcriptlen[isoname] += int(line[4]) - int(line[3])
            if isoname not in genes[genename]: genes[genename][isoname] = []
            else:
                if line[6] == '+': genes[genename][isoname].append((last, int(line[3])))
                else: genes[genename][isoname].append((int(line[4]), last)) ##editied, was .insert(0, before
            last = int(line[4]) if line[6] == '+' else int(line[3])
        elif line[2] == 'transcript': #and genename in genes:
            transcripts[line[8].split('; transcript_name "')[1].split('"')[0]] = (line[6], int(line[3]), int(line[4]), line[0])
# intron_nodes = {}
# # print(genes)
# for g in genes:
#     intron_nodes[g] = {}
#     for iso in genes[g]:
#         for node in genes[g][iso]:
#             if node not in intron_nodes[g]: intron_nodes[g][node] = []
#             intron_nodes[g][node].append(iso)
out = open('transcriptome_introns_to_genome_coords_' + '.'.join(sys.argv[1].split('/')[-1].split('.')[:-1]) + '.tsv', 'w') #gencode.vM32.primary_assembly.tsv', 'w')    #'/private/groups/brookslab/cafelton/fusions-code/FLAIR-fusion-v2.0/transcriptome_introns_to_genome_coords_gencode37.tsv', 'w')
for g in genes:
    for iso in genes[g]:
        coordlist = []
        isocoord = 0
        if transcripts[iso][0] == '+':
            genomecoord = transcripts[iso][1]
            ###allow for start of transcript as acceptable intron
            coordlist.append('.'.join([str(x) for x in [isocoord, transcripts[iso][1]-500, transcripts[iso][1]]]))
            for intron in genes[g][iso]:
                isocoord += intron[0]-genomecoord
                genomecoord = intron[1]
                coordlist.append('.'.join([str(isocoord)] + [str(x) for x in intron]))
            coordlist.append('.'.join([str(x) for x in [transcriptlen[iso], transcripts[iso][2], transcripts[iso][2] + 500]]))
        else:
            genomecoord = transcripts[iso][2]
            coordlist.append('.'.join([str(x) for x in [isocoord, transcripts[iso][2], transcripts[iso][2] + 500]]))
            for intron in genes[g][iso]:
                isocoord += genomecoord-intron[1]
                genomecoord = intron[0]
                coordlist.append('.'.join([str(isocoord)] + [str(x) for x in intron]))
            coordlist.append(
                '.'.join([str(x) for x in [transcriptlen[iso], transcripts[iso][1]-500, transcripts[iso][1]]]))
        out.write('\t'.join([g, iso, transcripts[iso][3], ','.join(coordlist)]) + '\n')
out.close()

# genes = {}
# for line in open('fusion-genes-drr.txt'):
#     genes[line.rstrip()] = {}
# transcripts = {}
# for line in open("/private/groups/brookslab/reference_annotations/gencode.v38.annotation.gtf"):
#     if line[0] != '#':
#         line = line.split('\t')
#         genename = line[8].split('; gene_name "')[1].split('"')[0]
#         if line[2] == 'exon' and genename in genes:
#             isoname = line[8].split('; transcript_name "')[1].split('"')[0]
#             if isoname not in genes[genename]: genes[genename][isoname] = []
#             else:
#                 if line[6] == '+': genes[genename][isoname].append((last, int(line[3])))
#                 else: genes[genename][isoname].append((int(line[4]), last)) ##editied, was .insert(0, before
#             last = int(line[4]) if line[6] == '+' else int(line[3])
#         elif line[2] == 'transcript' and genename in genes:
#             transcripts[line[8].split('; transcript_name "')[1].split('"')[0]] = (line[6], int(line[3]), int(line[4]))
# intron_nodes = {}
# print(genes)
# for g in genes:
#     intron_nodes[g] = {}
#     for iso in genes[g]:
#         for node in genes[g][iso]:
#             if node not in intron_nodes[g]: intron_nodes[g][node] = []
#             intron_nodes[g][node].append(iso)
# out = open('DRR-fusion-genes-introns-to-transcriptome-coords.txt', 'w')
# for g in genes:
#     for iso in genes[g]:
#         out.write('\t'.join([iso, transcripts[iso][0], transcripts[iso][1], transcripts[iso][2]]) + '\n')
#         if transcripts[iso][0] == '+':
#             genomecoord = transcripts[iso][1]
#             isocoord = 0
#             for intron in genes[g][iso]:
#                 isocoord += intron[0]-genomecoord
#                 genomecoord = intron[1]
#                 out.write('\t'.join([iso, str(isocoord), ','.join([str(x) for x in intron])]) + '\n')
#         else:
#             genomecoord = transcripts[iso][2]
#             isocoord = 0
#             for intron in genes[g][iso]:
#                 isocoord += genomecoord-intron[1]
#                 genomecoord = intron[0]
#                 out.write('\t'.join([iso, str(isocoord), ','.join([str(x) for x in intron])]) + '\n')
#
#
#
#
#
# for line in open("DRR-aligned-chimera-no-second-all.sam"):
#     if line[0] != '@':
#         line = line.split('\t')
#         iso = line[2].split('|')[4]
#         if iso in intron_nodes: #and int(line[4]) ==60:
#             start = int(line[3])
#             end = int(line[3])
#             cigar = line[5]
#             last = 0
#             while last < len(cigar):
#                 j = last+1
#                 while cigar[j].isnumeric():
#                     j += 1
#                 if cigar[j] in ['M', 'D', 'N', '=', 'X']:
#                     end += int(cigar[last:j])
#                 last = j+1
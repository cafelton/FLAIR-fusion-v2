
import sys

###TEKT2,chr1,36085562,36084093--CTAG1A,chrX,154586406,154586816	0	1879	TEKT2-201--CTAG1A-202_CTAG1A--TEKT2	1000	+	0	1879	0	4	56,208,126,291,	0,776,1069,1588,

##LETM1,chr4,1833896,1856156--USP2,chr11,119377497,119355214	0	44543	LETM1-201--USP2-201_LETM1--USP2	1000	+	0	44543	0	17	288,61,451,144,138,815,51,124,112,111,65,104,81,79,108,121,1708,	0,6947,14359,19584,21174,26236,39523,40097,40415,40623,40920,41505,41696,41922,42167,42450,42835,
# print('hi')
#'.'.join(sys.argv[1].split('.')[:-5]) +

isoreadsup = {}
freadsfinal = set()
for line in open('.'.join(sys.argv[1].split('.')[:-4]) + '.syntheticAligned-flair.collapse.combined.isoform.read.map.txt'):
    line = line.split('\t')
    isoreadsup[line[0]] = len(line[1].split(','))
    if len(line[1].split(',')) > 1:
        for i in line[1].split(','):
            freadsfinal.add(i)
            

if sys.argv[2].split('.')[-1] == 'fastq' or sys.argv[2].split('.')[-1] == 'fq':
    out5 = open('.'.join(sys.argv[2].split('.')[:-1]) + '-isoSupport.fastq', 'w')
    last = False
    for line in open(sys.argv[2]):
        if line[0] == '@':
            if line.lstrip('@').split(' ')[0] in freadsfinal:
                last = True
            else:
                last = False
        if last: out5.write(line)
elif sys.argv[2].split('.')[-1] == 'fasta' or sys.argv[2].split('.')[-1] == 'fa':
    out5 = open('.'.join(sys.argv[2].split('.')[:-1]) + '-isoSupport.fasta', 'w')
    last = False
    for line in open(sys.argv[2]):
        if line[0] == '>':
            if line.rstrip().lstrip('>').split(' ')[0] in freadsfinal:
                last = True
            else:
                last = False
        if last: out5.write(line)



out = open('.'.join(sys.argv[1].split('.')[:-4]) + '.genomeAligned-flair.collapse.isoforms.bed', 'w') #'sim-nice-10x-gencode38-fusion-sim-test-06-12-2023corrleft-fusionOnly.genomeAligned.flair.collapse.isoforms.bed', 'w')
for line in open(sys.argv[1]): #'sim-nice-10x-gencode38-fusion-sim-test-06-12-2023corrleft-fusionOnly.syntheticAligned-flair.collapse.isoforms.bed'):
    line = line.split('\t')
    iso, start, esizes, estarts = line[3], int(line[1]), [int(x) for x in line[10].split(',')[:-1]], [int(x) for x in line[11].split(',')[:-1]]
    if isoreadsup[iso] > 1:
        chimera = line[0].split('--')
        fivechr, fivebp, fiveouter = chimera[0].split('.')[-3:]
        threechr, threebp, threeouter = chimera[1].split('.')[-3:]
        fivebp, fiveouter, threebp, threeouter = int(fivebp), int(fiveouter), int(threebp), int(threeouter)
        breakpoint = abs(fiveouter-fivebp)
        # outline5 = [fivechr, None, None, iso, 1000, None, None, None, 0, 0, [], []]

        ###check if isoform actually crosses fusion breakpoint, don't convert coordinates otherwise
        if start < breakpoint and breakpoint < int(line[2]):
            introns5, exons5 = [], []
            introns3, exons3 = [0], []
            lastexonend = 0
            start3 = None
            for i in range(len(esizes)):
                if start + estarts[i] < breakpoint:
                    introns5.append(estarts[i]-lastexonend)
                    exons5.append(esizes[i])
                    lastexonend = estarts[i] + esizes[i]
                else:
                    if not start3: start3 = (start + estarts[i])-breakpoint
                    else: introns3.append(estarts[i]-lastexonend)
                    exons3.append(esizes[i])
                    lastexonend = estarts[i] + esizes[i]
            #print(iso, start3, exons3, introns3)
            if fivebp > fiveouter: ##5' gene is + direction
                tot5len, g5estarts = 0, []
                for i in range(len(exons5)):
                    g5estarts.append(tot5len + introns5[i])
                    tot5len += introns5[i] + exons5[i]
                outline5 = [fivechr, str(fiveouter+start), str(fiveouter+start+tot5len), iso, '1000', '+',
                            str(fiveouter+start), str(fiveouter+start+tot5len), '0',
                            str(len(exons5)), ','.join([str(x) for x in exons5]), ','.join([str(x) for x in g5estarts])]
            else: #5' gene is in - direction
                tot5len, g5estarts = 0, []
                for i in range(len(exons5)-1, -1, -1): #loop backwards through gene
                    g5estarts.append(tot5len)
                    tot5len += introns5[i] + exons5[i]
                outline5 = [fivechr, str(fiveouter - (start+tot5len)), str(fiveouter - start), iso, '1000', '-',
                            str(fiveouter - (start+tot5len)), str(fiveouter - start), '0',
                            str(len(exons5)), ','.join([str(x) for x in exons5[::-1]]), ','.join([str(x) for x in g5estarts])]
            if threebp < threeouter: #3' gene is in + direction
                tot3len, g3estarts = 0, []
                for i in range(len(exons3)):
                    g3estarts.append(tot3len + introns3[i])
                    tot3len += introns3[i] + exons3[i]
                outline3 = [threechr, str(threebp + start3), str(threebp + start3 + tot3len), iso, '1000', '+',
                            str(threebp + start3), str(threebp + start3 + tot3len), '0',
                            str(len(exons3)), ','.join([str(x) for x in exons3]), ','.join([str(x) for x in g3estarts])]
            else: #3' gene is in - direction
                tot3len, g3estarts = 0, []
                for i in range(len(exons3)-1, -1, -1): #loop backwards through gene
                    g3estarts.append(tot3len)
                    tot3len += introns3[i] + exons3[i]
                outline3 = [threechr, str(threebp - (start3+tot3len)), str(threebp - start3), iso, '1000', '-',
                            str(threebp - (start3+tot3len)), str(threebp - start3), '0',
                            str(len(exons3)), ','.join([str(x) for x in exons3[::-1]]), ','.join([str(x) for x in g3estarts])]
            out.write('\t'.join(outline5) + '\n')
            out.write('\t'.join(outline3) + '\n')
out.close()
#
# out = open('test.txt', 'w')
# out.write('hi')

import random
import scipy.stats as stats
import sys



lower, upper = 50, 99
mean, stdev = 95,2.5
readidentdist = stats.truncnorm(
    (lower - mean) / stdev, (upper - mean) / stdev, loc=mean, scale=stdev)

lower, upper = 10, 100000
mean, stdev = 1500,1300
lendist = stats.truncnorm(
    (lower - mean) / stdev, (upper - mean) / stdev, loc=mean, scale=stdev)


nuc = {'A', 'C', 'T', 'G'}
lastname, seq = None, ''
c, d = 0, 0
if len(sys.argv[2]) > 0 and sys.argv[2][-1] != '/': sys.argv[2] += '/'

out = open(sys.argv[2] + 'sim-avg-100x-' + sys.argv[1].split('/')[-1], 'w')#'mysim-avg-100x-gencode.vM32.transcripts.fa', 'w')
for line in open(sys.argv[1]): #'/private/groups/brookslab/cafelton/Vollmers-mouse-R10-R2C2-demultiplexed/gencode.vM32.transcripts.fa'):
    if line[0] == '>':
        if lastname:
            c += 1
            if c%1494 == 0: print(c/1494, '% done')
            theselens = lendist.rvs(100)
            theseidents = readidentdist.rvs(100)
            for i in range(100):
                d += 1
                mylen = int(theselens[i])
                myident = theseidents[i]
                if mylen >= len(seq):
                    subtrans = seq
                else:
                    startpos = random.randint(0, len(seq) - mylen + 1)
                    # print(len(seq), mylen, startpos, startpos + mylen)
                    subtrans = seq[startpos:startpos + mylen]
                subtrans = list(subtrans)
                # for j in random.sample(range(len(subtrans)), int(len(subtrans) * ((100 - myident)/100))):
                #     subtrans[j] = random.choice(list(nuc - {subtrans[j]}))
                ###error rate for substitutions, insertions, and deletions are about equal https://www.nature.com/articles/s41467-020-20340-8#Fig1
                for j in range(int(len(subtrans) * ((100 - myident)/100))):
                    chartochange = random.randint(0, len(subtrans)-1)
                    typeoferror = random.randint(0,2)
                    if typeoferror == 0: subtrans.pop(chartochange) ##deletion
                    elif typeoferror == 1: subtrans.insert(chartochange, random.choice(list(nuc))) ##insertion
                    else: subtrans[chartochange] = random.choice(list(nuc - {subtrans[chartochange]}))
                out.write('>' + lastname + '--len' + str(len(subtrans)) + '--ident' + str(round(myident, 2)) + '%--' + str(d) + '\n')
                out.write(''.join(subtrans) + '\n')
        lastname = line.split('|')[4]
        seq = ''
    else: seq += line.rstrip()



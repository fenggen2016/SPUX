###Part 1: Combines paired reads####
###Part 2: Finds cdrs, combines them into one sequence and checks Phred quality###
###Part 3: Translates CDR sequences
###Part 4: Calculates percentage of mutated CDRs and number of mutations per gene using the "perfect genes"###
###Part 5: Calculates the aminacid variation at each positions using###

import csv

quality_control = raw_input('Quality:') # Choose Phred quality threshold.
quality_control = int(quality_control)


#Quality dictionary
phred = {
    '!': 0, ',': 11, '7': 22, 'B': 33,
    '"': 1, '-': 12, '8': 23, 'C': 34,
    '#': 2, '.': 13, '9': 24, 'D': 35,
    '$': 3, '/': 14, ':': 25, 'E': 36,
    '%': 4, '0': 15, ';': 26, 'F': 37,
    '&': 5, '1': 16, '<': 27, 'G': 38,
    '': 6,  '2': 17, '=': 28, 'H': 39,
    '(': 7, '3': 18, '>': 29, 'I': 40,
    ')': 8, '4': 19, '?': 30, 'J': 41,
    '*': 9, '5': 20, '@': 31, 'K': 42,
    '+': 10, '6': 21, 'A': 32,

}


###Part 1: Combines paired reads####

def comp(sequence):
    '''Function that takes a sequence and complements it'''

    # define translation dictionary
    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'N': 'N',
        'R': 'Y',
        'Y': 'R',
        'K': 'M',
        'M': 'K',
        'B': 'V',
        'V': 'B',
        'D': 'H',
        'H': 'D',
    }

    # translate each nucleotide in the sequence
    complementary_sequence = "".join(
        [complement.get(nucleotide.upper(), '') for nucleotide in sequence]
    )

    return complementary_sequence


def revcomp(sequence):
    '''Function that takes a sequence and reversecomplements it'''

    # use the comp function to get complement
    complementary = comp(sequence)

    # reverse the sequence
    reverse_complementary = complementary[::-1]

    return reverse_complementary


def convert_fastqs(in_file_1, in_file_2):
    """
    This function reads a pair of gzipped fastqfiles and yields strings on fasta format
    the fasta format consist of the fastq read header and then the read one sequence followed by a stretch of N bases and then the reverse complement of the read two sequence
    the function is totally unavare of base calling qualities so some kind of filtering based on this might be needed later on to mask or remove low quality base calls
    """

    count = 0
    # a looop that reads all lines in the infiles
    while True:

        try:
            # reading the read one file
            header_1 = in_file_1.readline().rstrip()
            sequence_1 = in_file_1.readline().rstrip()
            junk_1 = in_file_1.readline().rstrip()
            quality_1 = in_file_1.readline().rstrip()

            # reading the read two file
            header_2 = in_file_2.readline().rstrip()
            sequence_2 = in_file_2.readline().rstrip()
            junk_2 = in_file_2.readline().rstrip()
            quality_2 = in_file_2.readline().rstrip()

            # check if we reached the end of the file
            if '' in [header_1, sequence_1, junk_1, quality_1, header_2, sequence_2, junk_2, quality_2]: raise EOFError

        # if we reached the end of file break the loop
        except EOFError:
            break

        # use this many N bases between the read sequences
        length_of_N_stretch = 20

        # make the long sequence that consist of both reads
        long_sequence = sequence_1 + 'N'.join('' for i in xrange(length_of_N_stretch)) + revcomp(sequence_2)
        
        # make the long sequence that consist of both qualities
        long_quality = quality_1 + 'N'.join('' for i in xrange(length_of_N_stretch)) + quality_2[::-1]

        # yield the long sequence on fasta format and continue with next iteration in the loop
        yield long_sequence + long_quality
        count = count + 1
        print count


#
# imports
#
import sys
import gzip

#
# open the infiles
#
in_file_1 = gzip.open(sys.argv[1])
in_file_2 = gzip.open(sys.argv[2])

#
# convert the infiles and print to stdout
#

###Part 2: Finds cdrs, combines them into one sequence and checks Phred quality###

ofile = open('cdrs_all.csv', 'w')
writer = csv.writer(ofile)
pfile = open('cdrs_correct.csv', 'w')
priter = csv.writer(pfile)

perfect_gene = 0
non_perfect_gene = 0

h1_start = 'CACCTTT'
h1_stop = 'TGGGTCC'

h1_good = 0
h1_bad = 0


h2_start = 'GGTCTCA'
h2_stop = 'TATGCAG'

h2_good = 0
h2_bad = 0

h3_start = 'TGCGCGC'
h3_stop = 'GACTATT'

h3_good = 0
h3_bad = 0

l1_start = 'CTATGGT'
l1_stop = 'CGTTCAC'

l1_good = 0
l1_bad = 0

l2_start = 'CCTGGGG'
l2_stop = 'TATCTAG'

l2_good = 0
l2_bad = 0

l3_start = 'GTTTTCA'
l3_stop = 'GACAACT'

l3_good = 0
l3_bad = 0

perfect_count = 0
some_fault = 0
super_count = 0


for long_sequence in convert_fastqs(in_file_1, in_file_2):
    
    if long_sequence.find(h1_stop) - long_sequence.find(h1_start) == 25:
        h1_read = long_sequence[long_sequence.find(h1_start) + 7:long_sequence.find(h1_stop)]
        
        h1_good = h1_good + 1
        h1_quality = long_sequence[long_sequence.find(h1_start) + len(long_sequence)/2 + 7:long_sequence.find(h1_stop) + + len(long_sequence)/2]
        h1 = ''
        count = 0
        for i in h1_quality:
            if phred[i] > quality_control:
                h1 = h1 + h1_read[count]
            else:
                h1 = h1 + 'X'

            count = count + 1

#print h1, h1_quality
    
    
    else:
        h1 = 'XXXXXXXXXXXXXXXXXX'
        h1_bad = h1_bad + 1
    #print h1
    
    if long_sequence.find(h2_stop) - long_sequence.find(h2_start) == 37:
        
        h2_read = long_sequence[long_sequence.find(h2_start) + 7:long_sequence.find(h2_stop)]
        h2_good = h2_good + 1
        h2_quality = long_sequence[long_sequence.find(h2_start) + len(long_sequence)/2 + 7:long_sequence.find(h2_stop) + + len(long_sequence)/2]
        h2 = ''
        count = 0
        for i in h2_quality:
            if phred[i] > quality_control:
                h2 = h2 + h2_read[count]
            else:
                h2 = h2 + 'X'

            count = count + 1
#print h2, h2_quality


    #print h2
    
    else:
        h2 = 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
        h2_bad = h2_bad + 1
    #print h2
    
    if long_sequence.find(h3_stop) - long_sequence.find(h3_start) == 37:
        h3_read = long_sequence[long_sequence.find(h3_start) + 7:long_sequence.find(h3_stop)]
        h3_good = h3_good + 1
        h3_quality = long_sequence[long_sequence.find(h3_start) + len(long_sequence)/2 + 7:long_sequence.find(h3_stop) + + len(long_sequence)/2]
        h3 = ''
        count = 0
        for i in h3_quality:
            if phred[i] > quality_control:
                h3 = h3 + h3_read[count]
            else:
                h3 = h3 + 'X'
            
            count = count + 1

#print h3, h3_quality




    else:
        h3 = 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
        h3_bad = h3_bad + 1
    #print h3




    reverse_sequence = long_sequence[::-1]

#print reverse_sequence.find(l1_start)
#print reverse_sequence.find(l1_stop)

    if reverse_sequence.find(l1_stop) - reverse_sequence.find(l1_start) == 40:
        l1_good = l1_good + 1
        l1 = reverse_sequence[reverse_sequence.find(l1_start) + 7:reverse_sequence.find(l1_stop)]
        l1_quality = reverse_sequence[reverse_sequence.find(l1_start) - len(reverse_sequence)/2 + 7:reverse_sequence.find(l1_stop) - len(reverse_sequence)/2]
        l1_read = l1[::-1]
        l1_quality = l1_quality[::-1]
        l1 = ''
        count = 0
        for i in l1_quality:
            if phred[i] > quality_control:
                l1 = l1 + l1_read[count]
            else:
                l1 = l1 + 'X'
            
            count = count + 1





    #print l1
    
    else:
        l1 = 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
        l1_bad = l1_bad + 1

#print reverse_sequence[::-1]
#print l1

    if reverse_sequence.find(l2_stop) - reverse_sequence.find(l2_start) == 28:
        l2_good = l2_good + 1
        l2 = reverse_sequence[reverse_sequence.find(l2_start) + 7:reverse_sequence.find(l2_stop)]
        l2_quality = reverse_sequence[reverse_sequence.find(l2_start) - len(reverse_sequence)/2 + 7:reverse_sequence.find(l2_stop) - len(reverse_sequence)/2]
        l2_read = l2[::-1]
        l2_quality = l2_quality[::-1]
        l2 = ''
        count = 0
        for i in l2_quality:
            if phred[i] > quality_control:
                l2 = l2 + l2_read[count]
            else:
                l2 = l2 + 'X'
            
            count = count + 1
    #print l2
    else:
        #print long_sequence
        l2 = 'XXXXXXXXXXXXXXXXXXXXX'
        l2_bad = l2_bad + 1
#print l2

#print reverse_sequence.find(l3_start), reverse_sequence.find(l3_stop)

    if reverse_sequence.find(l3_stop) - reverse_sequence.find(l3_start) == 25 :
        l3 = reverse_sequence[reverse_sequence.find(l3_start) + 7:reverse_sequence.find(l3_stop)]
        l3_good = l3_good + 1
        l3_quality = reverse_sequence[reverse_sequence.find(l3_start) - len(reverse_sequence)/2 + 7:reverse_sequence.find(l3_stop) - len(reverse_sequence)/2]
        l3_read = l3[::-1]
        l3_quality = l3_quality[::-1]
        l3 = ''
        count = 0
        for i in l3_quality:
            if phred[i] > quality_control:
                l3 = l3 + l3_read[count]
            else:
                l3 = l3 + 'X'
            
            count = count + 1

    #print l3
    else:
        l3 = 'XXXXXXXXXXXXXXXXXX'
        l3_bad = l3_bad + 1
#print l3
#print l3

    cdrs = h1 + h2 + h3 + l1 + l2 + l3
    
#print cdrs
    
    
    if cdrs.find('X') == -1:
        priter.writerow((cdrs,))
        writer.writerow((cdrs, ''))
        perfect_gene = perfect_gene + 1
    else:
        non_perfect_gene = non_perfect_gene + 1
        writer.writerow((cdrs, ''))


#good = no ambigous bases, bad = ambigous bases, perfect gene = no ambiguos bases in the CDRs, not perfect gene = at least one ambiguos base somewhere in a CDR

print 'quality_control', quality_control
print 'h1', h1_good, h1_bad
print 'h2', h2_good, h2_bad
print 'h3', h3_good, h3_bad
print 'l1', l1_good, l1_bad
print 'l2', l2_good, l2_bad
print 'l3', l3_good, l3_bad
print 'perfect genes', perfect_gene, 'not perfect genes', non_perfect_gene

###Part 3: Translates CDR sequences


#codon dictionary

codontable = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

######translates "perfect CDRs"

translation = ''
ofile  = open('aa_correct.csv', 'w')
writer = csv.writer(ofile)
count = 0

with open('cdrs_correct.csv', 'rb') as f:
    reader = csv.reader(f)
    for row in reader:
        translation = ''
        sequence = str(row)
        sequence = sequence[2:152]
        count = count + 1
        #print count
        
        if len(sequence) == 150:
            #print 'shit'
            #print sequence
            #print str(row)
        
        
            for n in range(0, 150, 3):
                if sequence[n:n + 3] in codontable:
                    translation += codontable[sequence[n:n + 3]]
                else:
                    translation += 'X'

            writer.writerow( (translation, '') )

        else:
            continue


######translates all CDRs

translation = ''
ofile  = open('aa_all.csv', 'w')
writer = csv.writer(ofile)
count = 0

with open('cdrs_all.csv', 'rb') as f:
    reader = csv.reader(f)
    for row in reader:
        translation = ''
        sequence = str(row)
        sequence = sequence[2:152]
        count = count + 1
        #print count
        
        
        if len(sequence) == 150:
            #print 'shit'
            #print sequence
            #print str(row)
            #print sequence, 'right'
            
            
            for n in range(0, 150, 3):
                if sequence[n:n + 3] in codontable:
                    translation += codontable[sequence[n:n + 3]]
                else:
                    translation += 'X'
    
            writer.writerow( (translation, '') )
        
        else:
            #print sequence, 'wrong'
            continue


###Part 4: Calculates percentage of mutated CDRs and number of mutations per gene using the "perfect genes"###

h1 = 0
h2 = 0
h3 = 0
l1 = 0
l2 = 0
l3 = 0

mutations = [0, 0, 0, 0, 0, 0, 0]

count = 0

with open('aa_correct.csv', 'rb') as f:
    reader = csv.reader(f)
    for row in reader:
        mut = 0
        if 'SSYAMS' not in row[0]:
            h1 = h1 + 1
            mut = mut + 1
        
        if 'AISGSGGSTY' not in row[0]:
            h2 = h2 + 1
            mut = mut +1
        
        if 'APYSGGYSSF' not in row[0]:
            h3 = h3 + 1
            mut = mut + 1
        
        if 'RASQSISSYLN' not in row[0]:
            l1 = l1 + 1
            mut = mut + 1
        
        if 'AASSLQS' not in row[0]:
            l2 = l2 + 1
            mut = mut + 1
        
        if 'GHYWLS' not in row[0]:
            l3 = l3 + 1
            mut = mut + 1
        
        
        mutations[mut] = mutations[mut] + 1
        
        count = count + 1

count = float(count)
h1 = float(h1)
h2 = float(h2)
h3 = float(h3)
l1 = float(l1)
l2 = float(l2)
l3 = float(l3)

print ''
print 'mutation in h1', h1/count*100, '%'
print 'mutation in h2', h2/count*100, '%'
print 'mutation in h3', h3/count*100, '%'
print 'mutation in l1', l1/count*100, '%'
print 'mutation in l2', l2/count*100, '%'
print 'mutation in l3', l3/count*100, '%'

print ''
print 'zero mutated', mutations[0]/count*100, '%'
print 'one mutated', mutations[1]/count*100, '%'
print 'two mutated', mutations[2]/count*100, '%'
print 'three mutated', mutations[3]/count*100, '%'
print 'four mutated', mutations[4]/count*100, '%'
print 'five mutated', mutations[5]/count*100, '%'
print 'six mutated', mutations[6]/count*100, '%'

###Part 5: Calculates the aminacid variation at each positions using###

from collections import Counter

supercount = 0

with open('aa_all.csv', 'rb') as f:                        # Calculates "row_length", the length of the sequences
    reader = csv.reader(f)
    for row in reader:
        row = str(row[0])
        row_length = len(row)
        break


positions = list(range(0, row_length))
lists = [[] for i in xrange(len(positions))]        # Makes a list of lists, one for each position in the sequence

template = 'SSYAMSAISGSGGSTYAPYSGGYSSFRASQSISSYLNAASSLQSGHYWLS'

listcount = -1
for i in positions:                                 # Adds all aminoacids for each position to its list.
    with open('aa_all.csv', 'rb') as f:
        reader = csv.reader(f)
        listcount += 1
        print listcount
        for row in reader:
            row = str(row[0])
            if len(row) == 50:
                if row[i] == 'X':
                    continue
                else:
                    lists[listcount].append(row[i])
            else:
                continue


with open('aa_all.csv', 'rb') as f:                        # Counts number of alignments
    reader = csv.reader(f)
    rowcount = 0
    for row in reader:
        rowcount += 1
rowcount = float(rowcount)


ofile = open('output_positioncounter.csv', 'w')
writer = csv.writer(ofile)
writer.writerow(('# of alignments in set',))
writer.writerow((rowcount, ))
writer.writerow(('Quality control',))
writer.writerow((quality_control, ))
writer.writerow(( ))


#writer.writerow( ('Position', 'Aminoacid', '#', '%' ))


variations = list()                                 # Determines if a position has variation, and if so calculates the variation
position = []
aminoacid = []
occurances = []
procentage = []
n = 0

aminoacids = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'C', 'G', 'P', '*' ]

for i in lists:
    counter = Counter(i)
    countersum = float(sum(counter.itervalues()))
    
    for i in aminoacids:
        #print i
        if i in counter:
            x = counter[i]/countersum * 100
            aminoacid.append(i)
            occurances.append(counter[i])
            procentage.append(x)
            position.append(positions[n]+1)
        else:
            aminoacid.append(i)
            occurances.append(0)
            procentage.append(0)
            position.append(positions[n]+1)
    n += 1

positions_correct = []

count = 0

while count < len(positions):
    count = count +1
    positions_correct.append(count)
    positions_correct.append(' ')
    positions_correct.append(' ')


template = ("  ".join(template))



writer.writerow(( ))
writer.writerow(( 'Variation histogram', ))
writer.writerow(( template))
writer.writerow(( positions_correct ))


count = 0
listrow = []
last = 0
nextaminoacid = aminoacid
nextprocentage = procentage
nextoccurance = occurances
nextposition = position

#most_varied_position = max(set(position), key=position.count)
#count_most_varied_position= position.count(most_varied_position)


loopcount = [0] * 25


for p in loopcount:
    aminoacid = nextaminoacid
    procentage = nextprocentage
    occurances = nextoccurance
    position = nextposition
    
    nextaminoacid = []
    nextprocentage = []
    nextoccurance = []
    nextposition = []
    aminorow = []
    last = 0
    
    
    while len(position) > 0:
        i = position[0]
        
        listofzeros = []
        
        space = (i-last-1) * 3
        listofzeros = [0] * space
        for n in listofzeros:
            listrow.append(n)
            aminorow.append('')
        listrow.append(i)
        aminorow.append(aminoacid[0])
        aminorow.append(procentage[0])
        aminorow.append(occurances[0])
        
        
        listofzeros = [0] * 2
        for n in listofzeros:
            listrow.append(n)

        position.pop(0)
        aminoacid.pop(0)
        occurances.pop(0)
        procentage.pop(0)
        
        while i in position:
            
            nextposition.append(position[0])
            position.pop(0)
            
            nextaminoacid.append(aminoacid[0])
            aminoacid.pop(0)
            
            nextoccurance.append(occurances[0])
            occurances.pop(0)
            
            nextprocentage.append(procentage[0])
            procentage.pop(0)
        
        last = i

    writer.writerow((aminorow))

print ''
print 'Done'
print ''
print 'Results are found in output_positioncounter.csv'
print ''



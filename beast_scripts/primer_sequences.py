import sys
sys.path.append("/home/dyap/INSTALL/argparse-1.2.1")
import argparse
import csv
import string

def read_sequences(fasta):
	id = None
	sequences = []
	for line in fasta:
		line = line.rstrip()
		if len(line) == 0:
			continue
		if line[0] == '>':
			if id is not None:
				yield (id, ''.join(sequences))
			id = line[1:].split()[0]
			sequences = []
		else:
			sequences.append(line)
	if id is not None:
		yield id, ''.join(sequences)
		
		
def reverse_complement(sequence):
	return sequence[::-1].translate(string.maketrans('ACTGactg','TGACtgac'))


def create_breakpoint_sequence(reference_sequences, entry, extend):
	breakpoint_sequences = ['', '']
	expected_strands = ('+', '-')
	for side in (0, 1):
		chromosome = entry['chromosome_' + str(side+1)]
		if entry['strand_' + str(side+1)] == '+':
			end = int(entry['break_' + str(side+1)])
			start = end - extend + 1
		else:
			start = int(entry['break_' + str(side+1)])
			end = start + extend - 1
		breakpoint_sequences[side] = reference_sequences[chromosome][start-1:end]
		if entry['strand_' + str(side+1)] != expected_strands[side]:
			breakpoint_sequences[side] = reverse_complement(breakpoint_sequences[side])
	return breakpoint_sequences[0] + '[]' + breakpoint_sequences[1]


def create_breakend_sequence(reference_sequences, entry, extend):
	breakend_sequences = ['', '']
	for side in (0, 1):
		chromosome = entry['chromosome_' + str(side+1)]
		if entry['strand_' + str(side+1)] == '+':
			breakend_left = int(entry['break_' + str(side+1)])
		else:
			breakend_left = int(entry['break_' + str(side+1)]) - 1
		breakend_sequences[side] = reference_sequences[chromosome][breakend_left-extend:breakend_left] + '[]' + reference_sequences[chromosome][breakend_left:breakend_left+extend]
	return tuple(breakend_sequences)


argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
argparser.add_argument('genome', help='genome fasta filename')
argparser.add_argument('length', help='extend sequence length', type=int)
args = argparser.parse_args()

reference_sequences = dict(read_sequences(open(args.genome, 'r')))


reader = csv.reader(sys.stdin, delimiter='\t')
header = next(reader)

print '\t'.join(['cluster_id', 'breakpoint_sequence', 'breakend1_sequence', 'breakend2_sequence'])
for row in reader:
	entry = dict(zip(header, row))
	if int(entry['num_split']) == 0:
		entry['break_1'] = (entry['start_1'], entry['end_1'])[entry['strand_1'] == '+']
		entry['break_2'] = (entry['start_2'], entry['end_2'])[entry['strand_2'] == '+']
	breakpoint_sequence = create_breakpoint_sequence(reference_sequences, entry, args.length)
	breakend1_sequence, breakend2_sequence = create_breakend_sequence(reference_sequences, entry, args.length)
	print '\t'.join([entry['cluster_id'], breakpoint_sequence, breakend1_sequence, breakend2_sequence])

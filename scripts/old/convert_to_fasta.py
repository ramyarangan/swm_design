import sys
import random

full_seq = "gcggauuuagcucaguugggagagcgccagacugaagaucuggagguccuguguucgauccacagaauucgcacca"

fasta_filename = sys.argv[1]
control_fasta_filename = sys.argv[2]
score_filename = sys.argv[3]

def convert_to_canonical(seq):
	new_seq = ""
	i = 0
	while i < len(seq):
		if seq[i] != 'X':
			new_seq += seq[i]
			i += 1
		else:
			new_seq += seq[i + 4].lower()
			i += 6
	return new_seq

def get_random_base():
	rnd_idx = int(random.random() * 4)
	if rnd_idx == 0:
		return 'c'
	if rnd_idx == 1:
		return 'g'
	if rnd_idx == 2:
		return 'a'
	return 'u'

def get_random_seq(len):
	char_arr = [get_random_base() for ii in range(len)]
	return ''.join(char_arr)

f = open(score_filename)
score_lines = f.readlines()
f.close()

control_f = open(control_fasta_filename, 'w')
f = open(fasta_filename, 'w')

for ii in range(len(score_lines)):
	score_line = score_lines[ii]
	score = float(score_line.split()[0])
	insert_seq = convert_to_canonical(score_line.split()[1])
	new_seq = full_seq[0:12] + insert_seq[0:10] + full_seq[22:45] + insert_seq[10:] + full_seq[65:]
	new_control = full_seq[0:18] + get_random_seq(2) + full_seq[20:53] + get_random_seq(7) + full_seq[60:]
	name_str = ">" + str(ii) + " " + str(score)
	f.write("%s\n" % name_str)
	f.write("%s\n" % new_seq)
	control_f.write("%s\n" % name_str)
	control_f.write("%s\n" % new_control)

f.close()
control_f.close()
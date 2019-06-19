import sys
import random
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd 

idx_filename = sys.argv[1]
score_filename = sys.argv[2]
alignment_filename = sys.argv[3]

## Stores the matching between the indices in a designed SWM sequence,
## in a full sequence for the benchmark PDB, and in an RFAM alignment
class IdxMatchings:
	def __init__(self, idx_lines): 

		self.idx_matchings = []
		self.swm_idxs = []
		self.full_idxs = []
		self.align_idxs = []

		for idx_line in idx_lines[1:]:
			idx_items = [int(x) for x in idx_line.split(',')]
			self.idx_matchings += [idx_items]
			self.swm_idxs += [idx_items[0]]
			self.full_idxs += [idx_items[1]]
			self.align_idxs += [idx_items[2]]

		self.update_full_to_swm()

	def update_full_to_swm(self):
		self.full_to_swm = {}
		for idx_matching in self.idx_matchings:
			self.full_to_swm[idx_matching[1]] = idx_matching[0]

## Data used to compute match of designs to alignment
## Later this might include secondary structure
class MatchParams:
	def __init__(self, idx_filename): 
		f = open(idx_filename)
		idx_lines = f.readlines()
		f.close()

		self.full_seq = idx_lines[0]
		self.idx_matching = IdxMatchings(idx_lines[1:])		

## Set of SWM designed sequences and their corresponding scores, and
## a matched length set of control sequences
class SWMControlSeqsScores:
	def __init__(self, score_filename, match_params):		
		score_file = open(score_filename)
		score_lines = score_file.readlines()
		score_file.close()

		self.scores = self.get_swm_scores(score_lines)
		self.swm_seqs = []
		self.control_seqs = []
		self.fill_swm_control_seqs(score_lines, match_params)

	## Utility function: convert sequence to all canonicals
	def convert_to_canonical(self, seq):
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

	## Utility function: get a random nucleotide
	def get_random_base(self):
		rnd_idx = int(random.random() * 4)
		if rnd_idx == 0:
			return 'c'
		if rnd_idx == 1:
			return 'g'
		if rnd_idx == 2:
			return 'a'
		return 'u'

	## Utility function: get a random set of nucleotides
	def get_random_seq(self, len):
		char_arr = [self.get_random_base() for ii in range(len)]
		return ''.join(char_arr)

	## Get SWM scores from design file
	def get_swm_scores(self, score_lines):
		return [float(score_line.split()[0]) for score_line in score_lines]

	## Place SWM designed sequence in full context
	def insert_swm_design_full_seq(self, match_params, insert_seq):
		full_seq = match_params.full_seq
		idx_matching = match_params.idx_matching
		new_seq = list(full_seq)
		for full_idx, design_idx in idx_matching.full_to_swm.items():
			new_seq[full_idx] = insert_seq[design_idx]
		return ''.join(new_seq)

	## Update data for all SWM and control sequences
	def fill_swm_control_seqs(self, score_lines, match_params):
		for ii in range(len(score_lines)):
			score_line = score_lines[ii]
			insert_seq = self.convert_to_canonical(score_line.split()[1])
			new_swm = self.insert_swm_design_full_seq(match_params, insert_seq)
			# Control sequence has random nucleotides at designed positions
			new_control = self.insert_swm_design_full_seq(match_params, self.get_random_seq(len(insert_seq)))
			self.swm_seqs += [new_swm]
			self.control_seqs += [new_control]

## Stores the alignment lines, processes from file
## Later this should include the consensus alignment / secondary structure
class Alignment:
	def __init__(self, alignment_filename):	
		alignment_file = open(alignment_filename)
		alignment_file_lines = alignment_file.readlines()
		alignment_file.close()

		self.alignment_lines = self.alignment_from_file(alignment_file_lines)

	def alignment_from_file(self, lines):
		alignment_lines = []

		for line in lines:
			if line == "\n" or line[0] == "#" or line[0] == "/":
				continue
			next_line = " ".join(line.split())
			next_line = next_line.split()[1]
			alignment_lines += [next_line]

		return alignment_lines

## Stores a PWM and computes match to it
class PWM:
	nt_idxs = {"A": 0, "U": 1, "C": 2, "G": 3}
	
	def __init__(self, alignment, match_params):
		align_idxs = match_params.idx_matching.align_idxs
		self.pwm = self.get_pwm(alignment, align_idxs)
		self.size = self.pwm.shape[1]
		self.pwm_df = pd.DataFrame(data=self.pwm.T)
		self.pwm_df.columns = ["A", "U", "C", "G"]

	def get_pwm(self, alignment, align_idxs):
		# Add in pseudocounts to prevent -inf in PWM
		counts = np.array([1] * len(align_idxs) * 4) 
		counts.resize((4, len(align_idxs)))

		for alignment_line in alignment.alignment_lines:
			for ii, idx in enumerate(align_idxs):
				new_char = alignment_line[idx]
				if new_char == "-":
					continue
				new_char = new_char.upper()
				counts[self.nt_idxs[new_char], ii] += 1

		counts = counts.T
		counts = counts/counts.sum(axis=1, keepdims=True)
		counts = counts.T
		counts = np.log(counts/0.25)

		return counts

	def score_pwm_match(self, seq, idxs):
		if (len(idxs) != self.size):
			raise TypeError
		align_score = 0
		for jj, idx in enumerate(idxs):
			new_char = seq[idx].upper()
			align_score += self.pwm[self.nt_idxs[new_char], jj]
		return align_score

## Computes and stores scores for SWM and control sequences, scored against a PWM
## Includes plotting utilities for displaying SWM / control scores
class SWMControlPWMScore:
	def __init__(self, swm_control_seqs_scores, pwm, match_params):
		self.idxs = match_params.idx_matching.full_idxs
		self.pwm = pwm
		self.swm_pwm_scores = self.get_pwm_scores(swm_control_seqs_scores.swm_seqs)
		self.control_pwm_scores = self.get_pwm_scores(swm_control_seqs_scores.control_seqs)
		self.scores = np.array(swm_control_seqs_scores.scores)

	## Plot histogram of SWM and Control sequence matches to PWM
	def plot_hist_scores(self, score_cutoff=0.1):
		swm_scores = self.filter_pwm_scores(self.swm_pwm_scores, score_cutoff)
		control_scores = self.filter_pwm_scores(self.control_pwm_scores, score_cutoff)

		plt.hist([swm_scores, control_scores], color=['blue', 'forestgreen'], \
			alpha=0.7, label=['swm', 'control'])#, rwidth=0.85)
		plt.legend(loc='upper left')
		plt.axvline(np.mean(swm_scores), color='blue', linestyle='dashed')
		plt.axvline(np.mean(control_scores), color='forestgreen', linestyle='dashed')
		plt.show()

	## Plot median score difference between SWM and control sequences as they vary by
	## SWM designd structure score cutoff
	def plot_median_score_diff(self):
		score_cutoff_range = np.arange(1, 0, -0.01)
		diff_medians = []
		for score_cutoff in score_cutoff_range:	
			align_scores = self.filter_pwm_scores(self.swm_pwm_scores, score_cutoff)
			control_scores = self.filter_pwm_scores(self.control_pwm_scores, score_cutoff)
			diff_medians += [np.median(align_scores) - np.median(control_scores)]
		plt.scatter(score_cutoff_range, diff_medians, color='black')
		plt.xlabel("% Top-scoring structures")
		plt.ylabel("Median SWM - control PWM score")
		plt.show()

	## Plot median score of SWM sequences as they vary by SWM designed 
	## structure score cutoff
	def plot_median_score(self):
		score_cutoff_range = np.arange(1, 0, -0.01)
		medians = []
		for score_cutoff in score_cutoff_range:	
			align_scores = self.filter_pwm_scores(self.swm_pwm_scores, score_cutoff)
			medians += [np.median(align_scores)]
		plt.scatter(score_cutoff_range, medians, color='black')
		plt.xlabel("% Top-scoring structures")
		plt.ylabel("Median SWM PWM score")
		plt.show()

	## Get PWM scores for a set of sequences
	def get_pwm_scores(self, seqs):
		align_scores = []
		for seq in seqs:
			align_scores += [self.pwm.score_pwm_match(seq, self.idxs)]
		return np.array(align_scores)

	## Filter PWM scores by a score cutoff
	def filter_pwm_scores(self, pwm_scores, score_cutoff):
		scores = self.scores
		idxs = np.argsort(scores)
		scores = scores[idxs]
		pwm_scores = pwm_scores[idxs]
		
		num_items = int(score_cutoff * len(scores))
		scores = scores[0:num_items]
		pwm_scores = pwm_scores[0:num_items]

		return pwm_scores


match_params = MatchParams(idx_filename)
swm_control_seqs_scores = SWMControlSeqsScores(score_filename, match_params)
alignment = Alignment(alignment_filename)
pwm = PWM(alignment, match_params)
print(pwm.pwm_df)
swm_control_pwm_score = SWMControlPWMScore(swm_control_seqs_scores, pwm, match_params)
swm_control_pwm_score.plot_hist_scores()
swm_control_pwm_score.plot_median_score_diff()

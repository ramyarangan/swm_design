import sys
import numpy as np 
import matplotlib.pyplot as plt
import logomaker as lm 
import pandas as pd 

alignment_filename = sys.argv[1]
swm_filename = sys.argv[2]
control_filename = sys.argv[3]

alignment_file = open(alignment_filename)
alignment_lines = alignment_file.readlines()
alignment_file.close()

swm_file = open(swm_filename)
swm_lines = swm_file.readlines()
swm_file.close()

control_file = open(control_filename)
control_lines = control_file.readlines()
control_file.close()


nt_idxs = {"A": 0, "U": 1, "C": 2, "G": 3}

def get_pwm(alignment_lines, align_idxs):
	# Add in pseudocounts to prevent -inf in PWM
	counts = np.array([1] * len(align_idxs) * 4) 
	counts.resize((4, len(align_idxs)))
	for alignment_line in alignment_lines:
		if alignment_line == "\n" or alignment_line[0] == "#" or \
			alignment_line[0] == "/":
			continue
		next_line = " ".join(alignment_line.split())
		next_line = next_line.split()[1]
		for ii, idx in enumerate(align_idxs):
			new_char = next_line[idx]
			if new_char == "-":
				continue
			new_char = new_char.upper()
			counts[nt_idxs[new_char], ii] += 1

	counts = counts.T
	counts = counts/counts.sum(axis=1, keepdims=True)
	counts = counts.T
	counts = np.log(counts/0.25)
	return counts

def get_scores(swm_lines, pwm, swm_idxs, score_cutoff=0.1):
	scores = []
	align_scores = []
	for ii in range(int(len(swm_lines)/2)):
		scores += [float(swm_lines[ii * 2].split()[1])]
		new_line = swm_lines[ii * 2 + 1]
		align_score = 0
		for jj, idx in enumerate(swm_idxs):
			new_char = new_line[idx].upper()
			align_score += pwm[nt_idxs[new_char], jj]
		align_scores += [align_score]
	
	scores = np.array(scores)
	align_scores = np.array(align_scores)
	idxs = np.argsort(scores)
	scores = scores[idxs]
	align_scores = align_scores[idxs]
	
	num_items = int(score_cutoff * len(scores))
	scores = scores[0:num_items]
	align_scores = align_scores[0:num_items]

	return [scores, align_scores]

def plot_hist_scores(swm_lines, control_lines, pwm, swm_idxs, score_cutoff=0.1):
	[scores, align_scores] = get_scores(swm_lines, pwm, swm_idxs, score_cutoff=score_cutoff)
	[scores, control_align_scores] = get_scores(control_lines, pwm, swm_idxs, score_cutoff=score_cutoff)

	plt.hist([np.array(align_scores), np.array(control_align_scores)], \
		color=['blue', 'forestgreen'], alpha=0.7, label=['swm', 'control'])#, rwidth=0.85)
	plt.legend(loc='upper left')
	plt.axvline(np.mean(align_scores), color='blue', linestyle='dashed')
	plt.axvline(np.mean(control_align_scores), color='forestgreen', linestyle='dashed')
	plt.show()

def plot_median_score_diff(swm_lines, control_lines, pwm, swm_idxs):
	score_cutoff_range = np.arange(1, 0, -0.01)
	diff_medians = []
	for score_cutoff in score_cutoff_range:	
		[_, align_scores] = get_scores(swm_lines, pwm, swm_idxs, score_cutoff=score_cutoff)
		[_, control_align_scores] = get_scores(control_lines, pwm, swm_idxs, score_cutoff=score_cutoff)
		diff_medians += [np.median(align_scores) - np.median(control_align_scores)]
	plt.scatter(score_cutoff_range, diff_medians, color='black')
	plt.xlabel("% Top-scoring structures")
	plt.ylabel("Median SWM - control PWM score")
	plt.show()

def plot_median_score(swm_lines, pwm, swm_idxs):
	score_cutoff_range = np.arange(1, 0, -0.01)
	medians = []
	for score_cutoff in score_cutoff_range:	
		[_, align_scores] = get_scores(swm_lines, pwm, swm_idxs, score_cutoff=score_cutoff)
		medians += [np.median(align_scores)]
	plt.scatter(score_cutoff_range, medians, color='black')
	plt.xlabel("% Top-scoring structures")
	plt.ylabel("Median SWM PWM score")
	plt.show()

# Later automate finding these indices
## HOW TO GET THIS RIGHT??
swm_idxs = [18, 19, 53, 54, 55, 56, 57, 58, 59]
align_idxs = [21, 22, 91, 93, 94, 95, 97, 101, 102]

pwm = get_pwm(alignment_lines, align_idxs)
df = pd.DataFrame(data=pwm.T)
df.columns = ["A", "U", "C", "G"]
print(df)

plot_hist_scores(swm_lines, control_lines, pwm, swm_idxs)
# plot_median_score_diff(swm_lines, control_lines, pwm, swm_idxs)
plot_median_score(swm_lines, pwm, swm_idxs)


import sys
from matplotlib import pyplot as plt

stats_filename = sys.argv[1]
score_filename = sys.argv[2]
control_filename = sys.argv[3]

f = open(stats_filename)
stats_lines = f.readlines()
f.close()

f = open(score_filename)
score_lines = f.readlines()
f.close()

f = open(control_filename)
control_lines = f.readlines()
f.close()

stat_idx = 6

pps = []
for stats in stats_lines:
	if stats[0] != '#':
		one_space = ' '.join(stats.split())
		pps += [float(one_space.split()[stat_idx])]

controls = []
for stats in control_lines:
	if stats[0] != '#':
		one_space = ' '.join(stats.split())
		controls += [float(one_space.split()[stat_idx])]

scores = []
for score_line in score_lines:
	scores += [float(score_line.split()[0])]

plt.scatter(pps, scores)
plt.show()

plt.hist(pps, alpha=0.5, label='designs')
plt.hist(controls, alpha=0.5, label='controls')
plt.legend(loc='upper right')
plt.show()
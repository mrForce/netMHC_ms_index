import sys
assert(len(sys.argv) == 2)
fpr = float(sys.argv[1])
"""
Each line should have three columns, seperated by whitespace: peptide  decoy/target  xcorr
"""
lines = sys.stdin.readlines()
data = []
from collections import defaultdict
xcorr_clusters = defaultdict(list)
for line in lines:
    sections = line.split()
    peptide = sections[0]
    decoy = sections[1]
    xcorr = float(sections[2])
    xcorr_clusters[xcorr].append((peptide, decoy))




num_target = 0
num_decoy = 0
targets = []
passing_index = -1
for xcorr_threshold  in sorted(xcorr_clusters.keys(), reverse=True):
    for peptide, decoy in xcorr_clusters[xcorr_threshold]:        
        if 'decoy' in decoy:
            num_decoy += 1
        elif 'target' in decoy:
            num_target += 1
            targets.append(peptide)
    false_positives = num_decoy*2
    fp_rate = 1.0*false_positives/(num_target + num_decoy)
    if fp_rate <= fpr:
        passing_index = len(targets)
        
if passing_index == -1:
    print('NO SOLUTION')
else:
    for peptide in targets[0:passing_index]:
        print(peptide)

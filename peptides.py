import subprocess
import os
import re

completed_netmhc = subprocess.run(["netmhc", "-listMHC"], stdout=subprocess.PIPE)
species_list = {'human':'HLA', 'mouse':'H-'}
i = 1
starter_selector = {}
print('Species list:')
for k,v in species_list.items():
    print('{0}: '.format(i) + k)
    starter_selector[i] = k
    i += 1

j = int(input('Select a species: '))


starter = species_list[starter_selector[j]]
species = starter_selector[j]
netmhc_output = list(filter(lambda x:  x.startswith(starter), completed_netmhc.stdout.decode().split('\n')))
output_table = {}
i = 1
selector = {}
print('MHC list:')
for j in range(0, len(netmhc_output)):
    mhc = netmhc_output[j]
    print('{0}: '.format(i) + mhc)
    selector[i] = mhc
    i += 1
k = int(input('Select an MHC variant: '))
print('You selected: ' + selector[k])
mhc_type = selector[k]
lengths = [int(x) for x in input('Type lengths seperated by commas: ').split(',')]
print('lengths:')
print(lengths)
keepGoing = True
fasta_location = ''
while keepGoing:
    fasta_location = input('Enter FASTA location: ')
    if os.path.isfile(fasta_location):
        keepGoing = False
    else:
        print('That is not a valid file. Try again')

cutoff_type= None
keepGoing = True
while keepGoing:
    cutoff_type = input('Type affinity for IC50 cutoff, rank for rank cutoff: ')
    if cutoff_type == 'affinity' or cutoff_type == 'rank':
        keepGoing = False
    else:
        print('Not a valid cutoff type')

cutoff_value = 0
if cutoff_type == 'affinity':
    keepGoing = True

    while keepGoing:
        try:
            cutoff_value = float(input('Enter the affinity cutoff: '))
        except ValueError:
            print('Not a valid affinity cutoff')
        else:
            keepGoing = False
elif cutoff_type == 'rank':
    keepGoing = True

    while keepGoing:
        try:
            cutoff_value = float(input('Enter the rank cutoff: '))
        except ValueError:
            print('Not a valid rank cutoff')
        else:
            keepGoing = False

            

output_prefix = species + '_' + mhc_type + '_' + '_'.join([str(x) for x in lengths]) + '_' + os.path.split(fasta_location)[1]
print('output prefix: ' + output_prefix)

r = re.compile("-{4,}\n(?P<header>.+)\n-{4,}\n(?P<results>(?:(?:.+\n))+)-{4,}")
rows = []
netmhc_output_file = open(output_prefix + '.tsv', 'w')
print('going to run netMHC')
completed_netmhc = subprocess.run(["netmhc", "-a", mhc_type, "-l", ','.join([str(x) for x in lengths]), '-f ' + fasta_location], stdout = netmhc_output_file)
print('netmhc output  can be found here: ' + output_prefix + '.tsv')
netmhc_output_file.close()
netmhc_output = open(output_prefix + '.tsv', 'r').read()

print('\n'.join(netmhc_output.split('\n')[0:20]))

for header, results in r.findall(netmhc_output):
    hla_column = 1
    peptide_column = 2
    affinity_column = 11
    rank_column = 12
    for line in results.split('\n'):
        if len(line) > 5:
            line_sections = re.split('\s{2,}', line)[1::]
            row = [line_sections[peptide_column], line_sections[hla_column], re.search('\d*\.\d*', line_sections[affinity_column]).group(0), re.search('\d*\.\d*', line_sections[rank_column]).group(0)]
            if cutoff_type == 'affinity' and float(row[2]) <= cutoff_value:
                rows.append(row)
            elif cutoff_type == 'rank' and float(row[3]) <= cutoff_value:
                rows.append(row)

row_file = open(output_prefix + '_rows.txt', 'w')
peptide_file = open(output_prefix + '_peptides.fasta', 'w')
i = 1
for row in rows:
    row_file.write(', '.join(row) + '\n')
    peptide_file.write('> {0}\n'.format(i))
    peptide_file.write(row[0] + '\n')
    i += 1
row_file.close()
peptide_file.close()

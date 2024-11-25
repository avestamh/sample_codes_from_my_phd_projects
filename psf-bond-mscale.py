import re

# User specification
inp_name = "mscale-i27-1x2ssra-1do2.psf"
out_name = "mscale-1do2-i27-1x2ssra_bonding.psf"
## match coarse-grained ClpY atoms: may need try & error several times
pattern_clpy = r'\s+\d+\sSEG.\s\d+\s+\w+\s{2}CA\s{3}CH\dE\s{4}0\.1*'

# Data containers
prev_bonds = []
old_bonds  = []
post_bonds = []
clpy_atoms = []

# Get bonds information of each entry
def process_bond_entry(inp_line):
    each_entry = inp_line.split()
    num_of_bonds = len(each_entry)/2
    return num_of_bonds, [(each_entry[i*2], each_entry[i*2+1]) for i in range(num_of_bonds)]

# Prev- and Post-bonds entries are trivial. BOND terms are recorded specifically.
with open(inp_name, 'r') as inp_file:
    nbonds = 0
    bound_counter = 0
    is_before_nbond = True
    is_after_nbond  = False
    for each_line in inp_file:
        if re.match(pattern_clpy, each_line):
            each_entry = each_line.split()
            clpy_atoms.append((int(each_entry[0]), each_entry[1]))
        if '!NBOND' in each_line:
            nbonds = int(each_line.split()[0])
            is_before_nbond = False
            continue
        if nbonds:
            nbond, each_bond = process_bond_entry(each_line)
            bound_counter += nbond
            nbonds, is_after_nbond = [(nbonds, False), (0, True)][bound_counter==nbonds]
            old_bonds.extend(each_bond)
        if is_before_nbond:
            prev_bonds.append(each_line)
        if is_after_nbond:
            post_bonds.append(each_line)

# Merge old bonds and new bonds derived from ClpY atom index
new_bonds = [(prev_atom[0], next_atom[0]) for prev_atom, next_atom in zip(clpy_atoms[:-1], clpy_atoms[1:]) if prev_atom[1]==next_atom[1]] + old_bonds

# Write new psf
with open(out_name, 'w') as out_file:
    out_file.write(''.join(prev_bonds))
    out_file.write('%8d!NBOND\n'%len(new_bonds))
    for i, nb in enumerate(new_bonds, 1):
        out_file.write("%8s%8s"%(nb[0], nb[1]))
        out_file.write("\n" if 0==i%4 else "")
    out_file.write('\n')
    out_file.write(''.join(post_bonds))

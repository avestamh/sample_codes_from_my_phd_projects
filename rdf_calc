

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import Series, DataFrame
import seaborn as sb
from matplotlib import rcParams
import re
import bio.physics
bio.physics.sucks(14)
rcParams['figure.figsize']=5,4
# sb.set_style('whitegrid')

rdf_eef1="addres/file1.dat"
rdf_go="address/file2.dat"
entry_pattern = '\s+[0-9]+\.[0-9]+\s+[0-9]\.[0-9]+'

def gofr_p(inp_name):

    gofr = []
    # flag, color = ('EEF1', 'green') if 'eef1' in inp_name else ('CG', 'red')
    flag, color = ('EEF1', 'green')  if inp_name == rdf_eef1 else ('CG', 'red')
    print flag
    with open(inp_name, "r") as inp_file:
        for each_line in inp_file:
            if re.match(entry_pattern, each_line):
                each_entry = each_line.split()
                # each_entry = map(float, each_entry)
                each_entry = [float(i) for i in each_entry]
                gofr.append(each_entry)

    rdf = pd.DataFrame(gofr)
    rdf.columns = ['bin_mid', 'norm_gr', 'tot_num_pairs']
    x = rdf['bin_mid']
    y = rdf['norm_gr']
    plt.plot(x, y, label=flag, c=color) 
    plt.xlabel('r ($\AA$)', fontsize = 16)
    plt.ylabel("g(r)", fontsize = 16)
    plt.legend()

gofr_p(rdf_eef1)
gofr_p(rdf_go)
img_name = "file.png"
plt.savefig(img_name, dpi=300, bbox_inches='tight')
plt.show()


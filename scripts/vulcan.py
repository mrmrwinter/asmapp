#!conda run -n vulcan python

import glob
import os
import pandas as pd


experiment_name = "VW4_HiFiasm.15.p_ctg.purge.SLR/"
out_dir = "outputs/pairs/vulcan/"
csv_path = "../" + experiment_name + "reports/blast/initial_blast.onlyPairs.tsv"

os.system("mkdir -p ../" + experiment_name + out_dir)

pairs = pd.read_csv(csv_path, sep='\t')

for index, value in pairs.iterrows():
    scaf_one = int(value[1])
    scaf_two = int(value[2])
    q = "../" + experiment_name + "tmp_initial/" + str(scaf_one) + ".fasta"
    h = "../" + experiment_name + "tmp_initial/" + str(scaf_two) + ".fasta"
    
    vulcan_cmd = "vulcan -i " + h + " -r " + q + " -o ../" + experiment_name + out_dir + str(index) + " -w ../" + experiment_name + out_dir + "vulcan_tmp -t 10 -p 0"
    
    os.system(vulcan_cmd)


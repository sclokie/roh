
import os
import sys
import subprocess
import more_itertools
from more_itertools import run_length
import pandas as pd
import csv
pd.set_option('display.expand_frame_repr', False)
from ruffus import *
import itertools
import collections
import pprint
import shutil
##
## python3.9 roh-bins-genome.py working_output_path analysis_name vcf1 vcf2 vcf3
##

# Set to your available cores
processes = 16

# open a temp file created by shiny that has the command arguments. This is becasue nq does not like more than one arg
parsed_args = []
with open('roh-bins-command', 'r+') as inputfile:
    for line in inputfile:
        parsed_args = line.strip().split(" ")
print(parsed_args)
dest_path = parsed_args[0]
output_name = parsed_args[1]
number_of_vcfs = len(parsed_args)
vcfs_to_process = []
for i in range(2,number_of_vcfs):
     vcfs_to_process.append(f"{dest_path}/{parsed_args[i]}")

# This script will:
# 1 - create genome wide bins (bedtools)
# 2 - sort vcf
# 3 - map vcf to windows, keeping vcf details. Separate each variant with th word 'SPLIT'
# 4 - merge runs of het and homs together, keeping GT details
# 5 - create a bed file of hets and homs, with summary details
# 6 - Create a karyotype plot of only homozygous runs
# 7 - Create a jaccard statistic for each homozygous run compared to each other.
# 8 - plot the jaccard stats as a heatmap.


# 1
bedtools = '<set-path>/bedtools'
genome = '<set-path>/hg19.chr.sizes.NoM.NoBlacklist.bed'
regions_to_filter = '<set-path>/roh-regions-to-filter.bed'
window_size = 250000000 # 250000000 gets all data on one chrom
window_bed = f"{genome.split('/')[-1]}-windows-{window_size}.bed"
os.system(f"{bedtools} makewindows -g {genome} -i srcwinnum -w {window_size} | {bedtools} sort | uniq > {dest_path}/{window_bed}")

# R scripts needed for plotting. Make sure they are available in cwd
# plot-ROH-karyotypes.R
# plot-ROH-jaccard_stats.R

print("STARTED")

blacklist = {}
with open(regions_to_filter,'r+') as blacklistfile:
    for line in blacklistfile:
        parts = line.strip().split("\t")
        blacklist[parts[0]] = [parts[1],parts[2]]

#2 sort vcfs
for i in vcfs_to_process:
    print("i:",i)
    line = subprocess.check_output(['tail', '-1', i])
    line = str(line.strip()[:3])

    if 'chr' in line:
        chr_status = 'alphanumeric_chr'
        os.system(f"{bedtools} sort -header -i {i} > {i[:-4]}.sorted.vcf")
    else:
        chr_status = 'numeric_chr'
        with open(i,'r+') as inputfile:
            with open(f'{i[:-4]}.modified.vcf', 'w+') as outputfile:
                for line in inputfile:
                    if line.startswith('#'):
                        outputfile.write(line)
                    else:
                        modified_line = 'chr' + line
                        outputfile.write(modified_line)
        os.system(f"{bedtools} sort -header -i {i[:-4]}.modified.vcf > {i[:-4]}.sorted.vcf")

#3 map variants to windows
for i in vcfs_to_process:
    os.system(f"{bedtools} map -c 2,10 -o collapse,collapse -delim \"SPLIT\" -a {window_bed} -b {i[:-4]}.sorted.vcf > {i[:-4]}.mapped-variants.txt")
###print(f"{bedtools} map -c 2,10 -o collapse,collapse -delim \"SPLIT\" -a {window_bed} -b {vcf1[:-4]}.sorted.vcf > {vcf1[:-4]}.mapped-variants.txt")

#4 summarise numbers of het and homs in each window
# could then attempt to stitch the smaller blocks together to make an ROH?
hets = set(['0/1','1/2',
            '0|1','0|2','0|3','0|4','0|5','0|6','0|7','0|8','0|9','0|10',
            '1|0','1|2','1|3','1|4','1|5','1|6','1|7','1|8','1|9','1|10',
            '2',
            '2|0','2|1','2|2','2|3','2|4','2|5','2|6','2|7','2|8','2|9','2|10',
            '3',
            '3|0','3|1','3|2','3|3','3|4','3|5','3|6','3|7','3|8','3|9','3|10',
            '4|1','4|2','4|3','4|4','4|4','4|5','4|6','4|7','4|8','4|9','4|10',
            '5|0','5|1','5|2','5|3','5|4','5|5','5|6','5|7','5|8','5|9','5|10',
            '6|0','6|1','6|2','6|3','6|4','6|5','6|6','6|7','6|8','6|9','6|10',
            '7|0','7|1','7|2','7|3','7|4','7|5','7|6','7|7','7|8','7|9','7|10',
            '8|0','8|1','8|2','8|3','8|4','8|5','8|6','8|7','8|8','8|9','8|10'])

homs = set(['1/1','1|1'])

results = {}
final_dfs = []
for vcf in vcfs_to_process:
    with open(f'{vcf[:-4]}.mapped-variants.txt', 'r+') as inputfile:
        with open(f'{vcf[:-4]}.results.bed', 'w+') as outputfile:
            #write the bed file header containing colour info
            outputfile.write(f'track name = "{vcf[:-4]}" visibility=2 itemRgb = "On"\n')
            for line in inputfile:
                parts = line.strip().split()
                chunk = parts[3]
                chrom = parts[0]
                start = parts[1]
                stop  = parts[2]
                GT_fields = parts[5].split('SPLIT') # SPLIT each conactenated variant
                POS_fields = parts[4].split('SPLIT') # SPLIT each conactenated variant
                #print('GT_fields:', len(GT_fields))
                #print('POS_fields:',len(POS_fields),POS_fields)
                GTs = [i.split(":")[0] for i in GT_fields] # split the 'GT_fields' data on ':' and select the first one, to get a list of GTs
                #print(chunk, len(GTs),GTs,POS_fields)
                # start a list to store result per window
                status_list = []
                for i in GTs:
                    if i in hets:
                        status = 0
                        status_list.append(status)

                    if i in homs:
                        status = 1
                        status_list.append(status)

                #print(chunk,POS_fields[1],status_list)
                stats_count = len(status_list)
                summary = list(run_length.encode(status_list))
                #print(len(summary),line[0], summary)

                # CALCULATE LENGTH OF RUNS
                # sum the second element of each tuple (which is the number of hets/homs in a row).
                # This should = the len of POS_Fields
                number_hets_and_homs = sum([pair[1] for pair in summary])
                running_total = 0 # keep a tab of the total # of observations in a run
                for i in summary:
                    # each i is a summary of a run of het/hom
                    #print(chunk,i,running_total,'len:',len(POS_fields))
                    # calculate the length of each run. Use the running total to slice the
                    # appropriate value from the POS_fields list
                    # output some useful info for a bed file
                    start = POS_fields[running_total ]
                    end = POS_fields[running_total + i[1] -1 ]
                    running_total += i[1]
                    length_in_bases = int(end) - int(start)
                    bases_per_snv = round(length_in_bases/i[1])
                    # Filter bad regions from bed files - i.e. centromeric and pericentric

                    key = {0: 'het', 1: 'hom'}
                    data = [chrom, start, end, length_in_bases, key[i[0]], i[1], bases_per_snv, chunk]
                    #print(blacklist[chrom][0],data)
                    df = pd.DataFrame([data], columns=['chrom', 'start', 'end', 'run', 'zygosity', 'snvs',
                                                       'bases_per_snv','chunk'])
                    df.insert(0,'sample', vcf[:-4])
                    colours = {1: '255,0,0', 0: '0,200,0'} # define the colours for hom or het using the 0 or 1 notation
                    final_dfs.append(df)
                    if i[1] >= 3:
                        bed_file = f"{chrom}\t{start}\t{end}\t{length_in_bases}|{key[i[0]]}|{i[1]}|{bases_per_snv}|{chunk}"\
                               f"\t0\t+\t{start}\t{end}\t{colours[i[0]]}\n"
                        outputfile.write(bed_file)

    #filter each bed file
    os.system(f"{bedtools} intersect -v -a {vcf[:-4]}.results.bed -b roh-regions-to-filter.bed > {vcf[:-4]}.roh.bed")

DF = pd.concat(final_dfs,axis=0)
DF.to_csv(f"FINAL-lengths.txt",sep="\t",header=True)

###create a final df.bed for bed file filtering out centromeric regions. Add all the info to an info column.
DF['info'] = DF[DF.columns[0:]].apply(
    lambda x: '|'.join(x.dropna().astype(str)),
    axis=1)

header = ['chrom','start','end','info']
DF.to_csv(f"{dest_path}/FINAL-lengths.bed",sep="\t",columns=header,index=False)
os.system(f"sed -i '1d' {dest_path}/FINAL-lengths.bed")
os.system(f"{bedtools} intersect -v -a {dest_path}/FINAL-lengths.bed -b roh-regions-to-filter.bed > {dest_path}/all-samples-roh.bed")
# Prepare a list of samples in the EXACT order kp will plot, in order to label correctly
sample_summary_tuples = list(run_length.encode(DF.iloc[:, 0])) # turn sample col into a list and summarise using counter
sample_order = ([pair[0] for pair in sample_summary_tuples])
with open(f"{dest_path}/sample_order.txt",'w+') as sampledatafile:
    for i in sample_order:
        sampledatafile.write(f"{i}\n")

sample_order = []
##Turn the filtered bed file into a text file ingestion to the R script.
with open(f'{dest_path}/all-samples-roh.bed','r+') as inputfile:
    with open(f'{dest_path}/all-samples-roh.txt', 'w+') as outputfile:
        for line in inputfile:
            line = line.strip().split("\t")
            info = line[3].split('|')
            data = f"{line[0]}\t{line[1]}\t{info[3]}\t{info[0]}\t{info[4]}\t{info[5]}\t{info[6]}\t{info[7]}\t{info[8]}\n"
            outputfile.write(data)

# plot the text file to make a karyogram
os.system(f"Rscript plot-ROH-karyotypes.R {dest_path}/all-samples-roh.txt {dest_path}/sample_order.txt {dest_path}/ROH-karyotype-plot.pdf")
#print("deleting the vcfs...")
for i in vcfs_to_process:
    os.remove(i)
# create a zip archive
analysis_files = os.listdir(f"{dest_path}/")
#with zipfile.ZipFile(f'outfile.zip', 'w') as zipper:
#    for file in analysis_files:
#        zipper.write(f"{dest_path}/{file}",compress_type=zipfile.ZIP_DEFLATED)
shutil.make_archive(dest_path, 'zip', dest_path)

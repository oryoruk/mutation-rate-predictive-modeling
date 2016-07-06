#!usr/bin/python

#example command:
#python binned_bootstrapped_rt_mr_corr.py whole acc_nc_auto rt_koren_lb snp_1kg_eur

#libraries:
import subprocess
import sys
import copy
import csv
import pandas as pd
import numpy as np
import scipy as scipy
from scipy.stats import pearsonr
import math
import matplotlib.pyplot as plt
#from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles
import itertools
from copy import deepcopy
from collections import defaultdict
#from datetime import datetime
from random import randint
import warnings
pd.options.display.mpl_style = 'default'
import pickle



#functions:
def bedtools_operation_cmd( query_regions_path, query_regions_file,\
                            regions_to_apply_path, regions_to_apply_file,\
                            kind_of_operation, files_are_sorted = True):
    cmd = ['bedtools']
    cmd.append(kind_of_operation)
    cmd.append('-a')
    cmd.append(query_regions_path + query_regions_file)
    cmd.append('-b')
    cmd.append(regions_to_apply_path + regions_to_apply_file)
    if files_are_sorted: cmd.append('-sorted')
    return cmd

#to-do: make sure this function works as expected
def total_size_of_regions(query_regions_path, query_regions_file):
    size_cmd = ['awk','{SUM += $3-$2} END {print SUM}']
    size_cmd.append(query_regions_path + query_regions_file)
    #note-to-self: handle exceptions here
    total_size_op = subprocess.Popen(size_cmd, stdout=subprocess.PIPE,\
                                                stderr=subprocess.PIPE)
    total_size, err = total_size_op.communicate()
    return int(total_size)

def no_of_mutations_in_regions( query_regions_path, query_regions_file,\
                                mutations_path, mutations_file):
    find_mut_cmd = bedtools_operation_cmd(  query_regions_path,\
                                            query_regions_file,\
                                            mutations_path, mutations_file,\
                                            'intersect')
    intersecting_mutations = subprocess.Popen(find_mut_cmd,\
                                                stdout=subprocess.PIPE,\
                                                stderr=subprocess.PIPE)
    no_of_intersecting_mutations = subprocess.check_output(('wc', '-l'),\
                                        stdin=intersecting_mutations.stdout)
    return int(no_of_intersecting_mutations)

def bedtools_subtract_and_save_results(query_regions_path, query_regions_file,\
                                        regions_to_subtract_path,\
                                        regions_to_subtract_file,\
                                        resulting_regions_path,\
                                        resulting_regions_file):
    bedtools_subtract_cmd = bedtools_operation_cmd(query_regions_path,\
                                                    query_regions_file,\
                                                    regions_to_subtract_path,\
                                                    regions_to_subtract_file,\
                                                    'subtract', False)
    apply_command_and_save_output(bedtools_subtract_cmd,\
                                    resulting_regions_path,\
                                    resulting_regions_file)

def bedtools_intersect_and_save_results(query_regions_path,\
                                        query_regions_file,\
                                        regions_to_intersect_path,\
                                        regions_to_intersect_file,\
                                        resulting_regions_path,\
                                        resulting_regions_file):
    bedtools_intersect_cmd = bedtools_operation_cmd(query_regions_path,\
                                                    query_regions_file,\
                                                    regions_to_intersect_path,\
                                                    regions_to_intersect_file,\
                                                    'intersect')
    apply_command_and_save_output(bedtools_intersect_cmd,\
                                    resulting_regions_path,\
                                    resulting_regions_file)

def apply_command_and_save_output(cmd, destination_path, destination_file):
    #shell way:
    cmd.append('>')
    cmd.append(destination_path + destination_file)
    subprocess.call(' '.join(cmd),shell=True)
    #python way:
    """
    op = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = op.communicate()
    with open(destination_path + destination_file, 'w') as output_file:
        output_file.write(output)
    """

def regions_file_name(regions_name):
    file_prefix = 'reg_'
    file_suffix = '.bed'
    return file_prefix + regions_name + file_suffix

#note: mutations and polymorphisms are point sites in the genome
def sites_file_name(sites_name):
    file_suffix = '.bed'
    return sites_name + file_suffix

#global variables:
no_of_rt_states = 4
list_of_rt_states = ['s'+str(i+1) for i in range(no_of_rt_states)]
no_of_resampling = 10**4
#no_of_resampling = 10**3
#resampling_ratio = 0.9
output_dir = './../../output/initial_analysis/binned_bootstrapped_rt_mr_corr/'
input_dir = './../../data/formatted_regions_sites/'


initial_regions = sys.argv[1]
filters = sys.argv[2]
rt_regions = sys.argv[3]
mutation_dataset = sys.argv[4]

"""
query_regions = 'encode'
rt_regions = 'uva_rt'
mutation_dataset = 'snp_pg_jw'
"""
analysis_folder_suffix =  ''
try:
    analysis_folder_suffix = sys.argv[5]
except:
    pass

analysis_folder = '__'.join([initial_regions, filters, \
                            rt_regions, mutation_dataset])
if analysis_folder_suffix:  analysis_folder += '__' + analysis_folder_suffix
analysis_dir = output_dir + analysis_folder + '/'

print 'creating analysis folder: ' + analysis_folder
subprocess.call('mkdir ' + analysis_dir, shell=True)

query_regions = initial_regions + '__' + filters
query_regions_file = regions_file_name(query_regions)
mutation_dataset_file = sites_file_name(mutation_dataset)

total_no_of_muts = total_size_of_regions(input_dir, \
                                                mutation_dataset_file)
total_sizes_of_all_rt_binned_regs = []
no_of_muts_in_all_rt_binned_regs = []

#binning the regions into n subsets for each rt state
#saving the resulting subsets in respective files
#measuring the total size of each subset
#counting the corresponding mutations in each subset
#note-to-self: might want to keep these regions ready for the script
for rt_state in list_of_rt_states:


    print 'intersecting query regions with ' + rt_state
    rt_state_regs_file = regions_file_name(rt_state + rt_regions.strip('rt') )
    rt_binned_query_regions_file = regions_file_name(query_regions + '__'\
                                        +rt_state+ rt_regions.strip('rt'))
    #using bedtools to find the intersections:
    bedtools_intersect_and_save_results(input_dir, query_regions_file, \
                                        input_dir, rt_state_regs_file, \
                                        analysis_dir, \
                                        rt_binned_query_regions_file)
    #example file name: 'reg_whole__acc_nc_auto__koren_s1.bed'
    print 'saving rt binned query regions for ' + rt_state
    #to-do: make sure this function works as expected
    total_size_of_rt_binned_regs = total_size_of_regions(analysis_dir, \
                                                         rt_binned_query_regions_file)
    no_of_muts_in_rt_binned_regs = no_of_mutations_in_regions(analysis_dir, \
                                                              rt_binned_query_regions_file, \
                                                              input_dir, \
                                                              mutation_dataset_file)
    total_sizes_of_all_rt_binned_regs.append(total_size_of_rt_binned_regs)
    no_of_muts_in_all_rt_binned_regs.append(no_of_muts_in_rt_binned_regs)

print total_no_of_muts
print total_sizes_of_all_rt_binned_regs
print no_of_muts_in_all_rt_binned_regs

#here an if total_no_of_muts > preset_amount can be used to 
#make 




#resampling the mutation dataset to identify the robustness of the analysis
#OPTION 1:
print 'bootstrapping'
#generating samples by bootstrapping:
samples =[]
for temp1 in range(no_of_resampling):
    sample = [0,0,0,0]
    for temp2 in range(total_no_of_muts):
        pick = randint(1,total_no_of_muts)
        for i in range(no_of_rt_states):
            pick -= no_of_muts_in_all_rt_binned_regs[i]
            if pick <= 0:
                sample[i] += 1
                break
    samples.append(sample)
samples =  np.array(samples)
samples =  np.transpose(samples)


with open(analysis_dir + analysis_folder + '__samples.csv', 'w') as samples_file:
    writer = csv.writer(samples_file,delimiter='\t')
    for i in range(no_of_rt_states):
        writer.writerow(samples[i])
#End of OPTION 1:

for i in range(no_of_rt_states):
    win_site_counts_filename = analysis_dir + 'counts_' + str(total_sizes_of_all_rt_binned_regs[i]) + '_' + list_of_rt_states[i] + '.pickle'
    fileObject = open(win_site_counts_filename, 'wb')
    pickle.dump(samples[i], fileObject)
    fileObject.close()
print 'counts for windows saved for rt state ' + rt_state

"""

#OPTION 2:
print 'bootstrapping'
#reading in generated bootstrapping samples from previous runs:
samples = np.zeros([no_of_rt_states, no_of_resampling], dtype=int)
with open(analysis_path + analysis_folder + '__samples.csv','r') as csvfile:
    reader  = csv.reader(csvfile, delimiter = '\t')
    for i,row in enumerate(reader):
        samples[i] = np.array([int(x) for x in row])
#End of OPTION 2:
"""


densities = np.empty([4,no_of_resampling])
for i in range(no_of_rt_states):
    densities[i] = (samples[i] / float(total_sizes_of_all_rt_binned_regs[i] ) )
list_of_densities = [val for sublist in list(densities) for val in sublist]
rt_states = []
for i in range(no_of_rt_states):
    rt_states +=no_of_resampling * [i+1]
norm_densities = densities / np.mean(densities[0])
list_of_norm_densities = [val for sublist in list(norm_densities) for val in sublist]



#reporting statistics on the bootstrapped results
print 'reporting results'
print 'saving 1st plot'
fig = plt.figure(figsize=(20,20))
#plt.boxplot(list(densities), labels=list_of_rt_states, showmeans=True)
ax = fig.add_subplot(111)

notes = ""

slope, intercept, r_value, p_value, slope_std_error = scipy.stats.linregress(rt_states,list_of_densities)

annotation = "Pearson R: "+ str(np.corrcoef(rt_states,list_of_densities)[0,1]) + "\n" \
            "Slope: "+ str(slope) + "\n" \
            "P-value: "+ str(p_value)
notes += annotation + '\n'
#plt.annotate(annotation,xy = (0,0), color = 'b',alpha=0.5,size= 20)

fitted_snp = intercept+ slope*np.array(rt_states)
plt.plot(rt_states, fitted_snp, ':',color='#5378C1', linewidth=10)


## add patch_artist=True option to ax.boxplot() 
## to get fill color
bp = ax.boxplot(list(densities), labels=[i.upper() for i in list_of_rt_states], showmeans=True, patch_artist=True)


## change outline color, fill color and linewidth of the boxes
for box in bp['boxes']:
    # change outline color
    box.set( color='#5C85D6', linewidth=3)
    # change fill color
    box.set( facecolor = '#D6E0F5' )

## change color and linewidth of the whiskers
for whisker in bp['whiskers']:
    whisker.set(color='#7570b3', linewidth=3)

## change color and linewidth of the caps
for cap in bp['caps']:
    cap.set(color='#7570b3', linewidth=3)

## change color and linewidth of the medians
for median in bp['medians']:
    median.set(color='#00B85C', linewidth=3)

## change the style of fliers and their fill
for flier in bp['fliers']:
    flier.set(marker='o', color='#800000', alpha=0.3)

for flier in bp['means']:
    flier.set(marker='o', color='#e7298a', alpha=1)

plot_title = analysis_folder.replace('__',' / ').replace('_',' ')
if plot_title.count('/') == 4 : plot_title = plot_title[:-6] + '(' + plot_title[-4:] + ')'

plt.title(plot_title, fontsize = 40)
plt.ylabel("Mutation Density",fontsize=40) #fontsizep
plt.xlabel("Replication Timing Segments",fontsize=40)
plt.yticks(fontsize = 30)
plt.xticks(fontsize = 30)
#plt.ylim((0,np.amax(densities)))


#supressing mac osx user warning
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plt.tight_layout()

plt.savefig(analysis_dir + analysis_folder)
#plt.show()

for i in range(no_of_rt_states):
    notes += list_of_rt_states[i].upper() + " SNP Density Mean:" +str(np.mean(densities[i])) + '\n'

for sx, sy in itertools.permutations(range(4),2):
    notes +=  "Wilcoxon Rank-Sum Test P-value (S"+str(sx+1)+" & S" + str(sy+1)+ "): " +str(scipy.stats.ranksums(densities[sx],densities[sy] )[1]) + '\n'
    notes +=  "Wilcoxon Signed-Rank Test P-value (S"+str(sx+1)+" & S" + str(sy+1)+ "): " +str(scipy.stats.wilcoxon(densities[sx],densities[sy] )[1]) + '\n'

with open(analysis_dir + analysis_folder + '__notes.txt', 'w') as notes_file:
    notes_file.write(notes)

print 'saving 2nd plot'
fig = plt.figure(figsize=(20,20))
#plt.boxplot(list(densities), labels=list_of_rt_states, showmeans=True)
ax = fig.add_subplot(111)


slope, intercept, r_value, p_value, slope_std_error = scipy.stats.linregress(rt_states,list_of_densities)

fitted_snp = intercept+ slope*np.array(rt_states)
plt.plot(rt_states, fitted_snp, ':',color='#5378C1', linewidth=10)


## add patch_artist=True option to ax.boxplot() 
## to get fill color
bp = ax.boxplot(list(densities), labels=[i.upper() for i in list_of_rt_states], showmeans=True, patch_artist=True)


## change outline color, fill color and linewidth of the boxes
for box in bp['boxes']:
    # change outline color
    box.set( color='#5C85D6', linewidth=3)
    # change fill color
    box.set( facecolor = '#D6E0F5' )

## change color and linewidth of the whiskers
for whisker in bp['whiskers']:
    whisker.set(color='#7570b3', linewidth=3)

## change color and linewidth of the caps
for cap in bp['caps']:
    cap.set(color='#7570b3', linewidth=3)

## change color and linewidth of the medians
for median in bp['medians']:
    median.set(color='#00B85C', linewidth=3)

## change the style of fliers and their fill
for flier in bp['fliers']:
    flier.set(marker='o', color='#800000', alpha=0.3)

for flier in bp['means']:
    flier.set(marker='o', color='#e7298a', alpha=1)

plot_title = analysis_folder.replace('__',' / ').replace('_',' ')
if plot_title.count('/') == 4 : plot_title = plot_title[:-6] + '(' + plot_title[-4:] + ')'

plt.title(plot_title, fontsize = 40)
plt.ylabel("Mutation Density",fontsize=40) #fontsizep
plt.xlabel("Replication Timing Segments",fontsize=40)
plt.yticks(fontsize = 30)
plt.xticks(fontsize = 30)
plt.ylim(0)


#supressing mac osx user warning
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plt.tight_layout()

plt.savefig(analysis_dir + analysis_folder + '__framed.png')
#plt.show()


print 'saving 3rd plot'
fig = plt.figure(figsize=(20,20))
#plt.boxplot(list(densities), labels=list_of_rt_states, showmeans=True)
ax = fig.add_subplot(111)

notes = ""

slope, intercept, r_value, p_value, slope_std_error = scipy.stats.linregress(rt_states,list_of_norm_densities)

annotation = "Pearson R: "+ str(np.corrcoef(rt_states,list_of_norm_densities)[0,1]) + "\n" \
            "Slope: "+ str(slope) + "\n" \
            "P-value: "+ str(p_value)
notes += annotation + '\n'
#plt.annotate(annotation,xy = (0,0), color = 'b',alpha=0.5,size= 20)

fitted_snp = intercept+ slope*np.array(rt_states)
plt.plot(rt_states, fitted_snp, ls=':',color='#5378C1', linewidth=10)


## add patch_artist=True option to ax.boxplot() 
## to get fill color
bp = ax.boxplot(list(norm_densities), labels=[i.upper() for i in list_of_rt_states], showmeans=True, patch_artist=True)


## change outline color, fill color and linewidth of the boxes
for box in bp['boxes']:
    # change outline color
    box.set( color='#5C85D6', linewidth=3)
    # change fill color
    box.set( facecolor = '#D6E0F5' )

## change color and linewidth of the whiskers
for whisker in bp['whiskers']:
    whisker.set(color='#7570b3', linewidth=3)

## change color and linewidth of the caps
for cap in bp['caps']:
    cap.set(color='#7570b3', linewidth=3)

## change color and linewidth of the medians
for median in bp['medians']:
    median.set(color='#00B85C', linewidth=3)

## change the style of fliers and their fill
for flier in bp['fliers']:
    flier.set(marker='o', color='#800000', alpha=0.3)

for flier in bp['means']:
    flier.set(marker='o', color='#e7298a', alpha=1)

plot_title = analysis_folder.replace('__',' / ').replace('_',' ')
if plot_title.count('/') == 4 : plot_title = plot_title[:-6] + '(' + plot_title[-4:] + ')'

plt.title(plot_title, fontsize = 40)
plt.ylabel("Mutation Density",fontsize=40) #fontsizep
plt.xlabel("Replication Timing Segments",fontsize=40)
plt.yticks(fontsize = 30)
plt.xticks(fontsize = 30)
#plt.ylim((0,np.amax(norm_densities)))

#supressing mac osx user warning
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plt.tight_layout()


plt.savefig(analysis_dir + analysis_folder + '__normalized.png')


for i in range(no_of_rt_states):
    notes += list_of_rt_states[i].upper() + " SNP Density Mean:" +str(np.mean(norm_densities[i])) + '\n'

for sx, sy in itertools.permutations(range(4),2):
    notes +=  "Wilcoxon Rank-Sum Test P-value (S"+str(sx+1)+" & S" + str(sy+1)+ "): " +str(scipy.stats.ranksums(norm_densities[sx],norm_densities[sy] )[1]) + '\n'
    notes +=  "Wilcoxon Signed-Rank Test P-value (S"+str(sx+1)+" & S" + str(sy+1)+ "): " +str(scipy.stats.wilcoxon(norm_densities[sx],norm_densities[sy] )[1]) + '\n'

with open(analysis_dir + analysis_folder + '__normalized_notes.txt', 'w') as notes_file:
    notes_file.write(notes)

print 'saving 4th plot'
fig = plt.figure(figsize=(20,20))
#plt.boxplot(list(densities), labels=list_of_rt_states, showmeans=True)
ax = fig.add_subplot(111)


slope, intercept, r_value, p_value, slope_std_error = scipy.stats.linregress(rt_states,list_of_norm_densities)


fitted_snp = intercept+ slope*np.array(rt_states)
plt.plot(rt_states, fitted_snp, ls=':',color='#5378C1', linewidth=10)


## add patch_artist=True option to ax.boxplot() 
## to get fill color
bp = ax.boxplot(list(norm_densities), labels=[i.upper() for i in list_of_rt_states], showmeans=True, patch_artist=True)


## change outline color, fill color and linewidth of the boxes
for box in bp['boxes']:
    # change outline color
    box.set( color='#5C85D6', linewidth=3)
    # change fill color
    box.set( facecolor = '#D6E0F5' )

## change color and linewidth of the whiskers
for whisker in bp['whiskers']:
    whisker.set(color='#7570b3', linewidth=3)

## change color and linewidth of the caps
for cap in bp['caps']:
    cap.set(color='#7570b3', linewidth=3)

## change color and linewidth of the medians
for median in bp['medians']:
    median.set(color='#00B85C', linewidth=3)

## change the style of fliers and their fill
for flier in bp['fliers']:
    flier.set(marker='o', color='#800000', alpha=0.3)

for flier in bp['means']:
    flier.set(marker='o', color='#e7298a', alpha=1)

plot_title = analysis_folder.replace('__',' / ').replace('_',' ')
if plot_title.count('/') == 4 : plot_title = plot_title[:-6] + '(' + plot_title[-4:] + ')'

plt.title(plot_title, fontsize = 40)
plt.ylabel("Mutation Density",fontsize=40) #fontsizep
plt.xlabel("Replication Timing Segments",fontsize=40)
plt.yticks(fontsize = 30)
plt.xticks(fontsize = 30)
plt.ylim(0)

#supressing mac osx user warning
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plt.tight_layout()


plt.savefig(analysis_dir + analysis_folder + '__normalized_framed.png')


print 'saving 5th plot'
fig = plt.figure(figsize=(20,20))
#plt.boxplot(list(densities), labels=list_of_rt_states, showmeans=True)
ax = fig.add_subplot(111)



slope, intercept, r_value, p_value, slope_std_error = scipy.stats.linregress(rt_states,list_of_norm_densities)


#plt.annotate(annotation,xy = (0,0), color = 'b',alpha=0.5,size= 20)

fitted_snp = intercept+ slope*np.array(rt_states)
plt.plot(rt_states, fitted_snp, ls='-',color='#5378C1', linewidth=100)


## add patch_artist=True option to ax.boxplot() 
## to get fill color
bp = ax.boxplot(list(norm_densities), labels=[i.upper() for i in list_of_rt_states], showmeans=True, patch_artist=True)


## change outline color, fill color and linewidth of the boxes
for box in bp['boxes']:
    # change outline color
    box.set( color='#5C85D6', linewidth=3)
    # change fill color
    box.set( facecolor = '#D6E0F5' )

## change color and linewidth of the whiskers
for whisker in bp['whiskers']:
    whisker.set(color='#7570b3', linewidth=3)

## change color and linewidth of the caps
for cap in bp['caps']:
    cap.set(color='#7570b3', linewidth=3)

## change color and linewidth of the medians
for median in bp['medians']:
    median.set(color='#00B85C', linewidth=3)

## change the style of fliers and their fill
for flier in bp['fliers']:
    flier.set(marker='o', color='#800000', alpha=0.3)

for flier in bp['means']:
    flier.set(marker='o', color='#e7298a', alpha=1)

plot_title = analysis_folder.replace('__',' / ').replace('_',' ')
if plot_title.count('/') == 4 : plot_title = plot_title[:-6] + '(' + plot_title[-4:] + ')'

plt.title(plot_title, fontsize = 40)
plt.ylabel("Mutation Density",fontsize=40) #fontsizep
plt.xlabel("Replication Timing Segments",fontsize=40)
plt.yticks(fontsize = 30)
plt.xticks(fontsize = 30)
plt.ylim((0.75,1.65))

                                                   

#supressing mac osx user warning
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plt.tight_layout()

plt.savefig(analysis_dir + analysis_folder + '__normalized_fixed_small.png')


print 'saving 6th plot'
fig = plt.figure(figsize=(20,20))
#plt.boxplot(list(densities), labels=list_of_rt_states, showmeans=True)
ax = fig.add_subplot(111)



slope, intercept, r_value, p_value, slope_std_error = scipy.stats.linregress(rt_states,list_of_norm_densities)


#plt.annotate(annotation,xy = (0,0), color = 'b',alpha=0.5,size= 20)

fitted_snp = intercept+ slope*np.array(rt_states)
plt.plot(rt_states, fitted_snp, ls=':',color='#5378C1', linewidth=10)


## add patch_artist=True option to ax.boxplot() 
## to get fill color
bp = ax.boxplot(list(norm_densities), labels=[i.upper() for i in list_of_rt_states], showmeans=True, patch_artist=True)


## change outline color, fill color and linewidth of the boxes
for box in bp['boxes']:
    # change outline color
    box.set( color='#5C85D6', linewidth=3)
    # change fill color
    box.set( facecolor = '#D6E0F5' )

## change color and linewidth of the whiskers
for whisker in bp['whiskers']:
    whisker.set(color='#7570b3', linewidth=3)

## change color and linewidth of the caps
for cap in bp['caps']:
    cap.set(color='#7570b3', linewidth=3)

## change color and linewidth of the medians
for median in bp['medians']:
    median.set(color='#00B85C', linewidth=3)

## change the style of fliers and their fill
for flier in bp['fliers']:
    flier.set(marker='o', color='#800000', alpha=0.3)

for flier in bp['means']:
    flier.set(marker='o', color='#e7298a', alpha=1)

plot_title = analysis_folder.replace('__',' / ').replace('_',' ')
if plot_title.count('/') == 4 : plot_title = plot_title[:-6] + '(' + plot_title[-4:] + ')'

plt.title(plot_title, fontsize = 40)
plt.ylabel("Mutation Density",fontsize=40) #fontsizep
plt.xlabel("Replication Timing Segments",fontsize=40)
plt.yticks(fontsize = 30)
plt.xticks(fontsize = 30)
plt.ylim((0.75,1.65))

                                                   

#supressing mac osx user warning
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plt.tight_layout()

plt.savefig(analysis_dir + analysis_folder + '__normalized_fixed.png')







"""
#initial plotting functions
########################################

fig = plt.figure(figsize=(20,20))
#plt.boxplot(list(densities), labels=list_of_rt_states, showmeans=True)
ax = fig.add_subplot(111)

notes = ""

slope, intercept, r_value, p_value, slope_std_error = scipy.stats.linregress(rt_states,list_of_densities)

annotation = "Pearson R: "+ str(np.corrcoef(rt_states,list_of_densities)[0,1]) + "\n" \
            "Slope: "+ str(slope) + "\n" \
            "P-value: "+ str(p_value)
notes += annotation + '\n'
#plt.annotate(annotation,xy = (0,0), color = 'b',alpha=0.5,size= 20)

fitted_snp = intercept+ slope*np.array(rt_states)
plt.plot(rt_states, fitted_snp, 'k:',color='#99B2E6', linewidth=10)

## add patch_artist=True option to ax.boxplot() 
## to get fill color
bp = ax.boxplot(list(densities), labels=[i.upper() for i in list_of_rt_states], showmeans=True, patch_artist=True)


## change outline color, fill color and linewidth of the boxes
for box in bp['boxes']:
    # change outline color
    box.set( color='#5C85D6', linewidth=3)
    # change fill color
    box.set( facecolor = '#D6E0F5' )

## change color and linewidth of the whiskers
for whisker in bp['whiskers']:
    whisker.set(color='#7570b3', linewidth=3)

## change color and linewidth of the caps
for cap in bp['caps']:
    cap.set(color='#7570b3', linewidth=3)

## change color and linewidth of the medians
for median in bp['medians']:
    median.set(color='#00B85C', linewidth=3)

## change the style of fliers and their fill
for flier in bp['fliers']:
    flier.set(marker='o', color='#800000', alpha=0.3)

for flier in bp['means']:
    flier.set(marker='o', color='#e7298a', alpha=1)

plot_title = analysis_folder.replace('__',' / ').replace('_',' ')
if plot_title.count('/') == 4 : plot_title = plot_title[:-6] + '(' + plot_title[-4:] + ')'

plt.title(plot_title, fontsize = 40)
plt.ylabel("Mutation Density",fontsize=40) #fontsize
plt.xlabel("Replication Timing Segments",fontsize=40)
plt.yticks(fontsize = 30)
plt.xticks(fontsize = 30)

#supressing mac osx user warning
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plt.tight_layout()

plt.savefig(analysis_path + analysis_folder)
#plt.show()


#saving a 2nd smaller one
fig = plt.figure(figsize=(20,20))
#plt.boxplot(list(densities), labels=list_of_rt_states, showmeans=True)
ax = fig.add_subplot(111)

slope, intercept, r_value, p_value, slope_std_error = scipy.stats.linregress(rt_states,list_of_densities)

#plt.annotate(annotation,xy = (0,0), color = 'b',alpha=0.5,size= 20)

fitted_snp = intercept+ slope*np.array(rt_states)
plt.plot(rt_states, fitted_snp, 'k-', linewidth=40)

## add patch_artist=True option to ax.boxplot() 
## to get fill color
bp = ax.boxplot(list(densities), labels=[i.upper() for i in list_of_rt_states], showmeans=True, patch_artist=True)


## change outline color, fill color and linewidth of the boxes
for box in bp['boxes']:
    # change outline color
    box.set( color='#5C85D6', linewidth=3)
    # change fill color
    box.set( facecolor = '#D6E0F5' )

## change color and linewidth of the whiskers
for whisker in bp['whiskers']:
    whisker.set(color='#7570b3', linewidth=3)

## change color and linewidth of the caps
for cap in bp['caps']:
    cap.set(color='#7570b3', linewidth=3)

## change color and linewidth of the medians
for median in bp['medians']:
    median.set(color='#00B85C', linewidth=3)

## change the style of fliers and their fill
for flier in bp['fliers']:
    flier.set(marker='o', color='#800000', alpha=0.3)

for flier in bp['means']:
    flier.set(marker='o', color='#e7298a', alpha=1)

plot_title = analysis_folder.replace('__',' / ').replace('_',' ')
if plot_title.count('/') == 4 : plot_title = plot_title[:-6] + '(' + plot_title[-4:] + ')'

plt.title(plot_title, fontsize = 40)
plt.ylabel("Mutation Density",fontsize=40) #fontsize
plt.xlabel("Replication Timing Segments",fontsize=40)
plt.yticks(fontsize = 30)
plt.xticks(fontsize = 30)

#supressing mac osx user warning
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plt.tight_layout()

plt.savefig(analysis_path + analysis_folder + '__small.png')


for i in range(no_of_rt_states):
    notes += list_of_rt_states[i].upper() + " SNP Density Mean:" +str(np.mean(densities[i])) + '\n'

for sx, sy in itertools.permutations(range(4),2):
    notes +=  "Wilcoxon Rank-Sum Test P-value (S"+str(sx+1)+" & S" + str(sy+1)+ "): " +str(scipy.stats.ranksums(densities[sx],densities[sy] )[1]) + '\n'
    notes +=  "Wilcoxon Signed-Rank Test P-value (S"+str(sx+1)+" & S" + str(sy+1)+ "): " +str(scipy.stats.wilcoxon(densities[sx],densities[sy] )[1]) + '\n'

with open(analysis_path + analysis_folder + '__notes.txt','w') as notes_file:
    notes_file.write(notes)

with open(analysis_path + analysis_folder + '__samples.csv','w') as samples_file:
    writer = csv.writer(samples_file,delimiter='\t')
    for i in range(no_of_rt_states):
        writer.writerow(samples[i])



########################################
"""


"""
matching_rt_dataset = {'whole': 'koren', 'encode': 'uva'}

analysis = 'test'

if analysis == 'snp':
    list_of_initial_regions = ['whole', 'encode']
    list_of_regions_to_subtract = ['ae', 'rep', 'cpg', '2tss', 'iss', 'cns']
    list_of_regions_to_intersect = ['acc', 'nc', 'auto']
    matching_mutation_datasets = {'whole': ['snp_1kg_eur','snp_1kg_afr','snp_1kg_asn'], 'encode': ['snp_pg_jw'] }

elif analysis == 'dm':
    list_of_initial_regions = ['whole']
    list_of_regions_to_subtract = []
    list_of_regions_to_intersect = ['acc', 'nc', 'auto']
    matching_mutation_datasets = {'whole': ['dm_decode','dm_gonl'] }

elif analysis == 'sm':
    list_of_initial_regions = ['whole']
    list_of_regions_to_subtract = []
    list_of_regions_to_intersect = ['acc', 'nc', 'auto']
    matching_mutation_datasets = {'whole': ['sm_koren'] } #note-to-self: later on add sm_tcga

elif analysis == 'divergence':
    list_of_initial_regions = ['whole']
    list_of_regions_to_subtract = []
    list_of_regions_to_intersect = ['acc', 'nc', 'auto']#note-to-self: add divergence filters
    matching_mutation_datasets = {'whole': ['snp_1kg_eur','snp_1kg_afr','snp_1kg_asn'], 'encode': ['snp_pg_jw'] }

elif analysis == 'test':
    list_of_initial_regions = ['whole', 'encode']
    #list_of_regions_to_subtract = ['ae', 'rep', 'cpg', '2tss', 'iss', 'cns']
    #list_of_regions_to_subtract = ['cpg']
    list_of_regions_to_subtract = ['']
    list_of_regions_to_intersect = ['acc', 'nc', 'auto']
    matching_mutation_datasets = {'whole': ['snp_1kg_eur'], 'encode': ['snp_pg_jw'] }


for initial_regions in list_of_initial_regions:
    print 'creating sub-analysis folder for ' + initial_regions + ' and ' + matching_rt_dataset[initial_regions]
    subprocess.call('mkdir ' + analysis+'/'+initial_regions+'_'+matching_rt_dataset[initial_regions], shell=True)
    current_regions = initial_regions
    print 'starting analysis of ' + current_regions
    subprocess.call('mkdir ' + analysis + '/' + initial_regions+ '_' + matching_rt_dataset[initial_regions]+'/'+current_regions,shell=True)
    #subprocess.call('cp ' file_prefix + initial_regions + file_suffix +' /'+analysis+'/'+file_prefix+initial_regions+file_suffix, shell=True)        
    list_of_rt_state_regions = [matching_rt_dataset[initial_regions] + '_' + sx for sx in rt_states]
    print 'binning ' + current_regions + 'into rt state bins using ' + matching_rt_dataset[initial_regions]
    for rt_state_regions in list_of_rt_state_regions:
	    print 'binning ' + current_regions + ' into ' + rt_state_regions 
        #bin
          
        for mutation_dataset in matching_mutation_datasets[initial_regions]:
        #count mutations in the bin
        #report
            print 'counting snps in the bin'
            pass
    for regions_to_subtract in list_of_regions_to_subtract:
        it_is_first_subtraction = list_of_regions_to_subtract.index(regions_to_subtract) == 0
        if it_is_first_subtraction: subtraction_prefix = "_f"
        else: subtraction_prefix = ''
        resulting_regions = current_regions + subtraction_prefix + str(list_of_regions_to_subtract.index(regions_to_subtract) + 1)
        target_folder =  analysis + '/' + initial_regions+ '_' + matching_rt_dataset[initial_regions]+'/'+resulting_regions + '/'
        subprocess.call('mkdir ' + target_folder ,shell=True)
        bedtools_subtract(current_regions, regions_to_substract, target_folder, resulting_regions)
        current_regions = resulting_regions
        for rt_state_regions in list_of_rt_state_regions:
            #binning into rt dataset
            resulting_regions = current_regions + '_' + rt_state_regions
            bedtools_intersect(current_regions, rt_state_regions, target_folder, resulting_regions)
            pass
            #count total size of the rt state & record the number
            #count snps for each snp dataset in the rt state & record the numbers
        #report using recorded numbers
    for regions_to_intersect in list_of_regions_to_intersect:
        intersection_prefix = '_'
        current_regions = current_regions + intersection_prefix + regions_to_intersect
        subprocess.call('mkdir ' + analysis + '/' + initial_regions+ '_' + matching_rt_dataset[initial_regions]+'/'+current_regions,shell=True)     	
        for rt_state in rt_states:
            pass

"""

#!usr/bin/python

#example command:
#python binned_rt_sr_corr.py whole acc_nc_auto rt_koren_lb snp_1kg_eur 50000


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
import time
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

#to-do
def filter_initial_regions():
    pass


"""
def no_of_sites_in_region(sites, reg_chrom, reg_chrom_start, reg_chrom_end):
    return len(sites[(sites.chrom == reg_chrom) & (sites.chrom_start>=reg_chrom_start)& (sites.chrom_end<=reg_chrom_end)].index)
"""

#this function is updated to handle regions from multiple chromosomes
#function assumes chromosome number to be of type string
def no_of_sites_in_window(sites, reg_chrom, reg_chrom_start, reg_chrom_end, chrom_order, window_starts_in_prev_chrom = False):
    if window_starts_in_prev_chrom:
        #no of sites in this chromosome
        no_of_sites_in_cur_chrom = len(sites[(sites.chrom == reg_chrom) & (sites.chrom_end<=reg_chrom_end)].index)
        prev_chrom = chrom_order[chrom_order.index(str(reg_chrom))-1]
        #no of sites in next chromosome
        no_of_sites_in_prev_chrom = len(sites[(sites.chrom == prev_chrom) & (sites.chrom_start>=reg_chrom_start)].index)
        return no_of_sites_in_prev_chrom + no_of_sites_in_cur_chrom
    else:
        return len(sites[(sites.chrom == reg_chrom) & (sites.chrom_start>=reg_chrom_start)& (sites.chrom_end<=reg_chrom_end)].index)


#global variables:
#order of chromosomes in bed files:
chrom_order = [str(c) for c in range(1,23)]+ ['X','Y']
chrom_order.sort()


#to-do: allow custom output folder
output_dir = './../../output/initial_analysis/rt_sr_corr/'
input_dir = './../../data/formatted_regions_sites/'





#get user input
#example command:
#python binned_rt_sr_corr.py whole acc_nc_auto koren_rt snp_1kg_eur 50000

initial_regions = sys.argv[1]
filters = sys.argv[2]
query_regions = initial_regions + '__' + filters
query_regions_file = regions_file_name(query_regions)
rt_regions = sys.argv[3]
site_dataset = sys.argv[4]
site_dataset_file = sites_file_name(site_dataset)
win_size = int(sys.argv[5])





analysis_folder_suffix = ''
try:
    analysis_folder_suffix = sys.argv[6]
except:
    pass
analysis_folder = '__'.join([initial_regions, filters, \
                             rt_regions, site_dataset])


"""
#e.g.

initial_regions =  'whole'
filters = 'acc_nc_auto'
query_regions = 'whole__acc_nc_auto'
rt_regions = 'koren_rt'
site_dataset = 'snp_1kg_eur'
analysis_folder_suffix = ''
analysis_folder = '__'.join([initial_regions, filters, \
                             rt_regions, site_dataset])

query_regions = initial_regions + '__' + filters
query_regions_file = regions_file_name(query_regions)
site_dataset_file = sites_file_name(site_dataset)
win_size = 50000



#e.g.:
#analysis_folder: 'whole__acc_nc_auto__koren_rt__snp_1kg_eur'



"""

#create analysis folder
if analysis_folder_suffix:  analysis_folder += '__' + analysis_folder_suffix
analysis_dir = output_dir + analysis_folder + '/'
print 'creating analysis folder: ' + analysis_folder
subprocess.call('mkdir ' + analysis_dir, shell=True)


total_no_of_sites = total_size_of_regions(input_dir, \
                                          site_dataset_file)
print 'total number of sites: ' + str(total_no_of_sites)

""""
I want to create the following example files:

reg_whole__acc_nc_auto__koren_rt.bed

snp_1kg_eur__whole__acc_nc_auto__koren_rt.bed
"""

print 'intersecting query regions with ' + rt_regions
rt_regs_file = regions_file_name(rt_regions)
rt_sites_file = sites_file_name(rt_regions)
rt_intersected_query_regions_file = regions_file_name(query_regions + '__'\
                                    + rt_regions)
#using bedtools to find the intersections:
bedtools_intersect_and_save_results(input_dir, query_regions_file, \
                                    input_dir, rt_regs_file, \
                                    analysis_dir, \
                                    rt_intersected_query_regions_file)
#example file name: 'reg_whole__acc_nc_auto__koren_s1.bed'
print 'saving rt intersected query regions for ' + rt_regions
#now intersect these rt binned query regions with sites
sites_in_rt_intersected_query_regions_file = sites_file_name(site_dataset + '__' + query_regions + '__' \
                                                 + rt_regions)
bedtools_intersect_and_save_results(analysis_dir, rt_intersected_query_regions_file, \
                                    input_dir, site_dataset_file, \
                                    analysis_dir, \
                                    sites_in_rt_intersected_query_regions_file)
print 'saving sites that fall into rt intersected query regions for ' + rt_regions

rt_meas_in_rt_intersected_query_regions_file = sites_file_name(rt_regions + '__' + query_regions + '__' \
                                                             + rt_regions)
bedtools_intersect_and_save_results(input_dir, rt_sites_file, \
                                    analysis_dir, rt_intersected_query_regions_file, \
                                    analysis_dir, \
                                    rt_meas_in_rt_intersected_query_regions_file)


print 'saving rt measurements that fall into rt intersected query regions for ' + rt_regions

print total_size_of_regions(analysis_dir, rt_intersected_query_regions_file)
print no_of_mutations_in_regions(analysis_dir, rt_intersected_query_regions_file, \
                                                          input_dir, \
                                                          site_dataset_file)



"""
# forcing the dtype of chrom to be of object
rt_query_sites = pd.read_csv(analysis_dir + sites_in_rt_binned_query_regions_file, delimiter='\t', header=None,
                             names=['chrom', 'chrom_start', 'chrom_end'], dtype={'chrom': object})
rt_query_regs = pd.read_csv(analysis_dir + rt_binned_query_regions_file, delimiter='\t', header=None,
                            names=['chrom', 'chrom_start', 'chrom_end'], dtype={'chrom': object})

cur_win_chrom = '1'
cur_win_len, cur_win_start, cur_win_end = [0] * 3
win_site_counts = []
window_starts_in_prev_chrom = False

for i, reg in rt_query_regs.iterrows():

    cur_reg_len = reg.chrom_end - reg.chrom_start
    if reg.chrom != cur_win_chrom:
        cur_win_chrom = reg.chrom
        window_starts_in_prev_chrom = True

    # we are extending the current window
    # there are 3 scenarios
    # most likely scenario: if the current region length is not exceeding our window size:
    if cur_reg_len <= (win_size - cur_win_len):
        # add the current region to the current window
        cur_win_len += cur_reg_len
        # update the end position of current window
        cur_win_end = reg.chrom_end
    # else if the current region is exceeding our window size:
    else:
        # add the current region to the current window
        # update the end position of current window
        cur_win_end = reg.chrom_start + (win_size - cur_win_len)
        # count sites in the current window
        # print cur_win_chrom, cur_win_start, cur_win_end, window_starts_in_prev_chrom
        win_site_count = no_of_sites_in_window(rt_query_sites, cur_win_chrom, cur_win_start, cur_win_end,
                                               chrom_order, window_starts_in_prev_chrom)
        # if the window started in previous chromosome update it for next window:
        if window_starts_in_prev_chrom:
            # print win_site_count
            window_starts_in_prev_chrom = False

        # append number of sites in this window
        win_site_counts.append(win_site_count)

        # if (len(win_site_counts) %100 )== 0:
        #   print len(win_site_counts)



        # print(win_site_count)

        # to-do: current region might be even larger than the win_size, tackle that scenario

        # update cur_win_len for next iteration:
        cur_win_len = cur_reg_len - (win_size - cur_win_len)
        # update cur_win_start and cur_win_end for next iteration:
        cur_win_start = cur_win_end
        cur_win_end = reg.chrom_end

# pickle win_site_counts + rt_state + win_size
win_site_counts_filename = analysis_dir + 'win_site_counts_' + str(win_size) + '_' + rt_state + '.pickle'
fileObject = open(win_site_counts_filename, 'wb')
pickle.dump(win_site_counts, fileObject)
fileObject.close()

"""
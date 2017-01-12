import sys
sys.path.append("/project/voight_subrate/avarun/Research/mutation_rate")
from modules import *
from function_wrapper import *
import pandas as pd
import itertools
import pickle

def get_seq_context_variant(chr_name,position,before,after):
    #Chr_name: chromosome name
    #position: position of variant
    #before after: nucleotides before or after the variant
    files = file_handles("/project/voight_subrate/avarun/Research/fasta_file_hg19")
    chr_name = str(chr_name)
    position = int(position)
    before = int(before)
    after = int(after)
    str_around = files.fetch_str_var(chr_name,position,before,after) #string around the variant including it
    #print "sequence context "+str(before)+" basepairs before and "+str(after)+" after the variant including the middle position"+''.join(str_around)
    #return ''.join(str_around)
    return ''.join(str_around)

def get_seq_context_interval(chr_name,start_p,end_p,before_after):
    #Chr_name: chromosome name
    #start_p: start position of interval
    #end_p: end postion of interval
    #before_after: nucleotides before or after the start/end (needed for padding)
    files = file_handles("/project/voight_subrate/avarun/Research/fasta_file_hg19")
    chr_name = str(chr_name)
    start_p = int(start_p)
    end_p = int(end_p)
    before_after = int(before_after)
    acceptable_chr_name = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]
    if chr_name in acceptable_chr_name:
        seq = ''.join(files.get_string(chr_name, start_p, end_p, before_after, before_after))
        return seq
    else:
        print "Chromosome entered is not valid"

def reverse_complement(seq):
    complement_code = dict( zip( "ATCGNatcgn" , "TAGCNtagcn" ) )
    #return "".join( complement_code[nuc] for nuc in reversed(seq) )
    return "".join( complement_code[nuc] for nuc in seq[::-1] )

def check_sub(snps,chrom,chrom_start,reverse_strand):
    snp_at_site = snps[(snps.chrom==chrom)&(snps.chrom_start==chrom_start)]
    if snp_at_site.empty:
        return (False, '')
    #if there is a snp:
    else:
        #return sub=True and the alt nucleotide
        #take the reverse complement if working on the reverse strand
        if reverse_strand:
            return (True, reverse_complement(snp_at_site.iloc[0].alt) )
        else:
            return (True, snp_at_site.iloc[0].alt )

def list_without_item(full_list,item):
    return [x for x in full_list if x != item]

def only_main_4_bases(seq):
    #this function returns True for empty strings
    return all(nuc in list('GATC') for nuc in seq)

regions_file_name = \
snp_file_name =
output_file =
total_no_of_jobs = int(sys.argv[1])
job_index = int(sys.argv[2])


context_size = 7
before = (context_size-1)/2
after = (context_size-1)/2

#cpg_file = './../../data/formatted_regions_sites/reg_whole__acc_nc_auto_cpg.bed'
regions_file = './../../data/formatted_regions_sites/reg_whole__nc2_cpg.bed'
snp_file = './../../data/mutation_datasets/snp/all_var_EUR_chr_loc_sorted'
output_file = './../../output/initial_analysis/cpg_ct/' + str(job_index)
#output_file = './../../output/initial_analysis/cpg_ct_rt/' + save_file_name


cpg_columns = ['chrom','chrom_start','chrom_end']

#initialization of dictionary to hold info on all contexts:
sub_rates_dct = {}
#here are the substitutions and corresponding fold I am recording:
subs   = ['CT','CA','CG','AT','AC','AG']
#folds = ['GA','GT','GC','TA','TG','TC']
nucs = list('GATC')
for sub in subs:
    for flanks in itertools.product(nucs,repeat = 6):
        ref = ''.join(flanks[:3]) + sub[0] + ''.join(flanks[3:])
        alt = ''.join(flanks[:3]) + sub[1] + ''.join(flanks[3:])
        sub_rates_dct[ref+alt] = [0,0]


bed_columns = ['chrom','chrom_start','chrom_end']
regions = pd.read_csv(regions_file, delimiter = '\t', header = None, names = bed_columns, dtype ={'chrom':object})

#subset of cpg sites that this job needs to consider:
cpg_sites = cpg_sites[job_index*10**4:(job_index+1)*10**4]

snps = pd.read_csv(snp_file, delimiter = '\t', header = None, names = ['chrom','chrom_start','ref','alt','alt2','id','freq'], dtype ={'chrom':object})
#snps.sort_values(['chrom','chrom_start'],ascending=[1,1], inplace=True)
#snps.reset_index(drop=True,inplace=True)
#subset of snps that this job needs to consider:
slice_start_chrom = cpg_sites.head(1).chrom.iloc[0]
slice_end_chrom = cpg_sites.tail(1).chrom.iloc[0]
slice_start_chrom_start = cpg_sites.head(1).chrom_start.iloc[0]
slice_end_chrom_end = cpg_sites.tail(1).chrom_end.iloc[0]
slice_start_index = snps[(snps.chrom == slice_start_chrom)&(snps.chrom_start >= slice_start_chrom_start)].head(1).index[0]
slice_end_index = snps[(snps.chrom == slice_end_chrom)&(snps.chrom_start < slice_end_chrom_end)].tail(1).index[0] + 1
snps = snps[slice_start_index : slice_end_index ]


exception_counter = 0

#for each region in regions from bed file:
for i, region in regions.iterrows():
    if i == 4:
        break
    #for each site in the region:
    #increment the context count
    #also increment respective substitution count if there is any
    for site in range(region.chrom_start, region.chrom_end):
        #7mer context of the site
        #set reverse strand flag to False, assume not reverse strand
        reverse_strand = False
        site_context = get_seq_context_variant(region.chrom, site, before, after)
        site_nuc = site_context[3]
        #continue if site context consists of only main 4 bases (G,A,T,C)
        if only_main_4_bases(site_context):
            #recording substitutions of C and A
            #considering reverse strands for G and T
            #so if it is not C or A, consider the reverse complement
            #set reverse strand to True
            #and take reverse complement of both the context and the site
            if site_nuc not in ['C','A']:
                reverse_strand = True
                site_context = reverse_complement(site_context)
                site_nuc = site_context[3]
            #incrementing the context count
            #in a loop because same context appears 3 times as alt sequence
            #using a dictionary bc. it is faster to lookup due to hashing
            #than a pandas dataframe which is not hashed
            for nuc  in list_without_item(nucs, site_nuc):
                sub_rates_dct[site_context+site_context[:3]+nuc+site_context[4:]][0] += 1
            #look up if there is a sub in this site
            sub, alt = check_sub(snps,region.chrom,site, reverse_strand)
            #if there is a substitution at this site, increment corresponding sub count
            if sub:
                sub_rates_dct[site_context+site_context[:3]+alt+site_context[4:]][1] += 1

print exception_counter
#save dictionary:
fileObject = open(output_file, 'wb')
pickle.dump(sub_rates_dct, fileObject)
fileObject.close()
#!usr/bin/python



#example command:
#python rt_sr_corr.py whole acc_nc_auto rt_koren_lb snp_1kg_eur 50000


import pandas as pd
import numpy as np
import pickle


output_dir = './../../output/initial_analysis/binned_rt_sr_corr/'
input_dir =  './../../data/formatted_regions_sites/'


#this variables should be user entered:
win_size = 6*10**6
sites = pd.read_csv(input_dir+'snp_1kg_eur.bed',delimiter = '\t',header=None, names = ['chrom','chrom_start','chrom_end'] , dtype={'chrom':object} )
rt_dataset = pd.read_csv(input_dir+'rt_fsu_hesc_bg01.bed',delimiter = '\t',header=None, names = ['chrom','chrom_start','chrom_end','rt_score'] , dtype={'chrom':object} )
chrom = '1'





def add_site_flanking_rt(avg_loc_rt_profile, site_flanking_rt):
    #add rt scores in the flanks of the site
    #add 0 for nan values, add the value if there is a measurement
    avg_loc_rt_profile.sum_rt_score = avg_loc_rt_profile.sum_rt_score + site_flanking_rt.fillna(0).rt_score
    #update the no_of_mes field:
    #if there is no value add 0, if there is a value add 1:
    avg_loc_rt_profile.no_of_mes = avg_loc_rt_profile.no_of_mes + np.logical_not(np.isnan(site_flanking_rt.rt_score)).astype(int)
    return avg_loc_rt_profile


#function gets rt scores within specified coordinates
#linearly interpolates scores for missing coordinates
#checks for large gaps (e.g. centromere) in rt measurements
def get_interpolated_rt_scores(chrom_rt, chrom_start, chrom_end):
    rt_slice = chrom_rt[(chrom_rt.chrom_start>=chrom_start ) & (chrom_rt.chrom_end<=chrom_end)]
    #this case holds if there are no rt measurements in the area
    if rt_slice.empty:
        i1 = chrom_rt[chrom_rt.chrom_start<chrom_start].tail(1).index.values[0]
        i2 = chrom_rt[chrom_rt.chrom_start>chrom_end].head(1).index.values[0]
    else:
        i1 = rt_slice.head(1).index.values[0]
        i2 = rt_slice.tail(1).index.values[0]
    #print i1, i2
    rt_interpolation_distance_threshold = 10**4
    if i1==0 or (chrom_rt.iloc[i1].chrom_start - chrom_rt.iloc[i1-1].chrom_start) > rt_interpolation_distance_threshold:
        there_is_near_rt_before_slice = False
        start_of_slice = chrom_start
    else :
        there_is_near_rt_before_slice = True
        start_of_slice = chrom_rt.iloc[i1-1].chrom_start
    if i2==chrom_rt.tail(1).index.values[0] or (chrom_rt.iloc[i2+1].chrom_start - chrom_rt.iloc[i2].chrom_start) > rt_interpolation_distance_threshold:
        there_is_near_rt_after_slice = False
        end_of_slice = chrom_end
    else :
        there_is_near_rt_after_slice = True
        end_of_slice = chrom_rt.iloc[i2+1].chrom_start
    #interpolate with start_of_slice & end_of_slice:
    flanks_range = pd.DataFrame(np.arange(start_of_slice,end_of_slice+1),columns=['chrom_start'])
    #then get the slice of measurements that fall around the site
    #but get 1 addititonal measurement on either flank, outside the flank
    #so that linear interpolation can work for positions at the start
    #and end of the region
    rt_slice = chrom_rt.iloc[i1-there_is_near_rt_before_slice: i2+1+there_is_near_rt_after_slice].copy()
    #print start_of_slice, end_of_slice
    #print there_is_near_rt_before_slice, there_is_near_rt_after_slice
    #get the relative coordinates
    #drop other columns
    rt_slice.drop(['chrom','chrom_end'],axis=1,inplace=True)
    rt_slice.reset_index(inplace=True,drop=True)
    #create a second dataframe with all positions we want to interpolate for
    #join these two dataframes effectively creating a dataframe
    #with all the coordinates. coordinates for which there is no
    #rt measurement will have NaN as their values.
    rt_slice = pd.merge(flanks_range, rt_slice,on = ['chrom_start'],how='outer')
    #interpolate missing values
    if not there_is_near_rt_after_slice:
        de_extrapolation_mask = rt_slice.fillna(method='backfill').isnull()
    rt_slice = rt_slice.interpolate()
    if not there_is_near_rt_after_slice:
        rt_slice[de_extrapolation_mask] = np.nan
    if there_is_near_rt_before_slice or there_is_near_rt_after_slice:
        #crop:
        #get the slice that we are really interested cropping start and end
        rt_slice = rt_slice[(rt_slice.chrom_start>=chrom_start)&(rt_slice.chrom_start<=chrom_end)]
        rt_slice.reset_index(inplace=True,drop=True)
    return rt_slice

then =  datetime.datetime.now()

avg_loc_rt_profile = pd.DataFrame(np.transpose(np.vstack([np.arange(-win_size/2,win_size/2+1) ,[0]*(win_size+1),[0.0]*(win_size+1)])),columns=['rel_coord','no_of_mes','sum_rt_score'])
#dataframe to hold the rt scores around a site
#initialized to all NaN s
site_flanking_rt = pd.DataFrame(np.transpose(np.vstack([np.arange(-win_size/2,win_size/2+1) ,[np.nan]*(win_size+1)])),columns=['rel_coord','rt_score'])
#chromosomal coordinate of the site
#initialized so that the window is just outside the chromosome
#and all set to NaN
chrom_start = - win_size/2 - 1
#for i, site in sites.iloc[0:1000].iterrows():
for i, site in sites[sites.chrom==chrom].iterrows():
    #how far is this next site: offset
    offset = site.chrom_start - chrom_start
    #update current chrom_start with sites coordinate
    chrom_start = site.chrom_start
    #shift site_flanking_rt by offset:
    site_flanking_rt.rel_coord = site_flanking_rt.rel_coord - offset
    #size of the overlap
    overlap = win_size - offset + 1
    #add new rt_scores to site_flanking_rt (how many?: offset many)
    #print chrom_start-win_size/2+overlap,chrom_start+win_size/2
    new_rt_scores = get_interpolated_rt_scores(chrom_rt, chrom_start-win_size/2+overlap,chrom_start+win_size/2)
    #calculate relative coordinates
    new_rt_scores.chrom_start = new_rt_scores.chrom_start - chrom_start
    #update name of the first column from chrom_start to rel_coord
    new_rt_scores.columns = ['rel_coord','rt_score']
    site_flanking_rt = site_flanking_rt.append(new_rt_scores, ignore_index=True)
    #crop unneeded bits outside the window
    site_flanking_rt = site_flanking_rt[(site_flanking_rt.rel_coord>= -win_size/2)&(site_flanking_rt.rel_coord<=win_size/2)]
    site_flanking_rt.reset_index(inplace=True,drop=True)
    #add this site's rt scores to the running sum dataframe
    add_site_flanking_rt(avg_loc_rt_profile, site_flanking_rt)
    print i

now  =  datetime.datetime.now()
print round((now-then).total_seconds(),2), "seconds passed"


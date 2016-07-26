#This file will find the divergence for a species for each of RT regions. Will save the results as each interval and  "Number of unchanged entries" "Number of divergent entries"

import os
import pdb
import sys
sys.path.append("/project/voight_subrate/avarun/Research/mutation_rate")
from function_wrapper import *
from modules import *

if __name__ == "__main__":
    species = sys.argv[1]
    rt_region = sys.argv[2] #Should be s1 or s2 or s3/s4
    chrom = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]
    #Find which regions to even consider in the first place
    handle = file("raw_data/reg_koren_"+rt_region+".bed_sorted_noncoding_merged_"+species,"r")
    regions_okay = dict()
    content = handle.readlines()
    for entry in content:
        entry = entry.rstrip('\n').split('\t')
	if entry[0] in chrom:
	    if not entry[0] in regions_okay:
		regions_okay[entry[0]] = {}
	    if not (int(entry[4]),int(entry[5])) in regions_okay[entry[0]]:
		regions_okay[entry[0]][ (int(entry[4]),int(entry[5])) ] = [] 
	    regions_okay[entry[0]][(int(entry[4]),int(entry[5]))].append( (int(entry[1]),int(entry[2])) )
    handle_write = file(rt_region+"_"+species+"_divergence_results","w")
    for chromosome in chrom:
    #Read the file and the counts of kmer
      handle = file("../../hg19_"+species+"_alignment/chr"+chromosome+".hg19."+species+".net.axt","r")
      content = handle.readlines()
      index = 0
      while(index<len(content)):
        if not content[index][0] == "#":
            #Current line is position
            vals = content[index].split(' ')
	    start = int(vals[2])
	    end = int(vals[3])
	    #Check if it even belongs to okay regions
	    if (start,end) in regions_okay[chromosome]:
		for (pos_start,pos_end) in regions_okay[chromosome][ (start,end) ] : #For each replicating region inside the aligned region
	            unchanged = 0.0
        	    changed = 0.0
		    counter = 0
        	    human_seq = content[index+1].rstrip('\n')
	            other_seq = content[index+2].rstrip('\n')
		    for position in range(start, end+1): #Each position in the aligned region. Check if it falls in replicating region
			if position >= pos_start and position <= pos_end: #it lies in the replicating region
	                    human_nuc = human_seq[counter].upper()
        	            other_nuc = other_seq[counter].upper()
                            if not "-" in human_nuc and not 'N' in human_nuc:
                                if human_nuc == other_nuc:
                                    unchanged += 1
                                else:
	                            changed += 1
			counter += 1
		    handle_write.write(chromosome+'\t'+str(pos_start)+'\t'+str(pos_end)+'\t'+str(unchanged)+'\t'+str(changed)+'\t'+str(changed/(unchanged+changed))+'\n')
            index += 4
        else:
            index += 1

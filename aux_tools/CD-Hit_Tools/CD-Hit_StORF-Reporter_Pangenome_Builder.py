import collections
import copy
import math
from collections import Counter
import sys
import matplotlib.pyplot as plt
import pandas as pd
from  matplotlib.ticker import PercentFormatter


#clusters_In = open('./E-coli/Escherichia_coli_PEP_UR-Con-StORFs_NO_STAR_CD_80_NEW.clstr','r') # Clusters for single Genera
#clusters_In = open('./E-coli/Escherichia_coli_PEP_UR_StORFs_NO_STAR_CD_80_NEW.clstr','r') # Clusters for single Genera
PEP_In = open('/home/nick/Documents/Single_Genome/E-coli/Filtered/Ecoli_Filtered_PEP_CD_c90_s60.clstr','r')
StORF_In = open('/home/nick/Documents/Single_Genome/E-coli/Filtered/Ecoli_PEP_Clustered_CD_c90_s60_With_Unclustered_UR_StORFs_CD_c90_s60.clstr','r')



clusters = collections.OrderedDict()
pangenome_clusters_PEP = collections.OrderedDict()
pangenome_clusters_PEP_SEQS = collections.OrderedDict()

count = 0
first = True
genome_dict = collections.defaultdict(int)
reps = collections.OrderedDict()
county = 0
clusters_With_Con_StORFs = []
singleton_cluster = []
## Load in all data for easier reuse later
for line in PEP_In:
    if line.startswith('>'):
        if first == False:
            cluster_size = len(clusters[cluster_id])
            reps.update({rep:[cluster_size,len(pangenome_clusters_PEP[cluster_id])]})
            if len(clusters[cluster_id]) == 1 and not singleton_cluster: # Stop at clusters smaller than 10
                singleton_cluster.append(cluster_id)
            #     pangenome_clusters_PEP.popitem()
            #     pangenome_clusters_PEP_SEQS.popitem()

            #     reps.popitem()
            #     break
        Ensem_genomes, Con_genomes = [], []
        cluster_id = line.strip('>')
        cluster_id = cluster_id.strip('\n')
        cluster_id = cluster_id.split(' ')[1]
        clusters.update({cluster_id: []})
        pangenome_clusters_PEP.update({cluster_id:[]})
        pangenome_clusters_PEP_SEQS.update({cluster_id:[]})

        first = False
    else:
        clustered = line.split('\t')[1]
        clustered = clustered.split('>')[1]
        clustered = clustered.split('...')[0]
        genome = clustered.split('|')[0]
        genome_dict[genome] +=1
        if '*' in line:
            rep = clustered
            reps.update({rep:[0,0]})
        if first == False:
            clusters[cluster_id].append(clustered)
            clustered_genome = clustered.split('|')[0]

            if clustered_genome not in pangenome_clusters_PEP[cluster_id]:
                pangenome_clusters_PEP[cluster_id].append(clustered_genome)
            pangenome_clusters_PEP_SEQS[cluster_id].append(clustered)

######################################

#penguins = sns.load_dataset("penguins")




######################################
Combined_pangenome_clusters_PEP = collections.OrderedDict()
Combined_pangenome_clusters_PEP_SEQS = collections.OrderedDict()

Combined_pangenome_clusters_StORF = collections.OrderedDict()
Combined_pangenome_clusters_StORF_SEQS = collections.OrderedDict()

Combined_pangenome_clusters_PEP_StORF_Clustered = collections.OrderedDict()

not_StORF_Only_Cluster_IDs = []

already_seen_PEP = []

Combined_clusters = collections.OrderedDict()
Combined_reps = collections.OrderedDict()
first = True
## We load in the combined PEP and StORF-Reporter data separately
for line in StORF_In:
    if line.startswith('>'):
        if first == False:
            cluster_size = len(Combined_clusters[cluster_id])
            Combined_reps.update({rep: cluster_size})
            for pep in Combined_pangenome_clusters_PEP_SEQS[cluster_id]:
                if pep != []:
                    if pep in already_seen_PEP:# and Combined_pangenome_clusters_PEP_SEQS[cluster_id] != []:
                        #print("Already Counted")
                        continue
                    else:
                        already_seen_PEP.append(pep)
          #  if len(Combined_pangenome_clusters_PEP[cluster_id]) > 1:
          #      print("Here")
            if len(Combined_pangenome_clusters_StORF_SEQS[cluster_id]) > 0 and len(Combined_pangenome_clusters_PEP_SEQS[cluster_id]) > 0:
                if len(Combined_pangenome_clusters_PEP_SEQS[cluster_id]) > 1: # If we have clustered >1 PEP family, we need to record 1 as key and all others are val
                    all_but_first = Combined_pangenome_clusters_PEP_SEQS[cluster_id][1:]
                    storfs_clustered = Combined_pangenome_clusters_StORF_SEQS[cluster_id]
                    VALUE = all_but_first+storfs_clustered
                else:
                    VALUE = Combined_pangenome_clusters_StORF_SEQS[cluster_id]
                KEY = Combined_pangenome_clusters_PEP_SEQS[cluster_id][0]
                Combined_pangenome_clusters_PEP_StORF_Clustered.update({KEY:VALUE})
            # if len(Combined_clusters[cluster_id]) > 1 : # Stop at end of file
            #     Combined_pangenome_clusters_PEP.popitem()
            #     Combined_pangenome_clusters_PEP_SEQS.popitem()
            #     Combined_pangenome_clusters_StORF_SEQS.popitem()
            #     Combined_pangenome_clusters_Con_StORF.popitem()
            #     Combined_reps.popitem()
            #     break
        cluster_id = line.strip('>')
        cluster_id = cluster_id.strip('\n')
        cluster_id = cluster_id.split(' ')[1]
        Combined_clusters.update({cluster_id: []})
        Combined_pangenome_clusters_PEP.update({cluster_id:[]})
        Combined_pangenome_clusters_PEP_SEQS.update({cluster_id:[]})
        Combined_pangenome_clusters_StORF_SEQS.update({cluster_id: []})
        Combined_pangenome_clusters_StORF.update({cluster_id: []})
        first = False
    else:
        clustered = line.split('\t')[1]
        clustered = clustered.split('>')[1]
        clustered = clustered.split('...')[0]
        genome = clustered.split('|')[0]
        genome_dict[genome] +=1
        if '*' in line:
            rep = clustered
            Combined_reps.update({rep:0})
        if first == False:
            Combined_clusters[cluster_id].append(clustered)
            clustered_genome = clustered.split('|')[0]
            if "StORF_Type=StORF" in line:
                if cluster_id not in clusters_With_Con_StORFs: # For counting?
                    clusters_With_Con_StORFs.append(cluster_id)
                if clustered_genome not in Combined_pangenome_clusters_StORF[cluster_id]:
                    Combined_pangenome_clusters_StORF[cluster_id].append(clustered_genome)
                Combined_pangenome_clusters_StORF_SEQS[cluster_id].append(clustered)
            else:
                if cluster_id not in not_StORF_Only_Cluster_IDs:
                    not_StORF_Only_Cluster_IDs.append(cluster_id)# Tell us which StORF-Reporter clustered are unmatched to a PEP
                if clustered_genome not in Combined_pangenome_clusters_PEP[cluster_id]:
                    Combined_pangenome_clusters_PEP[cluster_id].append(clustered_genome)
                Combined_pangenome_clusters_PEP_SEQS[cluster_id].append(clustered)

####################


ALL_PEPS = 0


StORF_Add_To_Singletons = []
Com_PEPs = 0
Com_PEPs_List = []
total_combined_pep_clusters = []
only_PEP = []
num_clustered_PEP = collections.defaultdict(list)
recorded_PEP = []
pangenome_clusters_Type = copy.deepcopy(pangenome_clusters_PEP)
list_of_reps = list(reps.keys())
pep_max_comb = [0,0]
for cluster, pep_genomes in pangenome_clusters_PEP.items():
    rep = list_of_reps[int(cluster)] # get the rep of the current pep cluster
    #if rep not in recorded_PEP:
    recorded_PEP.append(rep)
    Com_PEP_Genomes = 0
    StORFs = 0
    seen_StORFs = []
    Added_StORF_Genomes = 0
    try: # get the cluster from the storf clusters which contains this rep
        clustered_combined =  Combined_pangenome_clusters_PEP_StORF_Clustered[rep] # Not true clusters - I put a PEP as key myself
        seen_clust_Genomes = []
        num_clustered_PEP[cluster].append(rep+'_'+str(len(pep_genomes)))

        for clust in clustered_combined:
            if 'StORF_Type=StORF' not in clust: # Not good enough at the moment
                ### Need to get the number of pep genomes for each pep clustered into this
                Com_PEPs +=1 #
                recorded_PEP.append(clust)
                if (rep+'_'+str(len(pep_genomes))) not in Com_PEPs_List:
                    Com_PEPs_List.append(rep+'_'+str(len(pep_genomes)))
                clust_Genome = clust.split('|')[0]
                if clust_Genome not in seen_clust_Genomes:
                    seen_clust_Genomes.append(clust_Genome)
                    if clust_Genome not in pep_genomes:
                        Com_PEP_Genomes +=1
                num_clustered_PEP[cluster].append(clust+'_'+str(reps[clust][1]))
            elif 'StORF_Type=StORF' in clust:
                StORFs +=1
                clust_Genome = clust.split('|')[0]
                if clust_Genome not in seen_StORFs:
                    seen_StORFs.append(clust_Genome)
                if clust_Genome not in seen_clust_Genomes:
                    seen_clust_Genomes.append(clust_Genome)
                    if clust_Genome not in pep_genomes:
                        Added_StORF_Genomes +=1
            else:
                print("WHAT")

        if len(pep_genomes) == 1 and Added_StORF_Genomes >= 10:
            StORF_Add_To_Singletons.append([cluster,Added_StORF_Genomes])

        size_of_pep_clusters = []
        peps =  num_clustered_PEP[cluster]
        if len(peps) > 1:
            total_combined_pep_clusters.append(peps)
        for pep in peps:
            pep = pep.rsplit('_',1)
            #if pep[0] not in recorded_PEP: # Do not record PEPs are combined and on their own.
            size_of_pep_clusters.append(int(pep[1]))
            ALL_PEPS += int(pep[1])
            recorded_PEP.append(pep[0])
            # else:
            #     print("W")
        pangenome_clusters_Type[cluster] = [len(num_clustered_PEP[cluster]), sum(size_of_pep_clusters), size_of_pep_clusters, Added_StORF_Genomes, StORFs, len(seen_StORFs)]
        if len(peps) > 1:
            print("")
            if Added_StORF_Genomes > 10:
                print("")
        if len(peps) > 2:
            print("")
        if len(peps) > pep_max_comb[0]:
            pep_max_comb = [len(peps),peps]
        num_pep_combined = sum(size_of_pep_clusters)
        if len(peps) >1 and all(i >= 5 for i in size_of_pep_clusters):
            print("")
    except KeyError:
        ###Singleton
        num_pep_genomes = [len(pep_genomes)]
        ALL_PEPS += len(pep_genomes)
        pangenome_clusters_Type[cluster] = [1, len(pep_genomes), num_pep_genomes, Added_StORF_Genomes, StORFs, len(seen_StORFs)]
        only_PEP.append(cluster)

print(pep_max_comb)
# else:
#     print("PEP already recorded in another cluster: " + rep)
#####################
#### write out Cluster IDs for clusters of specific types?? To then be eggnogged
#sys.exit("Killed Ere")
Without_StORF = open('./Ensem_Clusters_Without_StORFs','w')
Without_StORF_min2 = open('./Ensem_Clusters_Without_StORFs_min2','w')
With_StORF_min2 = open('./Ensem_Clusters_With_StORFs_min2','w')
With_Extending_StORF = open('./Ensem_Clusters_With_Extending_StORFs','w')

# number of clustered PEP - Number of PEP genomes - Induvidual sizes of PEP clusters - Added StORF Genomes - StORF sequences - not sure
# [len(num_clustered_PEP[cluster]), sum(size_of_pep_clusters), size_of_pep_clusters, Added_StORF_Genomes, StORFs,len(seen_StORFs)]

for key, value in pangenome_clusters_Type.items():
    if value[4] == 0 and value[1] >=2:
        Without_StORF_min2.write(str(key)+',')
    elif value[4] == 0:
        Without_StORF.write(str(key)+',')
    elif value[3] != 0:
        With_Extending_StORF.write(str(key)+',')
    elif value[4] >=2:
        With_StORF_min2.write(str(key) + ',')
With_StORF_min2.close()
Without_StORF.close()
Without_StORF_min2.close()
With_Extending_StORF.close()



STORF_ONLY = []
ENSEM_ONLY = []
COMBINED_ONLY = []






#######################################

core_99 = math.floor(9.9/10 * len(genome_dict))
core_95 = math.floor(9.5/10 * len(genome_dict))
core_90 = math.floor(9/10 * len(genome_dict))
core_15 = math.floor(1.5/10 * len(genome_dict))

cores = collections.OrderedDict({'pep_core_99':0,'pep_core_95':0,'pep_core_15':0,'extended_99':0,'extended_95':0
         ,'extended_15':0,'comb_extended_99':0,'comb_extended_95':0,'comb_extended_15':0,'storf_core_99':0,'storf_core_95':0,'storf_core_15':0,
                                 'only_storf_core_99':0,'only_storf_core_95':0,'only_storf_core_15':0})

multi_PEP_Combined_By_StORFs = collections.OrderedDict()
multi_PEP_Combined_By_StORFs_1 = collections.OrderedDict()
multi_PEP_Combined_By_StORFs_num_of_PEP_Clusters = 0

StORF_Seqs_Extended = []
StORF_Genomes_Extended = []

record_all_pep_15 = []

core_list = []
soft_core_list = []
accessory_list = []

storf_core_only = []

list_all_pep_numbers = collections.defaultdict(int)


############################ Count PEP separately first to get TRUE Ensembl gene families
def calc_pep_only_core(pep_num):
    list_all_pep_numbers[pep_num] +=1
    if pep_num >= math.floor(core_99):# and StORF_num == 0:
        cores['pep_core_99'] += 1
    elif pep_num >= math.floor(core_95) and pep_num < math.floor(core_99):# and StORF_num == 0:
        cores['pep_core_95'] += 1
    # elif pep_num >= math.floor(core_90) and pep_num < math.floor(core_95):# and StORF_num == 0:
    #     cores['pep_core_90'] += 1
    if pep_num >= math.floor(core_15) and pep_num < math.floor(core_95):# and StORF_num == 0:  # this catch captures some from pep_core_90
        cores['pep_core_15'] += 1
        record_all_pep_15.append(pep_num)
    #####################
def calc_single_pep_extended_StORF_only_core(cluster,pep_num,storf_num): # Count gene families extended with StORFs
    if pep_num < math.floor(core_99) and pep_num != 0 and pep_num+storf_num >= math.floor(core_99):
        cores['extended_99'] +=1
        core_list.append(cluster)
    elif pep_num < math.floor(core_95) and pep_num != 0 and pep_num+storf_num >= math.floor(core_95) and pep_num+storf_num < math.floor(core_99):
        cores['extended_95'] +=1
        soft_core_list.append(cluster)
    # elif pep_num < math.floor(core_90) and pep_num != 0 and pep_num+storf_num >= math.floor(core_90) and pep_num+storf_num < math.floor(core_95):
    #     cores['extended_90'] +=1
    if pep_num < math.floor(core_15) and pep_num != 0 and pep_num+storf_num >= math.floor(core_15) and pep_num+storf_num < math.floor(core_95):
        cores['extended_15'] +=1
        accessory_list.append(cluster)
#####################################
def calc_multi_pep_extended_StORF_only_core(pep_num,storf_num): # Count seperately those gene families extended with StORF-Reporter but combined >1 PEP
    if pep_num < math.floor(core_99) and pep_num != 0 and pep_num+storf_num >= math.floor(core_99):
        cores['comb_extended_99'] +=1
    elif pep_num < math.floor(core_95) and pep_num != 0 and pep_num+storf_num >= math.floor(core_95) and pep_num+storf_num < math.floor(core_99):
        cores['comb_extended_95'] +=1
    # elif pep_num < math.floor(core_90) and pep_num != 0 and pep_num+storf_num >= math.floor(core_90) and pep_num+storf_num < math.floor(core_95):
    #     cores['comb_extended_90'] +=1
    if pep_num < math.floor(core_15) and pep_num != 0 and pep_num+storf_num >= math.floor(core_15) and pep_num+storf_num < math.floor(core_95):
        cores['comb_extended_15'] +=1
######################### StORFs Only >>><<<
def calc_StORF_only_core(storf_num):
    if storf_num >= math.floor(core_99):# and StORF_num == 0:
        cores['storf_core_99'] += 1
    elif storf_num >= math.floor(core_95) and storf_num < math.floor(core_99):# and StORF_num == 0:
        cores['storf_core_95'] += 1
    # elif storf_num >= math.floor(core_90) and storf_num < math.floor(core_95):# and StORF_num == 0:
    #     cores['storf_core_90'] += 1
    if storf_num >= math.floor(core_15) and storf_num < math.floor(core_95):# and StORF_num == 0:  # this catch captures some from pep_core_90
        cores['storf_core_15'] += 1
###########################
def calc_only_StORF_only_core(cluster,storf_num): # only count the true storf onlies
    if storf_num >= math.floor(core_99):# and StORF_num == 0:
        cores['only_storf_core_99'] += 1
        storf_core_only.append(cluster)
    elif storf_num >= math.floor(core_95) and storf_num < math.floor(core_99):# and StORF_num == 0:
        cores['only_storf_core_95'] += 1
    # elif storf_num >= math.floor(core_90) and storf_num < math.floor(core_95):# and StORF_num == 0:
    #     cores['only_storf_core_90'] += 1
    if storf_num >= math.floor(core_15) and storf_num < math.floor(core_95):# and StORF_num == 0:  # this catch captures some from pep_core_90
        cores['only_storf_core_15'] += 1

record_all_pep = []

counter = 0

Number_Of_StORF_Extending_But_Same_Genomes = 0


print("Running")
for cluster, numbers in pangenome_clusters_Type.items(): # put limits here to make sure storf and enembl only are in more than one genome.
    if numbers[3] >= 1:
        StORF_Genomes_Extended.append(numbers[3])
    if numbers[4] >= 1:
        StORF_Seqs_Extended.append(numbers[4])
############################### Calc PEP only
    ######### TO fix - Only loop through the first 1's to get the baseline pep numbs?
    if numbers[0] == 1 and numbers[1] >=2: # If StORFs did not combine PEP reps
        calc_pep_only_core(numbers[1])#,numbers[3])
        counter +=1
    elif numbers[0] >1 and numbers[1] >=2: # IF StORFs combined multiple PEP
        calc_pep_only_core(numbers[2][0])
        counter += 1
        for num in numbers[2]:
            multi_PEP_Combined_By_StORFs_num_of_PEP_Clusters +=1# ,numbers[3])
            multi_PEP_Combined_By_StORFs_1.update({cluster: numbers})
############################# Calc PEP and StORF-Reporter - M
    if numbers[0] == 1 and numbers[3] >= 1: # If StORFs did not combine PEP reps
        calc_single_pep_extended_StORF_only_core(cluster,numbers[1],numbers[3])
    elif numbers[0] >1 and numbers[3] >= 1: # IF unique StORFs combined multiple PEP
        #grouped_pep = sum(numbers[2])
        #for num in numbers[2]:
        calc_multi_pep_extended_StORF_only_core(numbers[1],numbers[3])
        multi_PEP_Combined_By_StORFs.update({cluster:numbers})
    elif numbers[4] >= 1:
        Number_Of_StORF_Extending_But_Same_Genomes +=1
        # for num in numbers[2]:
        #     multi_PEP_Combined_By_StORFs_num_of_PEP_Clusters +=1
###########################

############# Last bit is the storf only stuff - need to calc for all storfs and storfs which had a pep but..
#... dont calc pep towards total score

# c = 0
# for numbers in multi_PEP_Combined_By_StORFs_1.values():
#     for peps in numbers[3]:
#         c +=1
# print(str(c))


############## Typing for the StORF-Reporter-Data
Combined_pangenome_clusters_ONLY_StORF_Type = collections.defaultdict(list)
Combined_pangenome_clusters_StORF_Type = collections.defaultdict(list)
for cluster, genomes in Combined_pangenome_clusters_StORF.items():
    if cluster in not_StORF_Only_Cluster_IDs:
        Combined_pangenome_clusters_StORF_Type[cluster] = [cluster,len(genomes)]
    else:
        Combined_pangenome_clusters_ONLY_StORF_Type[cluster] = [cluster,len(genomes)]


StORF_Only_Out = open("./E-coli_StORF_Only_Clusters_min2",'w')
span = []
#
for cluster, data in Combined_pangenome_clusters_StORF_Type.items():
    #if data[1] >= 2:
    calc_StORF_only_core(data[1])  # ,numbers[3])multi_PEP_Combined_By_StORFs
for cluster, data in Combined_pangenome_clusters_ONLY_StORF_Type.items():
    if data[1] >= 2:
        STORF_ONLY.append(data[1])
        calc_only_StORF_only_core(cluster,data[1])  # ,numbers[3])
        StORF_Only_Out.write(str(cluster) + ',')
        span.append(data[1])


        #if data[1] >= core_99:
        #     os.system("python3 Extract_FASTA_From_Cluster.py -f ./E-coli/Escherichia_coli_PEP_Clustered_With_Unclustered_UR-StORFs.fa_CD_c90_s60.fa  "
        #              "-c ./E-coli/Escherichia_coli_PEP_Clustered_With_Unclustered_UR-StORFs.fa_CD_c90_s60.clstr -id "+ str(data[0]) + " -o ./E-coli/E-coli_Con-StORF_Only_Clusters_To_Be_Swissed.fa")
        #

print("End")


print(cores)

#print(extended)

print(dict(sorted(list_all_pep_numbers.items())))

#[len(num_clustered_PEP[cluster]), sum(size_of_pep_clusters), size_of_pep_clusters, Added_StORF_Genomes, StORFs,len(seen_StORFs)]

pangenome_clusters_Type_list = []
for clust, counts in pangenome_clusters_Type.items():
    if counts[1] >= 169:
        counts[1] = 169
    ENSEM_ONLY.append(counts[1])
    if counts[3] != 0:
        tmp = counts[1]+counts[3]
        COMBINED_ONLY.append(tmp)
    if (counts[1] >=1 and counts[3] >=1):# and counts[0] >1 and counts[1] >1:
        if counts[1] >= 169:
            counts[1] = 169 # Not fake - it controls the combined cluster problem

    clust[1] = clust[0]+clust[1]


    #writerow(COMBINED_ONLY)

Ensem_P = []
StORF_P = []
Combined_P = []

print(len(ENSEM_ONLY))
val = 1
try:
    while True:
        ENSEM_ONLY.remove(val)
except ValueError:
    pass

print(len(ENSEM_ONLY))




for i in range(170):
    if i > 1:
        E = round((ENSEM_ONLY.count(i)/len(ENSEM_ONLY))*100,2)
        Ensem_P.append(E)
        S = round((STORF_ONLY.count(i)/len(STORF_ONLY))*100,2)
        StORF_P.append(S)
        C = round((COMBINED_ONLY.count(i) / len(COMBINED_ONLY)) * 100,2)
        Combined_P.append(C)


Ensem_tmp = [Ensem_P[x:x+10] for x in range(0, len(Ensem_P),10)]
Ensem_Dec = []
for En in Ensem_tmp:
    Ensem_Dec.append(sum(En))


StORF_tmp = [StORF_P[x:x+10] for x in range(0, len(StORF_P),10)]
StORF_Dec = []
for St in StORF_tmp:
    StORF_Dec.append(sum(St))

import csv
import matplotlib.ticker as mtick
print("Here")

with open('./Ensembl_Then_Con-StORF.csv', 'w') as out_file:
    out_file.write("Ensembl:")
    for E in Ensem_P:
        out_file.write(str(E)+',')
    out_file.write("\nCombined:")
    for C in Combined_P:
        out_file.write(str(C)+',')
    out_file.write("\nStORF:")
    for S in StORF_P:
        out_file.write(str(S)+',')



print("DDD")


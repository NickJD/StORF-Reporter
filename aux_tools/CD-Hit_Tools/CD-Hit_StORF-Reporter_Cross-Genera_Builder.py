import collections
import copy
import math
import sys
import numpy as np
from itertools import chain

def get_Genus(clustered):
    clustered_genus = clustered.split('|')[0]
    if '_' in clustered_genus[0]:  # Remove name error
        clustered_genus = clustered_genus.split('_')[1]
    else:
        clustered_genus = clustered_genus.split('_')[0]
    return str(clustered_genus).capitalize()

def get_Species(clustered):
    clustered_species = clustered.split('|')[0]
    if '_' in clustered_species[0]:  # Remove name error
        clustered_species = clustered_species.split('_')[1]
    else:
        clustered_species = clustered_species.split('_')[:2]
    return str('_'.join(clustered_species)).capitalize()


PEP_In = open('/home/nick/Documents/Single_Genome/All_Ensembl_PEP_CD_Clustered_90_60.clstr','r')
StORF_In = open('/home/nick/Documents/Single_Genome/All_Ensem_PEP_CD_Clustered_90_60_Unclustered_UR_StORFs_AA_CD.clstr','r') # Clusters for single Genera

clusters = collections.OrderedDict()


pangenome_clusters_PEP_Genera = collections.OrderedDict()
pangenome_clusters_PEP_Species = collections.OrderedDict()
pangenome_clusters_PEP_Strains = collections.OrderedDict()
pangenome_clusters_PEP_SEQS = collections.OrderedDict()


max_storf_only_genera = 0




count = 0
first = True
genome_dict = collections.defaultdict(int)
reps = collections.OrderedDict()
county = 0
#singleton_cluster = "Null"
clusters_With_Con_StORFs = []
## Load in all data for easier reuse later
for line in PEP_In:
    if line.startswith('>'):
        if first == False:
            Ensem_Con = set(Ensem_genomes).intersection(Con_genomes)
            cluster_size = len(clusters[cluster_id])
            reps.update({rep: [cluster_size,len(pangenome_clusters_PEP_Genera[cluster_id])]}) # Add strains, species here if wanted
            #if len(clusters[cluster_id]) == 1 and "Null" not in singleton_cluster: # Stop at clusters smaller than 10
            #    singleton_cluster = cluster_id
            #if len(clusters[cluster_id]) < 10: # Stop at clusters smaller than 10
                # pangenome_clusters_PEP_Species.popitem()
                # pangenome_clusters_PEP_Genera.popitem()
                # pangenome_clusters_PEP_SEQS.popitem()
                # reps.popitem()
            # if len(clusters[cluster_id]) == 1:
            #    break    # REMEMBER
        Ensem_genomes, Con_genomes = [], []
        cluster_id = line.strip('>')
        cluster_id = cluster_id.strip('\n')
        cluster_id = cluster_id.split(' ')[1]
        clusters.update({cluster_id: []})
        pangenome_clusters_PEP_Genera.update({cluster_id: []})
        pangenome_clusters_PEP_Species.update({cluster_id:[]})
        pangenome_clusters_PEP_Strains.update({cluster_id: []})
        # pangenome_clusters_PEP_SEQS.update({cluster_id:[]})

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
            clustered_genus = get_Genus(clustered)
            clustered_species = get_Species(clustered)
            clustered_strain = clustered.split('|')[0]

            if clustered_genus not in pangenome_clusters_PEP_Genera[cluster_id]:
                pangenome_clusters_PEP_Genera[cluster_id].append(clustered_genus)
            #if clustered_species not in pangenome_clusters_PEP_Species[cluster_id]:
            #    pangenome_clusters_PEP_Species[cluster_id].append(clustered_species)
            if genome not in pangenome_clusters_PEP_Strains[cluster_id]:
                pangenome_clusters_PEP_Strains[cluster_id].append(genome)
            # pangenome_clusters_PEP_SEQS[cluster_id].append(clustered)
print("PEP DONE")
######################################
Combined_pangenome_clusters_PEP_Genera = collections.OrderedDict()
Combined_pangenome_clusters_PEP_Species = collections.OrderedDict()
Combined_pangenome_clusters_PEP_Strains = collections.OrderedDict()
Combined_pangenome_clusters_PEP_SEQS = collections.OrderedDict()

Combined_pangenome_clusters_StORF_Genera = collections.OrderedDict()
Combined_pangenome_clusters_StORF_Species = collections.OrderedDict()
Combined_pangenome_clusters_StORF_Strains = collections.OrderedDict()
Combined_pangenome_clusters_StORF_SEQS = collections.OrderedDict()

Combined_pangenome_clusters_PEP_StORF_Clustered_Genera = collections.OrderedDict()
Combined_pangenome_clusters_PEP_StORF_Clustered = collections.OrderedDict()

not_StORF_Only_Cluster_IDs = []

Combined_clusters = collections.OrderedDict()
Combined_reps = collections.OrderedDict()
first = True
###############
## We load in the combined PEP and StORF_Reporter data separately
for line in StORF_In:
    if line.startswith('>'):
        if first == False:
            cluster_size = len(Combined_clusters[cluster_id])
            Combined_reps.update({rep: cluster_size})
            # if len(Combined_pangenome_clusters_PEP_SEQS[cluster_id]) > 1:
            #     print("Here")
            if len(Combined_pangenome_clusters_StORF_SEQS[cluster_id]) > 0 and len(Combined_pangenome_clusters_PEP_SEQS[cluster_id]) > 0:
                if len(Combined_pangenome_clusters_PEP_SEQS[cluster_id]) > 1: # If we have clustered >1 PEP family, we need to record 1 as key and all others are val
                    all_but_first = Combined_pangenome_clusters_PEP_SEQS[cluster_id][1:]
                    storfs_clustered = Combined_pangenome_clusters_StORF_SEQS[cluster_id]
                    VALUE = all_but_first+storfs_clustered
                else:
                    VALUE = Combined_pangenome_clusters_StORF_SEQS[cluster_id]
                KEY = Combined_pangenome_clusters_PEP_SEQS[cluster_id][0]
                Combined_pangenome_clusters_PEP_StORF_Clustered.update({KEY:VALUE})


            ## Below needs to be rewritten. With >1 genus - be able to record multiple  PEPs for each combined...
            # if len(Combined_pangenome_clusters_StORF_Genera[cluster_id]) > 0 and len(Combined_pangenome_clusters_PEP_Genera[cluster_id]) > 0:
            #     KEY = Combined_pangenome_clusters_PEP_SEQS[cluster_id][0]
            #     VALUE = Combined_pangenome_clusters_StORF_SEQS[cluster_id]
            #     Combined_pangenome_clusters_PEP_StORF_Clustered_Genera.update({KEY:VALUE})
            if len(Combined_clusters[cluster_id]) == 1: # Stop at clusters smaller than 10
                print("First Singleton Cluster is: " +str(cluster_id))
                break
        cluster_id = line.strip('>')
        cluster_id = cluster_id.strip('\n')
        cluster_id = cluster_id.split(' ')[1]
        Combined_clusters.update({cluster_id: []})
        # Combined_pangenome_clusters_PEP_Genera.update({cluster_id:[]})
        # Combined_pangenome_clusters_PEP_Species.update({cluster_id: []})
        Combined_pangenome_clusters_PEP_Strains.update({cluster_id: []})
        Combined_pangenome_clusters_PEP_SEQS.update({cluster_id: []})
        #
        Combined_pangenome_clusters_StORF_Genera.update({cluster_id: []})
        # Combined_pangenome_clusters_StORF_Species.update({cluster_id: []})
        Combined_pangenome_clusters_StORF_Strains.update({cluster_id: []})
        Combined_pangenome_clusters_StORF_SEQS.update({cluster_id: []})
        first = False
    else:
        clustered = line.split('\t')[1]
        clustered = clustered.split('>')[1]
        clustered = clustered.split('...')[0]
        if '*' in line:
            rep = clustered
            Combined_reps.update({rep:0})
        if first == False:
            Combined_clusters[cluster_id].append(clustered)
            clustered_genus = get_Genus(clustered)
            clustered_species = get_Species(clustered)
            clustered_strain = clustered.split('|')[0]
            if '_' in clustered_strain[0]:  # Remove name error
                clustered_strain = clustered_strain.split('_')[1]

            if "StORF_Type" in line:
            #     if cluster_id not in clusters_With_Con_StORFs: # For counting?
            #         clusters_With_Con_StORFs.append(cluster_id)
                if clustered_genus not in Combined_pangenome_clusters_StORF_Genera[cluster_id]:
                    Combined_pangenome_clusters_StORF_Genera[cluster_id].append(clustered_genus)
            #     if clustered_species not in Combined_pangenome_clusters_StORF_Species[cluster_id]:
            #         Combined_pangenome_clusters_StORF_Species[cluster_id].append(clustered_species)
                if clustered_strain not in Combined_pangenome_clusters_StORF_Strains[cluster_id]:
                    Combined_pangenome_clusters_StORF_Strains[cluster_id].append(clustered_strain)
                Combined_pangenome_clusters_StORF_SEQS[cluster_id].append(clustered)
            #
            else:
            #     if clustered_genus not in Combined_pangenome_clusters_PEP_Genera[cluster_id]:
            #         Combined_pangenome_clusters_PEP_Genera[cluster_id].append(clustered_genus)
            #     if clustered_species not in Combined_pangenome_clusters_PEP_Species[cluster_id]:
            #         Combined_pangenome_clusters_PEP_Species[cluster_id].append(clustered_species)
                if clustered_strain not in Combined_pangenome_clusters_PEP_Strains[cluster_id]:
                    Combined_pangenome_clusters_PEP_Strains[cluster_id].append(clustered_strain)
                if cluster_id not in not_StORF_Only_Cluster_IDs:
                    not_StORF_Only_Cluster_IDs.append(cluster_id)  # Tell us which StORF_Reporter clustered are unmatched to a PEP
                Combined_pangenome_clusters_PEP_SEQS[cluster_id].append(clustered)



###HERE for tomorrow - copy the updated work from single to here and repeat for genus,species and strain
list_of_reps = list(reps.keys())
num_clustered_PEP_Genera = collections.defaultdict(list)
recorded_PEP = []
################################# Genera
pangenome_clusters_Type_Genera = copy.deepcopy(pangenome_clusters_PEP_Genera)
pangenome_clusters_Type_Strains = collections.defaultdict(list)

for cluster, pep_genomes in pangenome_clusters_PEP_Genera.items():
    recorded_PEP.append(cluster)
    rep = list_of_reps[int(cluster)]
    Com_PEPs = 0
    Com_PEP_Genomes = 0
    StORFs = 0
    Added_StORF_Genera = 0
    seen_clust_Strains = []

    PEP_Strains = pangenome_clusters_PEP_Strains[cluster]
    for clustered_strain in PEP_Strains:
        if '_' in clustered_strain[0]:  # Remove name error
            clustered_strain = clustered_strain[1:]
        if clustered_strain not in seen_clust_Strains:
            seen_clust_Strains.append(clustered_strain)


    try:
        clustered_combined = Combined_pangenome_clusters_PEP_StORF_Clustered[rep]
        seen_clust_Genera = []
        num_clustered_PEP_Genera[cluster].append(rep + '_' + str(len(pep_genomes)))
        for clust in clustered_combined:
            if 'StORF_Type' not in clust:
                ### Need to get the number of pep genomes for each pep clustered into this
                Com_PEPs += 1
                clustered_genus = get_Genus(clust)
                #clust_Genome = clust.split('|')[0]
                if clustered_genus not in seen_clust_Genera:
                    seen_clust_Genera.append(clustered_genus)
                    if clustered_genus not in pep_genomes:
                        Com_PEP_Genomes += 1
                try:
                    num_clustered_PEP_Genera[cluster].append(clust + '_' + str(reps[clust][1]))
                except TypeError:
                    sys.exit("Broken")

            elif 'StORF_Type' in clust:
                StORFs += 1
                clustered_genus = get_Genus(clust)
                #clust_Genome = clust.split('|')[0]
                if clustered_genus not in seen_clust_Genera:
                    seen_clust_Genera.append(clustered_genus)
                    if clustered_genus not in pep_genomes:
                        Added_StORF_Genera += 1
            else:
                print("WHAT")

        size_of_pep_clusters = []
        peps = num_clustered_PEP_Genera[cluster]
        for pep in peps:
            pep = pep.rsplit('_', 1)
            size_of_pep_clusters.append(int(pep[1]))
        pangenome_clusters_Type_Genera[cluster] = [len(num_clustered_PEP_Genera[cluster]), sum(size_of_pep_clusters),
                                            size_of_pep_clusters, Added_StORF_Genera, StORFs]
        pangenome_clusters_Type_Strains[cluster] = seen_clust_Strains
    except KeyError:
        ###Singleton
        num_pep_genomes = [len(pep_genomes)]
        pangenome_clusters_Type_Genera[cluster] = [1, len(pep_genomes), num_pep_genomes, Added_StORF_Genera, StORFs]
        pangenome_clusters_Type_Strains[cluster] = seen_clust_Strains

print("S")


#######################################

Without_StORF = open('./Ensem_Clusters_Without_StORFs_To_Be_Nogged_min2','w')
With_StORF = open('./Ensem_Clusters_With_StORFs_To_Be_Nogged','w')
#With_Extending_StORF = open('./Ensem_Clusters_With_Extending_Con-StORFs_To_Be_Nogged','w')

for key, value in pangenome_clusters_Type_Genera.items():
    pep_strains = pangenome_clusters_Type_Strains[key]
    if value[4] == 0 and len(pep_strains) >=2:
        Without_StORF.write(str(key)+',')
 #   elif value[3] != 0:
 #       With_Extending_StORF.write(str(key)+',')
 #       With_StORF.write(str(key) + ',')
    elif value[4] >=1:
        With_StORF.write(str(key) + ',')

With_StORF.close()
Without_StORF.close()
#With_Extending_StORF.close()

############## Typing for the StORF_Reporter-Data


multi_PEP_Combined_By_StORFs = collections.OrderedDict()

StORF_Seqs_Extended = []
StORF_Genomes_Extended = []

####################################
#cores = collections.OrderedDict({'pep_genera_single':[],'pep_genera_multi':[],'extended_genera':[],'comb_extended_genera_single':[],'comb_extended_genera_multi':[],'extended_genera_single':[],'extended_genera_multi':0,'storf_genera_single':0,'storf_genera_multi':0,
#                                 'only_storf_genera_single':0,'only_storf_genera_multi':0})

cores = collections.OrderedDict({'pep_genera':[],'extended_genera_single_pep':[],'many_extended_genera_pep':[],'extended_genera':[],'comb_extended_genera':[],'storf_genera':[],'only_storf_genera':[],'only_storf_genera_recording':[]})

extended = collections.OrderedDict()
############################

clsuters_to_be_validated = collections.defaultdict(list)


############################
def calc_pep_only(pep_num):
    cores['pep_genera'].append(pep_num)
    # if pep_num == 1:# and StORF_num == 0:
    #     cores['pep_genera_single'] += 1
    # elif pep_num > 1:# and StORF_num == 0:
    #     cores['pep_genera_multi'] += 1
##########################
def calc_pep_extended_StORF(cluster,pep_num,storf_num):
    if pep_num != 0 and storf_num >= 1:
        cores['extended_genera'].append(pep_num+storf_num)
        clsuters_to_be_validated['extended_genera'].append(cluster)
    if pep_num != 0 and storf_num >= 10:
        cores['many_extended_genera_pep'].append([cluster,pep_num+storf_num])

    if pep_num == 1 and storf_num >= 1:
        cores['extended_genera_single_pep'].append([cluster,pep_num + storf_num])
    #     cores['extended_genera_single'] +=1
    # if pep_num != 0 and storf_num > 1:
    #     cores['extended_genera_multi'] +=1
##########################
def calc_multi_pep_extended_StORF(cluster,number_of_pep_clustered,pep_num,storf_num):
    if pep_num !=0 and storf_num >= 1:
        cores['comb_extended_genera'].append(pep_num+storf_num)
        clsuters_to_be_validated['comb_extended_genera'].append(cluster)

#########################
def calc_StORF_only_when_with_pep(cluster,storf_num):
    cores['storf_genera'].append(storf_num)
    clsuters_to_be_validated['storf_genera'].append(cluster)
    # if storf_num == 1:# and StORF_num == 0:
    #     cores['storf_genera_single'] += 1
    # elif storf_num > 1:# and StORF_num == 0:
    #     cores['storf_genera_multi'] += 1
########################  What is the difference with these?
def calc_only_StORF(cluster,storf_num,max_storf_only_genera): # only count the true storf onlies
    cores['only_storf_genera'].append(storf_num)
    clsuters_to_be_validated['only_storf_genera'].append(cluster)
    if storf_num>=6:
        cores['only_storf_genera_recording'].append([cluster, storf_num])
    if storf_num > max_storf_only_genera:
        max_storf_only_genera = storf_num
    # if storf_num == 1:# and StORF_num == 0:
    #     cores['only_storf_genera_single'] += 1
    # elif storf_num > 1:# and StORF_num == 0:
    #     cores['only_storf_genera_multi'] += 1
    return max_storf_only_genera
#########################



###########################
print("Running")
check_all_calced = 0
for cluster, numbers in pangenome_clusters_Type_Genera.items():
    pep_strains = pangenome_clusters_Type_Strains[cluster]
    if numbers[3] >=1:
        StORF_Genomes_Extended.append(numbers[3])
    if numbers[4] >=1:
        StORF_Seqs_Extended.append(numbers[4])
############################### Calc PEP only
    if numbers[0] == 1 and len(pep_strains) >= 2: # If StORFs did not combine PEP reps
        calc_pep_only(numbers[1])#,numbers[3])
        check_all_calced +=1
    elif numbers[0] >1: # IF StORFs combined multiple PEP
        calc_pep_only(numbers[2][0])
        check_all_calced += 1
        # for num in numbers[2]:
        #     calc_pep_only(num)  # ,numbers[3])

############################# Calc PEP and StORF_Reporter
    if numbers[0] == 1 and numbers[3] >1: # If StORFs did not combine PEP reps
        calc_pep_extended_StORF(cluster,numbers[1],numbers[3])
        extended.update({cluster:numbers})
        check_all_calced += 1
    elif numbers[0] >1 and numbers[3] >1: # IF StORFs combined multiple PEP - Genera added
        #grouped_pep = sum(numbers[2])
        #for num in numbers[2]:
        calc_multi_pep_extended_StORF(cluster,numbers[2],numbers[1],numbers[3]) # same here
        print("combined: " + str(cluster))

        extended.update({cluster: numbers})
        check_all_calced += 1
    elif numbers[0] >1 and numbers[4] >1: # IF StORFs combined multiple PEP
        multi_PEP_Combined_By_StORFs.update({cluster: numbers})


import os
###########################
############################### Calc StORF_Reporter only
Combined_pangenome_clusters_ONLY_StORF_Type = collections.defaultdict(list)
Combined_pangenome_clusters_StORF_Type = collections.defaultdict(list)

biggest_genera = ""
big_genera = 0
biggest_strains = ""
big_strains = 0


#Without_StORF = open('./Ensem_Clusters_Without_Con-StORFs_To_Be_Nogged_min2','w')
#With_StORF = open('./Ensem_Clusters_With_Con-StORFs_To_Be_Nogged','w')
#With_Extending_StORF = open('./Ensem_Clusters_With_Extending_Con-StORFs_To_Be_Nogged','w')
StORF_Only = open("./StORF_Only_Clusters_To_Be_Nogged_min2",'w')

for cluster, genera in Combined_pangenome_clusters_StORF_Genera.items():
    storf_strains = Combined_pangenome_clusters_StORF_Strains[cluster]
    pep_strains = Combined_pangenome_clusters_PEP_Strains[cluster]
    if cluster in not_StORF_Only_Cluster_IDs:
        Combined_pangenome_clusters_StORF_Type[cluster] = [cluster,len(genera)]
        #if len(genera) >= 1:
        calc_StORF_only_when_with_pep(cluster,len(genera))  # ,numbers[3])
    else:
        if len(storf_strains) >= 2:
            StORF_Only.write(str(cluster) + ',')
            Combined_pangenome_clusters_ONLY_StORF_Type[cluster] = [cluster,len(genera)]
            max_storf_only_genera = calc_only_StORF(cluster,len(genera),max_storf_only_genera)
            if len(genera) > big_genera:
                big_genera = len(genera)
                biggest_genera = cluster
            if len(storf_strains) >= big_strains:
                big_strains = len(storf_strains)
                biggest_strains = cluster








print("Biggest: " +biggest_genera)



############################### Calc StORF_Reporter only
# for cluster, data in Combined_pangenome_clusters_StORF_Type.items():
#     if data[1] >=1:
#         calc_StORF_only_when_with_pep(data[1])  # ,numbers[3])
#
#


#################
print(cores)
#print(extended)

from collections import Counter

#print(Counter(cores['pep_genera']))
#print(Counter(cores['extended_genera']))
#print(Counter(cores['comb_extended_genera']))
#print(Counter(cores['storf_genera']))
print(Counter(cores['only_storf_genera']))
print(cores['only_storf_genera_recording'])

print("END")

### Emoty file ready for interesting storfs
# interesting_out = "./StORF_Only_Clusters_To_Be_Swissed.fa"
# with open(interesting_out, 'r+') as f:
#     f.truncate(4)
# for cluster, data in Combined_pangenome_clusters_ONLY_StORF_Type.items():
#     #if number >1:
#     if data[1] >=1:
#         calc_only_StORF(data[1])  # ,numbers[3])
#         # if data[1] >= 2:
#         #     print("Interesting:" + str(cluster))
#         #     os.system(
#         #         "python3 Extract_FASTA_From_Cluster.py -f ./All_Ensem_Filtered_PEP_Clustered_With_Unclustered_UR-StORFS_s.fa_CD_c90_s60.fa  "
#         #         "-c  ./All_Ensem_Filtered_PEP_Clustered_With_Unclustered_UR-StORFS_s.fa_CD_c90_s60.clstr -id " + str(
#         #             data[0]) + " -o "+ interesting_out)
# #
#




#
# ################################## Species
# pangenome_clusters_Type_Species = copy.deepcopy(pangenome_clusters_PEP_Species)
# for cluster, genomes in pangenome_clusters_PEP_Species.items():
#     print(str(len(genomes)) + '\t' + str(len(pangenome_clusters_StORF_Species[cluster])))
#     Con_StORFs = pangenome_clusters_StORF_Species[cluster]
#     unique_con = 0
#     all_con = 0
#     for con in Con_StORFs:
#         all_con +=1
#         if con not in genomes:
#             unique_con +=1
#     pangenome_clusters_Type_Species[cluster] = [len(genomes),all_con,unique_con]
# ################################# Strains
# pangenome_clusters_Type_Strains = copy.deepcopy(pangenome_clusters_PEP_Strains)
# for cluster, genomes in pangenome_clusters_PEP_Strains.items():
#     print(str(len(genomes))+'\t'+str(len(pangenome_clusters_StORF_Strains[cluster])))
#     Con_StORFs = pangenome_clusters_StORF_Strains[cluster]
#     unique_con = 0
#     all_con = 0
#     for con in Con_StORFs:
#         all_con +=1
#         if con not in genomes:
#             unique_con +=1
#     pangenome_clusters_Type_Strains[cluster] = [len(genomes),all_con,unique_con]
# ###################################
# Chris_Out = open('./Chris_Clusters.txt','w')
#
# clusters_For_Chris = collections.OrderedDict()
# clusters_For_Chris_PEP_0 = collections.OrderedDict()
#
# Chris_Out.write("Cluster\tSize\tEnsem_Genera_Num\tCon-StORF_Genera_Num\tCon-StORF_Only_Genera_Num\tEnsem_Species_Num\tCon-StORF_Species_Num\tCon-StORF_Only_Species_Num\tEnsem_Strain_Num\tCon-StORF_Strain_Num\tCon-StORF_Only_Strain_Num\n")
# #This for-loop will go through ALL Clusters allowing for the extraction of ALL different groupings
# for cluster, data in clusters.items():
#     genera_numbers = pangenome_clusters_Type_Genera[cluster]
#     species_numbers = pangenome_clusters_Type_Species[cluster]
#     strain_numbers = pangenome_clusters_Type_Strains[cluster]
#
#     Chris_Out.write(str(cluster)+'\t'+str(len(data))+'\t'+str(genera_numbers[0])+'\t'+str(genera_numbers[1])+'\t'+str(genera_numbers[2])+'\t'+str(species_numbers[0])+'\t'
#                                           +str(species_numbers[1])+'\t'+str(species_numbers[2])+'\t'+str(strain_numbers[0])+'\t'+str(strain_numbers[1])+'\t'+str(species_numbers[2])+'\n')

    # if cluster in clusters_With_Con_StORFs:
    #     print("Current")
    #     size_Of_Cluster = len(clusters[cluster])
    #     ensem_Num = 0
    #     con_StORF_Num = 0
    #     for i in clusters[cluster]:
    #         print(i)
    #         if 'Con-Stop' in i:
    #             con_StORF_Num +=1
    #         else:
    #             ensem_Num +=1
    #     clusters_For_Chris.update({cluster:[size_Of_Cluster,pep_Num,con_StORF_Num,numbers[0],numbers[1]]})
    #     ############# Add - Num of
    #     Chris_Out.write(str(cluster)+'\t'+str(size_Of_Cluster)+'\t'+str(pep_Num)+'\t'+str(ensem_Genera)+          str(con_StORF_Num)+'\t'+str(numbers[0])+'\t'+str(numbers[1])+'\n')
    #     if pep_Num == 0:
    #         clusters_For_Chris_PEP_0.update({cluster:[size_Of_Cluster,pep_Num,con_StORF_Num,numbers[0],numbers[1]]})

print("Da Da!!!!")


#
#
#
#
# ###################################
#
# core_99 = 9.9/10 * len(genome_dict)
# core_95 = 9.5/10 * len(genome_dict)
# core_90 = 9/10 * len(genome_dict)
# core_15 = 1.5/10 * len(genome_dict)
#
# pep_core_99 = 0
# pep_core_95 = 0
# pep_core_90 = 0
# pep_core_15 = 0
#
#
# extended_99 = 0
# extended_95 = 0
# extended_90 = 0
# extended_15 = 0
#             ############### Needs to be redone with new 'numbers'
# for cluster, numbers in pangenome_clusters_Type_Genera.items():
#     if numbers[0] >= math.floor(core_99) and numbers[1] == 0:
#         pep_core_99 +=1
#     elif numbers[0] >= math.floor(core_95) and numbers[0] < math.floor(core_99) and numbers[1] == 0:
#         pep_core_95 +=1
#     elif numbers[0] >= math.floor(core_90) and numbers[0] < math.floor(core_95) and numbers[1] == 0:
#         pep_core_90 +=1
#     if numbers[0] >= math.floor(core_15) and numbers[0] < math.floor(core_95) and numbers[1] == 0: # this catch captures some from pep_core_90
#         pep_core_15 +=1
#     ############ With Con-StORFs
#     if numbers[0] < math.floor(core_99) and numbers[0] != 0 and numbers[0]+numbers[1] >= math.floor(core_99):
#         extended_99 +=1
#     elif numbers[0] < math.floor(core_95) and numbers[0] != 0 and numbers[0]+numbers[1] >= math.floor(core_95) and numbers[0]+numbers[1] < math.floor(core_99):
#         extended_95 +=1
#     elif numbers[0] < math.floor(core_90) and numbers[0] != 0 and numbers[0]+numbers[1] >= math.floor(core_90) and numbers[0]+numbers[1] < math.floor(core_95):
#         extended_90 +=1
#     if numbers[0] < math.floor(core_15) and numbers[0] != 0 and numbers[0]+numbers[1] >= math.floor(core_15) and numbers[0]+numbers[1] < math.floor(core_95):
#         extended_15 +=1
#
# print("Out")
# print(pep_core_99)
# print(pep_core_95)
# print(pep_core_90)
# print(pep_core_15)
#
# print(extended_99)
# print(extended_95)
# print(extended_90)
# print(extended_15)
#

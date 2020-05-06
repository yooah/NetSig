#!/usr/bin/python
# ---------------------------------------------------------------------------

"""
Run target (phenotype) permutation. The phenotype is permuted across samples.

usage: target_permute_NETPHIX.py [-h] [-sol ILP_SOL_FILE]
                                  [-map GENE_MAP_FILE] [-idx INDEX_NAME]
                                  [-sep SEP] [-add_file ADD_FILE]
                                  alt_file target_file correlation k
                                  filter_high filter_low density norm
                                  netfile_name permutations permutefile_name

positional arguments:
  alt_file              alteration data file (alteration matrix, see
                        data/Alterataion_data/ directory for an example)
  target_file           target profile data file (target (phenotype) see
                        data/Target_data/ directory for an example)
  correlation           positive or negative association
  k                     size of modules
  filter_high           upperbound mutation frequency of a gene to be included
  filter_low            lowerbound mutation frequency of a gene to be included
  density               module density
  norm                  target normalization method, z or zlog
  netfile_name          network file
  permutations          num of permutations
  permutefile_name      file name to write the permutation results

optional arguments:
  -h, --help            show this help message and exit
  -sol ILP_SOL_FILE, --ILP_sol_file ILP_SOL_FILE
                        file with optimal ILP solution with target permutation
                        pvalue.
  -map GENE_MAP_FILE, --gene_map_file GENE_MAP_FILE
                        file with ENSG ID to gene name mapping
  -idx INDEX_NAME, --index_name INDEX_NAME
                        index name to use in the target file
  -sep SEP, --sep SEP   separator in file (\t for default)
  -add_file ADD_FILE, --add_file ADD_FILE
                        additional mutation file. take OR with the original
                        mutations

"""


import cplex
from cplex.exceptions import CplexSolverError
import pandas as pd
import time
import sys
import os
import random
import networkx as nx
import argparse
import numpy as np

import netphix_utils as netphix
import permute_utils as perm

timer_start = time.time()

# fixed params
num_itr = 1                             # NEPHLIX finds multiple modules iteratively if num_itr > 1
min_permute = 100                        # minimum # permutations for dynamic test
min_pv = 0.5                            # stop permutations if pv > 0.5 after 100 runs
penalty = "p"                           # use penalty option for mutual exclusivity)

#################
# read parameters


parser = argparse.ArgumentParser()
# required arguments
parser.add_argument("alt_file", help="alteration data file (alteration matrix, see data/Alterataion_data/ directory for an example)")
parser.add_argument("target_file", help="target profile data file (target (phenotype) see data/Target_data/ directory for an example)")
parser.add_argument("correlation", help="positive or negative association")
parser.add_argument("k", type=int, help="size of modules")
parser.add_argument("filter_high", type=float, help="upperbound mutation frequency of a gene to be included")
parser.add_argument("filter_low", type=float, help="lowerbound mutation frequency of a gene to be included")
parser.add_argument("density", type=float, help="module density")
parser.add_argument("norm", help="target normalization method, z or zlog")
parser.add_argument("netfile_name", type=str, help="network file", default= "data/HumanStringNet.txt")
parser.add_argument("permutations", type=int, help="num of permutations")
parser.add_argument("permutefile_name", help="file name to write the permutation results")
parser.add_argument("-sol", "--ILP_sol_file", type=str,
                    help="file with optimal ILP solution with target permutation pvalue.")
parser.add_argument("-map", "--gene_map_file",
                    help="file with ENSG ID to gene name mapping")
parser.add_argument("-idx", "--index_name",
                    help="index name to use in the target file")
parser.add_argument("-sep", "--sep",
                    help="separator in file (\\t for default)")
parser.add_argument("-add_file", "--add_file",
                    help="additional mutation file. take OR with the original mutations")
parser.add_argument("-count_file", "--mut_count_file",
                    help="filename with mutation count info")
parser.add_argument("-mut_th", "--hypermutated_th",type=int,
                    help="mutation count threshold to remove hypermutated samples")
parser.add_argument("-subtype_file", "--subtype_file",
                    help="filename with subtype info")
parser.add_argument("-subtype", "--subtype",
                    help="brca subtype -- 'LumA', 'LumB', 'Basal', or 'Her2'")


args = parser.parse_args()


alteration_file = args.alt_file
target_file = args.target_file
correlation = args.correlation
k = args.k
filter_high = args.filter_high
filter_low = args.filter_low
density = args.density
norm = args.norm
permutations = args.permutations
net_file = args.netfile_name
permutefile_name = args.permutefile_name
solutionfile_name=args.ILP_sol_file
gene_map_file = args.gene_map_file
idx = args.index_name
add_file = args.add_file
subtype=args.subtype
hypermutated_th = args.hypermutated_th

if args.mut_count_file is not None:
    mut_count_file = args.mut_count_file
else:
    mut_count_file = "data/total_mut_count_per_sample.tsv"

if args.subtype_file is not None:
    subtype_file = args.subtype_file
else:
    subtype_file = "data/BRCA-EU.Expression.Subtyping.txt"

if args.sep is not None:
    sep = args.sep
else:
    sep ="\t"

#################
# read files
timer_start = time.clock()
sys.stdout.write("reading target, alteration, net files.. %f\n" %timer_start)
target_df = pd.read_csv(target_file, delimiter=sep, index_col=0)
alt_df = pd.read_csv(alteration_file, delimiter=sep, index_col=0)
all_net = nx.read_edgelist(net_file, data=(('weight',float),))

if add_file is not None:
    brcaness = pd.read_csv(add_file, sep=sep, index_col=0)
    alt_df = alt_df.combine_first(brcaness)
    alt_df.loc[brcaness.index] = np.maximum(alt_df.loc[brcaness.index], brcaness)

if gene_map_file is not None: # map Ensembl ID to gene name
    map_df = pd.read_csv(gene_map_file, sep="\t", index_col=0)
    map_dic = dict(zip(map_df.Ensembl, map_df.Name))
    alt_df = alt_df.rename(map_dic)


#################
# processing data
sys.stdout.write("processing data.. \n" )

print(target_df.shape)

if hypermutated_th is not None:
    print("threshold is "+str(hypermutated_th))
    target_df = netphix.remove_hypermutated(target_df, mut_count_file, hypermutated_th)
if subtype is not None:
    target_df = netphix.filter_subtype(target_df, subtype_file, subtype)

print(target_df.shape)

target_df, alt_df, samples, num_samples  = netphix.preproc_data(target_df, alt_df, filter_high, filter_low)

# exit if there are genes appearing multiple times
# duplicates = alt_df.index[alt_df.index.duplicated()]
# if len(duplicates) > 0:
#     sys.exit("multiple entries for "+",".join(duplicates))

# todo: check. remove duplicate elements
alt_df = alt_df.loc[~alt_df.index.duplicated()]

# normalize
norm_target_df = netphix.norm_target(target_df, norm, correlation)


#################
# prepare to construct ILP model
sys.stdout.write("preparing to create cplex model.. \n" )

mutated_gene_sets, num_alterations =  netphix.proc_alt(alt_df)

edge_lists, mut_lists, num_genes = netphix.proc_net(alt_df.index, all_net, num_alterations)

if idx is None:
    weights = norm_target_df.iloc[0, :].values
else:
    print(idx)
    print(norm_target_df.index)
    weights = norm_target_df.loc[idx, :].values

penalties, average_j = netphix.comp_penalties(weights, num_samples, penalty)


#################
# run ILP with the original solution


# Create a new model and populate it below. (no network constraints)
model = netphix.create_ILP_model(k, num_samples, num_alterations, num_genes, weights, penalties, mut_lists, alt_df.index,
                         mutated_gene_sets)

# add network constraints
model = netphix.add_density_constraints(model, num_genes, edge_lists, mut_lists, k, density)


ILP_Solutions = []
for itr in range(num_itr):
    print("start solving ILP: itr " + str(itr))
    # Solve
    try:
        model.solve()
    except CplexSolverError as e:
        print("Exception raised during solve: ")
        print(e)
    else:
        solution = model.solution
        selected_muts, TotCost, selected_idx, selected_values = netphix.proc_solution(solution, alt_df, k)
        timer_end = time.clock() - timer_start

        solution_dic = {}
        solution_dic["selected_muts"] = selected_muts
        solution_dic["TotCost"] = TotCost
        solution_dic["selected_values"] = selected_values
        solution_dic["time"] = timer_end
        if itr == 0:  # optimal solution
            OptCost = TotCost
        ILP_Solutions.append((solution_dic))

        if len(selected_muts) > k:  # condition violated --> rerun
            itr -= 1
            continue
        if len(selected_muts) == 0:  # if no modules are selected
            break
        # remove the selected nodes, and find the next modules
        selected_nodes_constraint = cplex.SparsePair(ind=["x" + str(i) for i in selected_idx],
                                                     val=[1] * int(len(selected_idx)))
        model.linear_constraints.add(lin_expr=[selected_nodes_constraint], senses=["E"], rhs=[0])

#################
# permutation test
permuted_weights = [x for x in weights]  # deep copy
PermTotCosts = []

# append results if the file exists
if os.path.isfile(permutefile_name):
    # read file and extract TotCosts
    PermTotCosts += perm.read_permute_file(permutefile_name)
    pstart = len(PermTotCosts)
else: # new file
    pstart = 0
    if permutations > 0:
        Label_list = ["selected_muts", "TotCost", "time", "selected_values"]
        params = [target_df.index[0], k, filter_high, filter_low, density, penalty, norm]
        netphix.write_label(permutefile_name, params, Label_list)

for permute in range(pstart, permutations):
    # adaptive test: STOP condition
    pv = (len(list(filter(lambda x: x >= OptCost, PermTotCosts))) + 1) / (len(PermTotCosts) + 1)
    print(str(permute) + "-th iteration pv is " + str(pv))
    if (permute >= min_permute) & (pv > min_pv):  # stop if p-values not significant
        break

    # permute the target
    random.shuffle(permuted_weights)

    penalties, average_j = netphix.comp_penalties(permuted_weights, num_samples, penalty)

    # Create a new model and populate it below. (no network constraints)
    model = netphix.create_ILP_model(k, num_samples, num_alterations, num_genes, permuted_weights, penalties, mut_lists,
                                      alt_df.index, mutated_gene_sets)

    # add network constraints
    model = netphix.add_density_constraints(model, num_genes, edge_lists, mut_lists, k, density)

    # Solve
    try:
        model.solve()
    except CplexSolverError as e:
        print("Exception raised during solve: ")
        print(e)
    else:
        solution = model.solution

        # Display solution.
        selected_muts, TotCost, selected_idx, selected_values = netphix.proc_solution(solution, alt_df, k)
        timer_end = time.clock() - timer_start
        PermTotCosts.append(TotCost)

        # write permutation solution (write solutions as computed)
        solution_dic = {}
        solution_dic["selected_muts"] = selected_muts
        solution_dic["TotCost"] = TotCost
        solution_dic["selected_values"] = selected_values
        solution_dic["time"] = timer_end
        netphix.write_solutionline(permutefile_name, solution_dic)


# write optimal solution with pv
Label_list = ["selected_muts", "TotCost", "time", "selected_values", "pv"]


if solutionfile_name is not None:
    if idx is None:
        params = [target_df.index[0], k, filter_high, filter_low, density, penalty, norm]
    else:
        params = [idx, k, filter_high, filter_low, density, penalty, norm]
    netphix.write_label(solutionfile_name, params, Label_list)
    for solution_dic in ILP_Solutions:
        if len(PermTotCosts) > 0:
            pv = (len(list(filter(lambda x: x >= OptCost, PermTotCosts))) + 1) / (len(PermTotCosts) + 1)
            print(pv)
            solution_dic["pv"] = pv
        netphix.write_solutionline(solutionfile_name, solution_dic)


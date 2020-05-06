#!/usr/bin/python
# ---------------------------------------------------------------------------

"""
Run NETPHIX. Epirical p-values are computed if target and network permutation files are provided.

This is a simplified version where positive/negative associations are computed separately.
"""

import cplex
from cplex.exceptions import CplexSolverError
import pandas as pd
import time
import sys
import networkx as nx
import os
import argparse
import numpy as np

import netphix_utils as netphix
import permute_utils as perm


timer_start = time.time()

# fixed params
num_itr = 1                             # NEPHLIX finds multiple modules iteratively if num_itr > 1
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
parser.add_argument("net_file", type=str, help="network file", default= "data/HumanStringNet.txt")
parser.add_argument("solutionfile_name", help="file name to write the solution")
# optional arguments
parser.add_argument("-tp", "--target_perm_file",
                    help="file with target permutation results")
parser.add_argument("-np", "--net_perm_file",
                    help="file with net permutation results")
parser.add_argument("-ap", "--alt_perm_file",
                    help="file with alteration table permutation results")
parser.add_argument("-map", "--gene_map_file",
                    help="file with ENSG ID to gene name mapping")
parser.add_argument("-idx", "--index_name",
                    help="index name to use in the target file")
parser.add_argument("-count_file", "--mut_count_file",
                    help="filename with mutation count info")
parser.add_argument("-mut_th", "--hypermutated_th",type=int,
                    help="mutation count threshold to remove hypermutated samples")
parser.add_argument("-subtype_file", "--subtype_file",
                    help="filename with subtype info")
parser.add_argument("-subtype", "--subtype",
                    help="brca subtype -- 'LumA', 'LumB', 'Basal', or 'Her2'")
parser.add_argument("-sep", "--sep",
                    help="separator in file (\\t for default)")
parser.add_argument("-restricted", "--restricted_gene_file",
                    help="file containing restricted gene modules. compute the objective only with genes in the restricted module")
parser.add_argument("-add_file", "--add_file",
                    help="additional mutation file. take OR with the original mutations")
# parser.add_argument("-brca_method", "--brca_method",
#                     help="brca inactivation merging method (BRCA or OR)")
parser.add_argument('--append', action='store_true', help="add solution to existing file")
parser.add_argument('--recompute', action='store_true', help="recompute solution and write to the existing file")

args = parser.parse_args()


alteration_file = args.alt_file
target_file = args.target_file
correlation = args.correlation
k = args.k
filter_high = args.filter_high
filter_low = args.filter_low
density = args.density
norm = args.norm
net_file = args.net_file
solutionfile_name = args.solutionfile_name
target_perm_file = args.target_perm_file
net_perm_file = args.net_perm_file
alt_perm_file = args.alt_perm_file
gene_map_file = args.gene_map_file
idx = args.index_name
subtype=args.subtype
hypermutated_th = args.hypermutated_th
append = args.append
recompute = args.recompute
restricted_gene_file = args.restricted_gene_file
add_file = args.add_file
# brca_method = args.brca_method

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
    sep="\t"

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


if restricted_gene_file is not None:
    restricted_module = open(restricted_gene_file).readlines()[3].split()[1].split(",")
    alt_df = alt_df[alt_df.index.isin(restricted_module)]

#################
# processing data
sys.stdout.write("processing data.. \n")

print(target_df.shape)

if hypermutated_th is not None:
    print("threshold is "+str(hypermutated_th))
    target_df = netphix.remove_hypermutated(target_df, mut_count_file, hypermutated_th)
if subtype is not None:
    target_df = netphix.filter_subtype(target_df, subtype_file, subtype)

print(target_df.shape)

target_df, alt_df, samples, num_samples  = netphix.preproc_data(target_df, alt_df, filter_high, filter_low)

print(target_df.shape)
# exit if there are genes appearing multiple times
# duplicates = alt_df.index[alt_df.index.duplicated()]
# if len(duplicates) > 0:
# todo: check
alt_df = alt_df.loc[~alt_df.index.duplicated()]
    # sys.exit("multiple entries for "+",".join(duplicates))


# normalize
norm_target_df = netphix.norm_target(target_df, norm, correlation)
run_ILP = True
# read solution file if already exists (no need to recompute)
if (append is False) and os.path.isfile(solutionfile_name) and (recompute is False):
        Solution_dics, OptCost, pv = netphix.read_solutionfile(solutionfile_name)
        if len(Solution_dics) > 0:
            run_ILP = False

if run_ILP is True:
    #################
    # prepare to construct ILP model
    sys.stdout.write("preparing to create cplex model.. \n" )

    mutated_gene_sets, num_alterations =  netphix.proc_alt(alt_df)
    edge_lists, mut_lists, num_genes = netphix.proc_net(alt_df.index, all_net, num_alterations)
    if idx is None:
        weights = norm_target_df.iloc[0,:].values
    else:
        weights = norm_target_df.loc[idx, :].values
    penalties, average_j = netphix.comp_penalties(weights, num_samples, penalty)

    #################
    # create and run ILP
    # Create a new model and populate it below. (no network constraints)
    model = netphix.create_ILP_model(k, num_samples, num_alterations, num_genes, weights, penalties, mut_lists, alt_df.index,
                             mutated_gene_sets)

    # add network constraints
    model = netphix.add_density_constraints(model, num_genes, edge_lists, mut_lists, k, density)

    Solution_dics = []
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
            if itr == 0:  # optimal solution
                OptCost = TotCost
            solution_dic = {}
            solution_dic["selected_muts"] = selected_muts
            solution_dic["TotCost"] = TotCost
            solution_dic["selected_values"] = selected_values
            solution_dic["time"] = timer_end
            Solution_dics.append(solution_dic)

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
# compute p-values and write the results

Label_list = ["selected_muts", "TotCost", "time", "selected_values"]
PermTotCosts = []
if (target_perm_file is not None) and os.path.isfile(target_perm_file):
    # read file and extract TotCosts
    PermTotCosts += perm.read_permute_file(target_perm_file)
    Label_list.append("pv")
print(PermTotCosts)

NetPermTotCosts = []
if (net_perm_file is not None) and os.path.isfile(net_perm_file):
    # read file and extract TotCosts
    NetPermTotCosts += perm.read_permute_file(net_perm_file)
    Label_list.append("net_pv")

AltPermTotCosts = []
if (alt_perm_file is not None) and os.path.isfile(alt_perm_file):
    # read file and extract TotCosts
    AltPermTotCosts += perm.read_permute_file(alt_perm_file)
    Label_list.append("alt_pv")

print(Label_list)
# open file
if append: # if append, no label
    file = open(solutionfile_name, 'a' )
else: # else write the parameters and labels
    if idx is None:
        params = [target_df.index[0], k, filter_high, filter_low, density, penalty, norm]
    else:
        params = [idx, k, filter_high, filter_low, density, penalty, norm]
    netphix.write_label(solutionfile_name, params, Label_list)

# write the solution
for solution_dic in Solution_dics:
    TotCost = solution_dic["TotCost"]
    print(TotCost)
    if len(PermTotCosts) > 0:
        pv = (len(list(filter(lambda x: x >= TotCost, PermTotCosts)))+1)/(len(PermTotCosts)+1)
        solution_dic["pv"] = pv
    if len(NetPermTotCosts) > 0:
        net_pv = (len(list(filter(lambda x: x >= TotCost, NetPermTotCosts))) + 1) / (len(NetPermTotCosts) + 1)
        solution_dic["net_pv"]  = net_pv
    if len(AltPermTotCosts) > 0:
        alt_pv = (len(list(filter(lambda x: x >= TotCost, AltPermTotCosts))) + 1) / (len(AltPermTotCosts) + 1)
        solution_dic["alt_pv"] = alt_pv
    print(solution_dic)
    netphix.write_solutionline(solutionfile_name, solution_dic)





import cplex
import sys
import random
import pandas as pd
import scipy.stats as ss
import numpy as np
import itertools

# fixed params
max_size = 100  # maximum module size


def preproc_data(target_df, alt_df, filter_high=2, filter_low=-1):
    """
    remove NA from target profile
    keep only common samples in target and alterations
    filter genes based on mutation frequency

    :param target_df: target profile
    :param alt_df: alteration matrix
    :param filter_high: mutation frequency upperbound (excluding the val)
    :param filter_low: mutation frequency lowerbound (excluding the val)
    :return: target_df, alt_df, samples
    """
    # drop NA
    target_df.dropna(axis=1, inplace=True)

    # keep common samples only
    samples = list(set(target_df.columns).intersection(alt_df.columns))
    num_samples = len(samples)
    alt_df = alt_df.loc[:, samples]
    target_df = target_df.loc[:, samples]

    # filtering with mut frequency
    column_sum = alt_df.sum(axis=1)
    Count_Col = alt_df.shape[1]
    up = Count_Col * filter_high
    down = Count_Col * filter_low
    alt_df = alt_df[((column_sum < up) & (column_sum > down))]

    return target_df, alt_df, samples, num_samples

def remove_hypermutated(target_df, count_file, mut_count_th):
    """
    :param target_df:
    :param count_file: "../data/count_mutations_per_cluster_sig_per_sample_total_mut_count.tsv"
    :param mut_count_th: remove samples with mutations >= mut_count_th
    :return:
    """
    count = pd.read_csv(count_file, sep="\t", index_col=0)
    filtered_samples = count.index[count["count"] < mut_count_th]
    filtered_samples = list(set(target_df.columns).intersection(filtered_samples))
    return target_df[filtered_samples]

def filter_subtype(target_df, subtype_file, subtype):
    """
    :param target_df:
    :param subtype_file: "../data/BRCA-EU.Expression.Subtyping.txt"
    :param subtype: 'LumA', 'LumB', 'Basal', 'Her2' or 'Lum' ('LumA'+'LumB')
    :return:
    """
    print(subtype_file)
    subtype_df = pd.read_csv(subtype_file, sep="\t", index_col=0)
    if subtype == "Lum":
        filtered = (subtype_df["AIMS.subtype"] == "LumA")|(subtype_df["AIMS.subtype"] == "LumB")
        filtered_samples = subtype_df.index[filtered]
    else:
        filtered_samples = subtype_df.index[subtype_df["AIMS.subtype"] == subtype]
    filtered_samples = list(set(target_df.columns).intersection(filtered_samples))
    return target_df[filtered_samples]

def norm_target(target_df, norm, correlation):
    """
    normalize target profile and take negative if correlation is negative
    :param target_df: target profile
    :param norm:  normalization method z (z score) or zlog (zscore of log10)
    :param correlation: negative or positive
    :return: norm_target_df
    """
    # normalize
    norm_target_df = pd.DataFrame()
    if norm == "z":
        norm_target_df = target_df.apply(ss.zscore, axis=1, result_type='broadcast')
    elif norm == "zlog":
        norm_target_df = np.log10(target_df + 1).apply(ss.zscore, axis=1, result_type='broadcast')

    if correlation == 'negative':
        norm_target_df = -norm_target_df

    return norm_target_df


def proc_alt(alt_df):
    """
    create list each way (sample-> mutated genes, gene-> mutated samples)
    :param alt_df: alteration matrix
    :return: mutated_gene_sets, num_alterations
    """
    mutated_sample_sets = []  # mutated samples for each gene
    for i in range(alt_df.shape[0]):
        row = alt_df.iloc[i].values
        mutated_sample_sets.append(np.where(row == 1)[0])

    mutated_gene_sets = []  # mutated genes for each sample
    for i in range(alt_df.shape[1]):
        col = alt_df.iloc[:, i].values
        mutated_gene_sets.append(np.where(col == 1)[0])

    num_alterations = len(mutated_sample_sets)
    return mutated_gene_sets, num_alterations


def proc_net(mutated_genes, all_net, num_alterations):
    """
    create
    1) edge_lists (for each gene, list of neighbors)
    2) mut_list (for each gene, list of mutation idxs -- tp53_del, tp53_mut etc)
    :param mutated_genes: list of genes with mutation types
    :param all_net: network
    :param num_alterations: # alterations
    :return: edge_lists, mut_lists
    """
    # Create edge list
    genes = list(set([g.split("_")[0] for g in mutated_genes]))
    net = all_net.subgraph(genes)
    num_genes = len(genes)
    sys.stdout.write("network size: %d nodes and %d edges\n" % (len(net), len(net.edges())))
    gene_idx = dict([(genes[j], j) for j in range(num_genes)])
    edge_lists = {}
    for j in range(num_genes):
        if genes[j] not in net:
            edge_lists[j] = []
            continue
        edge_lists[j] = [gene_idx[neigh] for neigh in net.neighbors(genes[j])]

    # create a mut_gene list dict = {(gene_idx, idxs of mutations for the gene)}
    mut_lists = {}
    for j in range(num_genes):
        mut_lists[j] = []
    for i in range(num_alterations):
        mutated_gene = mutated_genes[i].split("_")[0]
        if mutated_gene not in gene_idx:
            continue
        mut_lists[gene_idx[mutated_gene]].append(i)
    return edge_lists, mut_lists, num_genes


def comp_penalties(weights, num_samples, penalty):
    """
    compute weights/penalties in the objective function
    :param target_df: target profile
    :param num_samples: number of samples
    :param penalty: p or np or value given for average_j
    :return weights, penalties
    """


    if penalty == "p":
        sys.stdout.write("\timpose penalty..")
        # average_i = target_df.iloc[0][target_df.iloc[0] > 0].mean()
        average_j = np.average(list(filter(lambda w: w > 0, weights)))
    elif penalty == "np":
        sys.stdout.write("\tno penalty..")
        average_j = 0  # no penalty

    penalties = []
    for l in range(0, num_samples):
        if weights[l] > 0:
            penalties.append(average_j)
        else:
            penalties.append(-weights[l])
    return penalties, average_j


def create_ILP_model(k, num_samples, num_alterations, num_genes, weights, penalties, mut_lists, mutated_genes, mutated_gene_sets):
    """
    create ILP and populate (no network information)

    :param k: size of module
    :param num_samples: number of patients
    :param num_alterations: number of alterations
    :param num_genes: number of genes (can be different from num_alterations when a gene has different types of alterations
    :param weights: phenotype
    :param penalties: penalty for mutual exclusivity
    :param mut_lists
    :param mutated_genes: a mut_gene list dict = {(gene_idx, idxs of mutations for the gene)}
    :param mutated_gene_sets: mutated genes for each sample
    :return: ILP model

    """
    # Create a new (empty) model and populate it below.
    model = cplex.Cplex()
    # Set up cplex parameters to use
    model.parameters.threads.set(4)
    model.parameters.mip.limits.auxrootthreads.set(4)

    # Add variables
    model.variables.add(names=["z" + str(j) for j in range(num_samples)],
                        obj=[weights[j] + penalties[j] for j in range(num_samples)], lb=[0] * num_samples,
                        ub=[1] * num_samples,
                        types=["B"] * num_samples)

    model.variables.add(names=["y" + str(j) for j in range(num_samples)],
                        obj=[-penalties[j] for j in range(num_samples)], lb=[0] * num_samples,
                        ub=[100] * num_samples,
                        types=["I"] * num_samples)

    model.variables.add(names=["x" + str(i) for i in range(num_alterations)], lb=[0] * num_alterations,
                        ub=[1] * num_alterations,
                        types=["B"] * num_alterations)

    # for mutation/gene mapping
    model.variables.add(names=["g" + str(i) for i in range(num_genes)], lb=[0] * num_genes,
                        ub=[1] * num_genes,
                        types=["B"] * num_genes)

    # module size
    model.variables.add(names=["m"], lb=[0],
                        ub=[max_size],
                        types=["I"])

    # Set the type of each variables
    # for i in range(num_alterations):
    # 		model.variables.set_types("x"+str(i), model.variables.type.binary)
    for j in range(num_samples):
        model.variables.set_types("y" + str(j), model.variables.type.integer)

    # Add constraints!
    sets_constraint = cplex.SparsePair(ind=["x" + str(i) for i in range(num_alterations)],
                                       val=[1.0] * num_alterations)
    # size no more than k
    model.linear_constraints.add(lin_expr=[sets_constraint], senses=["L"], rhs=[k])

    for j in range(num_samples):
        number_constraint = cplex.SparsePair(ind=["y" + str(j), "z" + str(j)], val=[1.0, -1.0])
        model.linear_constraints.add(lin_expr=[number_constraint], senses=["G"], rhs=[0])

    for j in range(num_samples):
        number2_constraint = cplex.SparsePair(ind=["y" + str(j), "z" + str(j)], val=[1.0, -k])
        model.linear_constraints.add(lin_expr=[number2_constraint], senses=["L"], rhs=[0])

    # m = sum(xi)
    index = ["m"]
    value = [1]
    for j in range(num_alterations):
        index.append("x" + str(j))
        value.append(-1)
    size_constraint = cplex.SparsePair(index, value)
    model.linear_constraints.add(lin_expr=[size_constraint], senses=["E"], rhs=[0])

    # loss/LOH and gain cannot be selected at the same time
    for i in range(num_genes):
        gain, loss, loh, am, de = -1, -1, -1, -1, -1
        for j in mut_lists[i]:
            if mutated_genes[j].endswith("gain"):
                gain = j
            elif mutated_genes[j].endswith("loss"):
                loss = j
            elif mutated_genes[j].endswith("LOH"):
                loh = j
            elif mutated_genes[j].endswith("amp"):
                am = j
            elif mutated_genes[j].endswith("del"):
                de = j
        if (gain != -1) and (loss != -1):
            cnv_constraint = cplex.SparsePair(ind=["x" + str(gain), "x" + str(loss)], val=[1.0, 1.0])
            model.linear_constraints.add(lin_expr=[cnv_constraint], senses=["L"], rhs=[1])
        if (gain != -1) and (loh != -1):
            cnv_constraint = cplex.SparsePair(ind=["x" + str(gain), "x" + str(loh)], val=[1.0, 1.0])
            model.linear_constraints.add(lin_expr=[cnv_constraint], senses=["L"], rhs=[1])
        if (am != -1) and (de != -1):
            cnv_constraint = cplex.SparsePair(ind=["x" + str(am), "x" + str(de)], val=[1.0, 1.0])
            model.linear_constraints.add(lin_expr=[cnv_constraint], senses=["L"], rhs=[1])
    sys.stdout.write("\nadding xy constraints..")
    for j in range(num_samples):
        if j % 100 == 0:
            sys.stdout.write("%d\t" %j)
        index = ["y" + str(j)]
        value = [1.0]
        for i in mutated_gene_sets[j]:
            index.append("x" + str(i))
            value.append(-1.0)
        number3_constraint = cplex.SparsePair(ind=index, val=value)
        model.linear_constraints.add(lin_expr=[number3_constraint], senses=["E"], rhs=[0])

    # Our objective is to minimize cost. Fixed and variable costs
    # have been set when variables were created.
    model.objective.set_sense(model.objective.sense.maximize)
    return model


def add_density_constraints(model, num_genes, edge_lists, mut_lists, k, density):
    """
    add density constraints to existing ILP model
    :param model: existing ILP model
    :param num_genes:
    :param edge_lists:
    :param mut_lists:
    :param k: module size
    :param density: density of a module (connectivity)
    :return: model (updated model)
    """

    # mutation/gene mapping
    # make sure gi is 1 iff one or more xj is 1
    # gi <= sum_{j in mlist(i)}(xj)
    for i in range(num_genes):
        index = ["g" + str(i)]
        value = [1]
        for i1 in mut_lists[i]:
            index.append("x" + str(i1))
            value.append(-1)
        density_constraint = cplex.SparsePair(ind=index, val=value)
        rhs = 0
        model.linear_constraints.add(lin_expr=[density_constraint], senses=["L"], rhs=[rhs])
    # mutation/gene mapping gi >= (xj) {all j in mlist(i)}
    for i in range(num_genes):
        for i1 in mut_lists[i]:
            index = ["g" + str(i)]
            value = [1]
            index.append("x" + str(i1))
            value.append(-1)
            density_constraint = cplex.SparsePair(ind=index, val=value)
            rhs = 0
            model.linear_constraints.add(lin_expr=[density_constraint], senses=["G"], rhs=[rhs])

    # density constraints (size less than k)
    sum = (k - 1.0) * density
    count = 0
    sys.stdout.write("\nadding density contraints..")
    for i in range(num_genes):
        if i % 1000 == 0:
            sys.stdout.write("%d\t" % i)
        index = ["g" + str(i)]
        value = [-sum]
        for i1 in edge_lists[i]:
            index.append("g" + str(i1))
            value.append(1.0)
        index.append("m")
        value.append(-density)
        density_constraint = cplex.SparsePair(ind=index, val=value)
        rhs = -(sum + density)
        model.linear_constraints.add(lin_expr=[density_constraint], senses=["G"], rhs=[rhs])
    return model


def run_bootstrap(alt_df, target_df, num_random):
    """
    perform bootstrap
    :param alt_df: alteration matrix, columns are samples and rows are genes
    :param target_df: target file, columns are samples
    :param num_random: how many to select?
    :return: new_alt_df, new_target_df
    """
    orig_sample_size = target_df.shape[1]
    bootstrap = [random.choice(range(orig_sample_size)) for i in range(num_random)]
    return alt_df.iloc[:, bootstrap], target_df.iloc[:,bootstrap]


def proc_solution(solution, alt_df, k):
    """
    given ILP solution, extract necessary information
    :param solution: ILP solution
    :param alt_df: alteration matrix, columns are samples and rows are genes
    :param k: module size
    :return: selected_muts, TotCost, selected_idx, selected_values
    """

    # extract information from solution
    num_alterations = alt_df.shape[0]
    TotCost = solution.get_objective_value()
    m = solution.get_values("m")
    selected_muts = []
    selected_values = []
    selected_idx = []  # for iteration
    for i in range(num_alterations):
        if solution.get_values("x" + str(i)) > 0.5:
            selected_values.append(solution.get_values("x" + str(i)))
            mut_name = alt_df.index[i]
            selected_muts.append(mut_name)
            selected_idx.append(i)
    # display solution
    print("Solution status = ", solution.get_status(), ":", end=' ')
    print(solution.status[solution.get_status()])
    print("Total cost = ", TotCost)
    print(",".join(selected_muts))
    print("size of module is " + str(m))
    if len(selected_muts) > k:  # condition violated?
        print(selected_values)
    return selected_muts, TotCost, selected_idx, selected_values


def sum_bootstrap_results(solution_dics):
    """
    count # appearance of each gene_mutations in bootstrapping
    :param solutions: list of solution dics from bootsrapping
    :return: count_dic
    """
    all_genes = []
    all_edges = []
    for sdic in solution_dics:
        selected_muts  = sdic["selected_muts"]
        all_genes.extend(selected_muts)
        selected_muts.sort()  # order doesn't matter
        all_edges.extend(list(itertools.combinations(selected_muts, 2)))
    gene_set = set(all_genes)
    gene_count_dic = {}
    for gene in gene_set:
        gene_count_dic[gene] = (all_genes.count(gene))
    edge_set = set(all_edges)
    edge_count_dic = {}
    for edge in edge_set:
        edge_count_dic[edge] = (all_edges.count(edge))
    return gene_count_dic, edge_count_dic


def read_solutionfile(ILP_file):
    """
    read ILP file and return solutions
    file format:
        required: selected_muts, TotCost, time
        optional: selected_values, pv, net_pv, alt_pv

    :param ILP_file:
    :return: Solutions: list of solution dic
    :return: OptCost: optimal cost
    :return opt_pv: pv (None if not given)

    """
    lines = open(ILP_file).readlines()
    # sys.stderr.write("%s" % lines[0])  # target & params
    labels = lines[1].split()

    Solutions = []
    # solution = selected_muts, selected_values, TotCost, pv, net_pv, timer_end
    OptCost=-1
    opt_pv = -1
    for i in range(2, len(lines)):
        solution_dic = {}
        tkns = tuple(lines[i].split())
        idx = labels.index("selected_muts")
        solution_dic["selected_muts"] = tkns[idx].split(",")
        idx = labels.index("TotCost")
        solution_dic["TotCost"] = float(tkns[idx])
        idx = labels.index("time")
        solution_dic["time"] = float(tkns[idx])
        if "selected_values" in labels:
            idx = labels.index("selected_values")
            solution_dic["selected_values"] = tkns[idx].split(",")
        if "pv" in labels:
            idx = labels.index("pv")
            solution_dic["pv"] = float(tkns[idx])
        if "net_pv" in labels:
            idx = labels.index("net_pv")
            solution_dic["net_pv"] = float(tkns[idx])
        if "alt_pv" in labels:
            idx = labels.index("alt_pv")
            solution_dic["alt_pv"] = float(tkns[idx])
        Solutions.append(solution_dic)
        if i == 2:
            OptCost = solution_dic["TotCost"]
            if "pv" in labels:
                opt_pv = solution_dic["pv"]

    return Solutions, OptCost, opt_pv


def write_label(filename, params, label_list):
    """
    write parameter in the first row
    label in the second
    """
    file = open(filename, 'w')
    file.write("%s\n" %"\t".join([str(x) for x in params]))
    file.write("\t%s\n" % "\t".join(label_list))
    file.close()


def write_solutionline(filename, sdic):
    """
    write a solution in one line
    required: selected_muts, TotCost, time
    optional: selected_values, pv, net_pv, alt_pv
    :param filename:
    :param sdic: solution dic
    :return:
    """
    file = open(filename, 'a')
    file.write('\t%s' % ",".join(sdic["selected_muts"]))
    file.write('\t%f' % sdic["TotCost"])
    file.write('\t%f' % (sdic["time"]))
    if "selected_values" in sdic:
        file.write('\t%s' % ",".join([str(v) for v in sdic["selected_values"]]))
    if "pv" in sdic:
        file.write('\t%f' % sdic["pv"])
    if "net_pv" in sdic:
        file.write('\t%f' % sdic["net_pv"])
    if "alt_pv" in sdic:
        file.write('\t%f' % sdic["alt_pv"])
    file.write("\n")
    print("writing solution..")
    file.close()






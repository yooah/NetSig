# NETPHIX

NETPHIX (NETwork-to-PHenotype association wIth eXclusivity) is a computational tool
to identify mutated subnetworks that are associated with a continuous phenotype. 
Using an integer linear program with properties for mutual exclusivity and interactions among genes,
NETPHIX finds an optimal set of genes maximizing the association. For more details, please see [1].


### Requirements

* Linux/Mac OS/Windows
* python version 3.6.2/pandas version 0.24.2
* CPLEX version 12.9.0.0 (set PYTHONPATH to the cplex location);
    - install the CPLEX-Python modules using the script "setup.py" located in "yourCplexhome/python/"
        $ python setup.py install
    - install python CPLEX interface
        $ pip install cplex

### How to use it

#### run_NETPHIX.py

Run NETPHIX. Epirical p-values are computed if target and network permutation files are provided.
This is a simplified version where positive/negative associations are computed separately.

usage: run_NETPHIX.py [-h] [-tp TARGET_PERM_FILE] [-np NET_PERM_FILE]
                      [-ap ALT_PERM_FILE] [-map GENE_MAP_FILE]
                      [-idx INDEX_NAME] [-count_file MUT_COUNT_FILE]
                      [-mut_th HYPERMUTATED_TH] [-subtype_file SUBTYPE_FILE]
                      [-subtype SUBTYPE] [-sep SEP]
                      [-restricted RESTRICTED_GENE_FILE] [-add_file ADD_FILE]
                      [--append] [--recompute]
                      alt_file target_file correlation k filter_high
                      filter_low density norm net_file solutionfile_name

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
  net_file              network file
  solutionfile_name     file name to write the solution

optional arguments:
  -h, --help            show this help message and exit
  -tp TARGET_PERM_FILE, --target_perm_file TARGET_PERM_FILE
                        file with target permutation results
  -np NET_PERM_FILE, --net_perm_file NET_PERM_FILE
                        file with net permutation results
  -ap ALT_PERM_FILE, --alt_perm_file ALT_PERM_FILE
                        file with alteration table permutation results
  -map GENE_MAP_FILE, --gene_map_file GENE_MAP_FILE
                        file with ENSG ID to gene name mapping
  -idx INDEX_NAME, --index_name INDEX_NAME
                        index name to use in the target file
  -count_file MUT_COUNT_FILE, --mut_count_file MUT_COUNT_FILE
                        filename with mutation count info
  -mut_th HYPERMUTATED_TH, --hypermutated_th HYPERMUTATED_TH
                        mutation count threshold to remove hypermutated
                        samples
  -subtype_file SUBTYPE_FILE, --subtype_file SUBTYPE_FILE
                        filename with subtype info
  -subtype SUBTYPE, --subtype SUBTYPE
                        brca subtype -- 'LumA', 'LumB', 'Basal', or 'Her2'
  -sep SEP, --sep SEP   separator in file (\t for default)
  -restricted RESTRICTED_GENE_FILE, --restricted_gene_file RESTRICTED_GENE_FILE
                        file containing restricted gene modules. compute the
                        objective only with genes in the restricted module
  -add_file ADD_FILE, --add_file ADD_FILE
                        additional mutation file. take OR with the original
                        mutations
  --append              add solution to existing file
  --recompute           recompute solution and write to the existing file



Examples: 

- Run NETPHIX with Mutational Signature 2C
(k=5, mutation frequency > 0.01, network density 0.5, normalization method=zlog)

    $ python run_NETPHIX.py data/Sig_alt_data/Alterations2_1.txt data/Sig_count/Target2_1.txt positive 5 1 0.01 0.5 zlog data/HumanStringNet.txt results/netphix_results2_1.txt

Run NETPHIX with Mutational Signature 2C with the same parameters when target permutation results are available
    
    $ python run_NETPHIX.py data/Sig_alt_data/Alterations2_1.txt data/Sig_count/Target2_1.txt positive 5 1 0.01 0.5 zlog data/HumanStringNet.txt results/netphix_results2_1.txt -tp results/netphix_perm2_1.txt

### target_permute_NETPHIX.py

Run target (phenotype) permutation. The phenotype is permuted across samples.

Usage:

    target_permute_NETPHIX.py [-h]
                                  alt_file target_file correlation k
                                  filter_high filter_low density norm
                                  permutations permutefile_name

Arguments:

    alt_file          alteration data file (alteration matrix, see data/Alterataion_data/ directory for an example)
    target_file       target profile data file (target (phenotype) see data/Target_data/ directory for an example)
    correlation       positive or negative association
    k                 size of modules
    filter_high       upperbound mutation frequency of a gene to be included
    filter_low        lowerbound mutation frequency of a gene to be included
    density           module density
    norm              target normalization method, z or zlog
    permutations      num of permutations
    permutefile_name  file name to write the permutation results


Examples: 

Run target permutation 20 times with Mutational Signature 2C
(k=5, mutation frequency > 0.01, network density 0.5, normalization method=zlog)

    $ python target_permute_NETPHIX.py data/Sig_alt_data/Alterations2_1.txt data/Sig_count/Target2_1.txt positive 5 1 0.01 0.5 zlog data/HumanStringNet.txt 20 results/netphix_perm2_1.txt
    


## References
[1] Identifying Drug Sensitivity Subnetworks with NETPHIX. https://www.biorxiv.org/content/10.1101/543876v1?rss=1
## Authors

* **Yoo-Ah Kim** (kimy3@ncbi.nlm.nih.gov) 
* **Rebecca Sarto Basso** (rebeccasarto@berkeley.edu)





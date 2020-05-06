# Network Based Analysis for Mutational Signatures

To investigate the genetic aberrations associated with mutational signatures, we took network-based approach aiming to answer the following two complementary questions: 
	

 1.  what are functional pathways whose gene expression activities   
    correlate with the strengths of mutational signatures? 	
 2. are there pathways whose genetic alterations might have led to specific 
    mutational signatures?

Please see [1] for the details of the methods and the analysis results.

## Expression Correlation Modules

#### Jupyter notebook files for analysis
    expr_analysis:
    comp_corr_sigma.ipynb 
    	--> read expression and mutation count data files and compute correlation
    corr_analysis.ipynb  
    	--> filtering and clustering genes based correlation patterns
#### data files   
    expr_analysis/data:
    c5.bp.v6.1.symbols.gmt  
    	--> GO terms downloaded from https://www.gsea-msigdb.org/gsea/downloads.jsp
    	    c5: gene ontology (GO) gene sets,version 6.1
    common_brca_expr_pid.tsv  
    	--> gene expression
    common_sigma_count.tsv
    	--> mutation count for each signature
    
#### results files
      
    expr_analysis/results:
    all_cluster_id.tsv  
    	--> clustering ID for selected genes 
    meta_cluster_id.tsv
    	--> clustering ID for selected genes from DNA metabolic genes

## Mutation modules with NetPhix

#### NetPhix source files (see NETPHIX_README.md in the directory)

    netphix:
    run_NETPHIX.py       
    	--> run NETPHIX
    target_permute_NETPHIX.py
    	--> run NETPHIX on permuted instances
    select_sig_modules.ipynb
    	--> select the best module for each signature
    permute_utils.py  
    netphix_utils.py  
    	--> utility functions
    
#### data files for NetPhix analysis
    netphix/data:
    HumanStringNet.txt  
    	--> functional interaction network downloaded from string-db.org
    	    high confidence edges only (900/1000)
    biallelic_inactivation_genes_per_sample.csv.gz
    	--> biallelic inactivation data 
    	    
    
    netphix/data/Sig_alt_data:
    Alterations$sig_id$_0.txt
    Alterations$sig_id$_1.txt
    	--> gene alteration data for clustered (0) and dispersed (1) signatures
    
    netphix/data/Sig_count:
    Target$sig_id$_0.txt
    Target$sig_id$_1.txt
    	--> mutational signature count data for clustered (0) and dispersed (1) signatures
    
 #### NetPhix analysis results
   
    netphix/results:
    all_sig_results_brca_combined.tsv  
    	--> netphix results for all parameters (biallelic data combined)

## References
[1] Network-based approaches elucidate differences within APOBEC and clock-like signatures inbreast cancer, Genome Medicine 2020. 

## Authors

* **Yoo-Ah Kim** (kimy3@ncbi.nlm.nih.gov) 

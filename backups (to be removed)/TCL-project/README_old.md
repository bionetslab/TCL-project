![scAnalyzer](/readme-images/scAnalyzer-logo.png)

# scAnalyzer
The source code for scAnalyzer (**s**ingle and **c**ell **Analyzer**). It performs clustering and differential analyses on spatial proteomics data.

## Tutorial
### I. Clustering
#### i. Perform hierarchical clustering
- Agglomerative clustering performed using the [schist nested model](https://schist.readthedocs.io/en/latest/clustering_pbmc.html#clustering-pbmc).
- The schist nested model achieves a more fine-grained clustering as compared to Leiden or schist planted model, achieved over several clustering levels--thereby giving the user higher control over cluster analyses, and how many clusters to stop the clustering algorithm at.
###### Running the clustering algorithm:
Open command prompt/ terminal, then run:
```bash
python3 1_compute_list_of_essential_proteins_for_clustering.py <adata_pickle_path> --user_selected_cluster_level <user_selected_cluster_level>
```

The positional arguments are:
```
[1] adata_pickle_path                   Description: Specify path to anndata pickle data; type=str
```

The optional arguments are:
```
[1] --user_selected_cluster_level       Description: Specify cluster level to get essential proteins list from; type=str [read more below]

{any_integer_value, "complete"} accepted, default is 'complete':
    - If value is integer positive: That represents the cluster level
    - If value is integer negative: That represents the cluster number in reverse order (for example: -1 represents (n-1)-th cluster)
    - If value is 'complete': That represents the last (n-th) cluster. Please note that if your input is "complete", just type it in normally without inverted commas/ quotation marks -- the program is capable of reading strings directly without any quotation symbols.
    - For all other values (i. negative integers; ii. fractional values; iii. any strings other than "complete"; iv. integer value is out of range of [1, n] where n is the total number of clusters): the default is n.
```

The outputs are:
- List of essential proteins per cluster per patient
- Dotplot of potential marker genes vs. cell type clusters (example below for **patient id: 1** of the [TNBC MIBI dataset](https://www.science.org/doi/full/10.1126/sciadv.aax5851))
<br/><br/><br/><br/>
![dotplot_____PatientID-1_____ClusteringLevel-2](readme-images/dotplot_____PatientID-1_____ClusteringLevel-2.png)


#### ii. Protein importance:
- In this section, we compute and visualize importances of essential proteins that drive clustering results.

###### Running the protein importance computation:
Open command prompt/ terminal, then run:
```bash
python3 2_visualize_importances_of_essential_proteins_for_clustering.py <adata_pickle_path> <essential_proteins_per_patient_pickle_path>
```

The positional arguments are:
```
[1] adata_pickle_path                                       Description: Specify path to anndata pickle data; type=str
[2] essential_proteins_per_patient_pickle_path              Description: Specify path to essential proteins dictionary generated from the clustering (previous step); type=str
```


The outputs are:
- Feature importance scores calculated using **permutation importance** on **one-vs-all classifier** for **patient id: 1** in the [TNBC MIBI dataset](https://www.science.org/doi/full/10.1126/sciadv.aax5851):
<br/><br/><br/><br/>
![feature_imp_scores_patient1_PermImportance_OneVsAllClassifier](readme-images/feature_imp_scores_patient1_PermImportance_OneVsAllClassifier.png)

- Feature importance scores calculated using **Gini index** on **random forest (RF) classifier** for **patient id: 1** in the [TNBC MIBI dataset](https://www.science.org/doi/full/10.1126/sciadv.aax5851):
<br/><br/><br/><br/>
![feature_imp_scores_patient1_Gini_RFClassifier](readme-images/feature_imp_scores_patient1_Gini_RFClassifier.png)

- Feature importance scores calculated using **permutation importance** on **random forest (RF) classifier** for **patient id: 1** in the [TNBC MIBI dataset](https://www.science.org/doi/full/10.1126/sciadv.aax5851):
<br/><br/><br/><br/>
![feature_imp_scores_patient1_PermImportance_RFClassifier](readme-images/feature_imp_scores_patient1_PermImportance_RFClassifier.png)

### II. Cluster analyses (deciding on the optimal cluster level)
#### i. Maximum bipartite matching between clusterings:
- The similarity scores between clusterings have been obtained for all patients, and their distributions have been plotted against cluster level in the schist nested model.

###### Running the protein importance computation:
Open command prompt/ terminal, then run:
```bash
python3 3_plot_of_robustness_vs_intercluster_similarity.py <adata_pickle_path> --method <method>
```

The positional arguments are:
```
[1] adata_pickle_path                   Description: Specify path to anndata pickle data; type=str
```

The optional arguments are:
```
[1] --method                            Description: Specify cluster agreement method between clusterings; type=str; options={"top-down", "bottom-up"}; default='top-down' [read more below]

Note: The schist hierarchical clustering method is agglomerative, meaning bottom-up.
```

The outputs are:
- Plot of robustness vs. inter-cluster similarity (example below for **patient id: 1** of the [TNBC MIBI dataset](https://www.science.org/doi/full/10.1126/sciadv.aax5851)):
<br/><br/><br/><br/>
(A) Top-down cluster mapping:
<br/><br/><br/><br/>
![Plot_of_robustness_vs_inter_cluster_similarity_method_topDown](readme-images/Plot_of_robustness_vs_inter_cluster_similarity_method_topDown.png)

<br/><br/><br/><br/>
(B) Bottom-up cluster mapping:
<br/><br/><br/><br/>
![Plot_of_robustness_vs_inter_cluster_similarity_method_bottomUp](readme-images/Plot_of_robustness_vs_inter_cluster_similarity_method_bottomUp.png)

<br/><br/><br/><br/>
<br/><br/><br/><br/>

- You can decide on the optimal cluster level based on these plots.
- For our case study on the [TNBC MIBI dataset](https://www.science.org/doi/full/10.1126/sciadv.aax5851) dataset, we will just set optimal cluster level $n = N-1$, where $N$ is the highest cluster level in the schist agglomerative model, where there is just one cluster.





#### ii. Measure no. of clusters per cluster level:

- This may further help some users decide how many cluster levels to set for the schist agglomerative clustering model.

###### Running the computation for no. of clusters per cluster level:
Open command prompt/ terminal, then run:
```bash
python3 4_plot_no_of_cluster_levels_vs_no_of_clusters_per_cluster_level.py <adata_pickle_path> <clusterings_patientLevel_dict_path>
```

The positional arguments are:
```
[1] adata_pickle_path                                                     Description: Specify path to anndata pickle data; type=str
[2] clusterings_patientLevel_dict_path                                    Specify path to clustering combinations calculated in the previous step (saved as 'clusterings_patientLevel_dict.pkl' in step **II.i. Maximum bipartite matching between clusterings**)
```

The outputs are:
- Plot of no. of clusters per cluster level across all samples in the [TNBC MIBI dataset](https://www.science.org/doi/full/10.1126/sciadv.aax5851):
<br/><br/><br/><br/>
![No_of_cluster_levels_vs_No_of_clusters_per_cluster_level](readme-images/No_of_cluster_levels_vs_No_of_clusters_per_cluster_level.png)






### III. Differential analyses (before cell type annotations)
#### i. Protein correlations
- Generate protein correlation matrix.
- Subsequently, perform MWU-test on protein correlation values between conditions.
- Retain correlations with p-values<0.05 as important protein-protein correlations.

###### Obtaining the essential protein correlations:
Open command prompt/ terminal, then run:
```bash
python3 5_DA_protein_correlations.py <adata_pickle_path> <dependent_variable_name>
```

The positional arguments are:
```
[1] adata_pickle_path                   Description: Specify path to anndata pickle data; type=str
[2] dependent_variable_name             Description: Specify the name of the dependent variable; type=str
    
Note: The aforementioned variable <dependent_variable_name> should be present under observation metadata (obsm) of the anndata onject containing gene/ protein expression.
```

The outputs are:
- Bar plot of significant protein correlation p-values
<br/><br/><br/><br/>
![protein_correlation_pvalues_impProteins_allPatients.png](readme-images/protein_correlation_pvalues_impProteins_allPatients.png)



#### ii. Check multiple protein profiles
- Count all cells containing all N-combinaions of proteins present above the threshold.
- Subsequently, perform MWU-test on no. of cells protein profiles expressed in, across conditions.
- Retain protein profiles with p-values<0.05 as significant protein profiles.

###### Obtaining the essential protein correlations:
Open command prompt/ terminal, then run:
```bash
python3 6_DA_multiple_protein_profiles.py <adata_pickle_path> <dependent_variable_name> --N <N> --threshold <threshold>
``` 

The positional arguments are:
```
[1] adata_pickle_path                   Description: Specify path to anndata pickle data; type=str
[2] dependent_variable_name             Description: Specify the name of the dependent variable; type=str
    
Note: The aforementioned variable <dependent_variable_name> should be present under observation metadata (obsm) of the anndata onject containing gene/ protein expression.
```

The optional arguments are:
```
[1] --N                         Description: Specify number of proteins to be analyzed within protein profile; type=int, default=2
[2] --threshold                 Description: Specify minimum value over which protein expressions are considered for further analyses; type=float, default=0.5
```

The outputs are:
- Bar plot of significant, multiple protein coexpression p-values
<br/><br/><br/><br/>
![multi_protein_coexpression_pvalues_impProteins_allPatients.png](readme-images/multi_protein_coexpression_pvalues_impProteins_allPatients.png)
<br/><br/><br/><br/>
Nothing found for N=2, threshold=0.5, across all patients in the [TNBC MIBI dataset](https://www.science.org/doi/full/10.1126/sciadv.aax5851).



#### iii. PCA MWU analysis for finding important proteins (removed)
- First principal component is derived from each of the cells. MWU-test is performed on the original data.
- Each of the proteins are randomly shuffled, and the first step is repeated.
- Provided that the p-value of the actual protein expression between conditions<0.05: Then, for proteins where the p-value decreases on random shuffling are considered unimportant. The proteins where p-value goes up on random shuffling, are considered important.

###### Obtaining the essential protein correlations:
Open command prompt/ terminal, then run:
```bash
python3 7_DA_PCA_MWU_analysis.py <adata_pickle_path> <dependent_variable_name>
```

The positional arguments are:
```
[1] adata_pickle_path                   Description: Specify path to anndata pickle data; type=str
[2] dependent_variable_name             Description: Specify the name of the dependent variable; type=str
    
Note: The aforementioned variable <dependent_variable_name> should be present under observation metadata (obsm) of the anndata onject containing gene/ protein expression.
```

The outputs are:
- Bar plot of p-values of significant proteins (p-value<0.05) [Note: there is a bug in the code; the p-value<0.05 check has not been done.]
<br/><br/><br/><br/>
![DA_significant_protein_expressions_pvalues.png](readme-images/DA_significant_protein_expressions_pvalues.png)


##
##
##

# Under development:

### IV. Differential analyses (after cell type annotations)
#### i. Neighborhood enrichment, MWU-test
#### ii. Co-occurence, MWU-test
#### iii. Centrality scores
#### iv. Ripley's statistics
#### v. Spatial autocorrelation
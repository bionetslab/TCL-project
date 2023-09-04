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
![Plot of robustness vs. inter-cluster similarity method top-down](readme-images/Plot of robustness vs. inter-cluster similarity method top-down.png)





    
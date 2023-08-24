import PyWGCNA
import pandas as pd

# geneExp = '/Users/surya/Documents/GITHUB-REPOSITORIES/spatial_proteomics/2. feature_extraction/5xFAD_paper/Patient_20751_relativeAreaWise_cell_expression_count.csv'
# # geneExp = '/Users/surya/Documents/GITHUB-REPOSITORIES/spatial_proteomics/2. feature_extraction/5xFAD_paper/expressionList.csv'
# pyWGCNA_data = PyWGCNA.WGCNA(species='homo sapiens', 
#                               geneExpPath=geneExp,
#                               TPMcutoff=0.00,
#                               outputPath='',
#                               save=True)

# print(pyWGCNA_data.geneExpr.to_df().head(5))

# # ---------------------------------------------

# # Preprocessing workflow:
# pyWGCNA_data.preprocess()

# # ---------------------------------------------

# # pyWGCNA_data.findModules()

# ==============================================
# '/Users/surya/Documents/GITHUB-REPOSITORIES/spatial_proteomics/2. feature_extraction/5xFAD_paper/Patient_20751_relativeAreaWise_cell_expression_count.csv'


tpm = pd.read_csv('/Users/surya/Documents/GITHUB-REPOSITORIES/spatial_proteomics/2. feature_extraction/5xFAD_paper/sampleInfo.csv', index_col=0)

test = PyWGCNA.WGCNA(species='Humans',
                     TPMcutoff=0.00,
                     geneExp=tpm,
                     TOMType='signed',
                     # powers=[1],
                     networkType="unsigned",
                     save=True)

test.preprocess()

test.findModules()

